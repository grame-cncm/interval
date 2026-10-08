#include <algorithm>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <limits>
#include <mpfr.h>
#include <random>

#include "interval/check.hh"
#include "interval/interval_algebra.hh"
#include "interval/affine_ops.hh"
#include "interval/directed_rounding.hh"

using namespace itv;

namespace {

// A 256-bit directed oracle is independent of the native kernel's exact dyadic
// comparisons and bounded series. Checking both sides of its enclosure avoids trusting one rounded
// high-precision approximation, especially for enormous or subnormal results.
class Reference {
   public:
    mpfr_t a, b, lo, hi;
    Reference() { mpfr_inits2(256, a, b, lo, hi, (mpfr_ptr)nullptr); }
    ~Reference() { mpfr_clears(a, b, lo, hi, (mpfr_ptr)nullptr); }
    Reference(const Reference&) = delete;
    Reference& operator=(const Reference&) = delete;

    bool enclosed(const interval& result) const
    {
        return !result.isEmpty() && mpfr_cmp_d(lo, result.lo()) >= 0 &&
               mpfr_cmp_d(hi, result.hi()) <= 0;
    }
};

using BinaryReference = int (*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
using BinaryMethod = interval (interval_algebra::*)(const interval&, const interval&) const;
using UnaryReference = int (*)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
using UnaryMethod = interval (interval_algebra::*)(const interval&) const;

bool checkBinary(const interval_algebra& algebra, Reference& ref, BinaryMethod method,
                 BinaryReference operation, double x, double y)
{
    mpfr_set_d(ref.a, x, MPFR_RNDN);
    mpfr_set_d(ref.b, y, MPFR_RNDN);
    operation(ref.lo, ref.a, ref.b, MPFR_RNDD);
    operation(ref.hi, ref.a, ref.b, MPFR_RNDU);
    return ref.enclosed((algebra.*method)(interval(x, x, -24), interval(y, y, -24)));
}

bool checkUnary(const interval_algebra& algebra, Reference& ref, UnaryMethod method,
                UnaryReference operation, double x)
{
    mpfr_set_d(ref.a, x, MPFR_RNDN);
    operation(ref.lo, ref.a, MPFR_RNDD);
    operation(ref.hi, ref.a, MPFR_RNDU);
    return ref.enclosed((algebra.*method)(interval(x, x, -24)));
}

// Finite binary64 bit patterns exercise every exponent, unlike a uniform real
// distribution which almost never reaches cancellation or gradual underflow.
double finiteSample(std::mt19937_64& random)
{
    uint64_t bits;
    do { bits = random(); } while ((bits & UINT64_C(0x7ff0000000000000)) == UINT64_C(0x7ff0000000000000));
    return std::bit_cast<double>(bits);
}

}  // namespace

int main()
{
    programPrecision() = 2;
    libmCompensation() = false; // test the analyzer's reference, not its target margin
    interval_algebra algebra;
    const double tiny = std::numeric_limits<double>::denorm_min();
    const double largest = std::numeric_limits<double>::max();
    const double nextOne = std::nextafter(1.0, HUGE_VAL);

    checkExact("double: half-ulp addition is enclosed, not collapsed",
               algebra.Add(interval(1, 1, -24), interval(0x1p-54)),
               interval(1, nextOne, -54));
    checkExact("double: subtraction encloses a negative half-ulp",
               algebra.Sub(interval(-1, -1, -24), interval(0x1p-54)),
               interval(-nextOne, -1, -54));
    check("double: a literal stays a singleton", true, algebra.FloatNum(0.1).is(0.1));
    check("double: an integral-looking float literal does not wrap as int32", true,
          algebra.Mul(algebra.FloatNum(65536), algebra.FloatNum(65536)).is(4294967296.0));
    check("double: an exact computed point stays a singleton", true,
          algebra.Mul(interval(0.5), interval(0.25)).is(0.125));
    check("double: nonpoint times zero remains exactly zero", true,
          algebra.Mul(interval(1, 2), interval(0)).isZero());
    check("double: underflow of a nonpoint does not become a false singleton", true,
          algebra.Mul(interval(tiny, 2 * tiny, -1074), interval(0.125, 0.25))
              == interval(0, tiny, -1098));
    const interval overflow = algebra.Mul(interval(largest, largest, -24), interval(2, 2, -24));
    check("double: outward overflow keeps max-finite below +inf", true,
          overflow.lo() == largest && overflow.hi() == HUGE_VAL);
    check("double: integral-looking huge bounds never enter an int cast", true,
          algebra.Add(interval(largest), interval(largest)).has(HUGE_VAL));
    check("double: an indeterminate infinite corner does not erase numeric values", true,
          algebra.Add(interval(HUGE_VAL), interval(-HUGE_VAL, 0)).has(HUGE_VAL));
    const interval quotient = algebra.Div(interval(1, 1, -24), interval(10, 10, -24));
    check("double: 1/10 is bracketed by adjacent doubles", true,
          quotient.lo() == std::nextafter(0.1, 0.0) && quotient.hi() == 0.1);
    const interval squareRoot = algebra.Sqrt(interval(2, 2, -24));
    check("double: irrational sqrt stays an enclosure", true,
          squareRoot.lo() < squareRoot.hi() && squareRoot.has(std::sqrt(2.0)));
    check("double: exact sqrt stays a singleton", true, algebra.Sqrt(interval(4)).is(2));
    check("double: sin evaluated at a point is not an emitted constant", true,
          !algebra.SinBounds(interval(0.5)).isconst());
    const double a = 1 + 0x1p-27, b = 1 - 0x1p-27;
    const interval contracted = algebra.Add(
        algebra.Mul(algebra.FloatNum(a), algebra.FloatNum(b)), algebra.FloatNum(-1));
    check("double: directed composition encloses the FMA cancellation result", true,
          contracted.has(std::fma(a, b, -1)) && contracted.lo() < 0);

    // Compare all elementary endpoints against a certified wider enclosure.
    Reference ref;
    std::mt19937_64 random(0x1a7e4ba1);
    for (const auto& entry : {
             std::pair<BinaryMethod, BinaryReference>{&interval_algebra::Add, mpfr_add},
             {&interval_algebra::Sub, mpfr_sub}, {&interval_algebra::Mul, mpfr_mul},
             {&interval_algebra::Div, mpfr_div}}) {
        bool enclosed = true;
        for (int i = 0; i < 4000 && enclosed; ++i) {
            const double x = finiteSample(random), y = finiteSample(random);
            if (y == 0 && entry.first == &interval_algebra::Div) continue;
            enclosed = checkBinary(algebra, ref, entry.first, entry.second, x, y);
        }
        check("double: elementary operation encloses 4000 MPFR references", true, enclosed);
    }

    const struct UnaryCase { const char* name; UnaryMethod method; UnaryReference reference; } cases[] = {
        {"sqrt", &interval_algebra::Sqrt, mpfr_sqrt},
        {"acos", &interval_algebra::AcosBounds, mpfr_acos},
        {"acosh", &interval_algebra::AcoshBounds, mpfr_acosh},
        {"asin", &interval_algebra::AsinBounds, mpfr_asin},
        {"asinh", &interval_algebra::AsinhBounds, mpfr_asinh},
        {"atan", &interval_algebra::AtanBounds, mpfr_atan},
        {"atanh", &interval_algebra::AtanhBounds, mpfr_atanh},
        {"cos", &interval_algebra::CosBounds, mpfr_cos},
        {"cosh", &interval_algebra::CoshBounds, mpfr_cosh},
        {"exp", &interval_algebra::ExpBounds, mpfr_exp},
        {"log", &interval_algebra::LogBounds, mpfr_log},
        {"log10", &interval_algebra::Log10Bounds, mpfr_log10},
        {"sin", &interval_algebra::SinBounds, mpfr_sin},
        {"sinh", &interval_algebra::SinhBounds, mpfr_sinh},
        {"tan", &interval_algebra::TanBounds, mpfr_tan},
        {"tanh", &interval_algebra::TanhBounds, mpfr_tanh}
    };
    for (const auto& entry : cases) {
        bool enclosed = true;
        for (double x : {-largest, -1e100, -10.0, -1.0, -0.5, -tiny, 0.0,
                         tiny, 0.5, 1.0, 2.0, 10.0, 1e100, largest}) {
            mpfr_set_d(ref.a, x, MPFR_RNDN);
            entry.reference(ref.lo, ref.a, MPFR_RNDD);
            if (mpfr_nan_p(ref.lo)) continue; // check numeric values on valid domains
            enclosed = enclosed && checkUnary(algebra, ref, entry.method, entry.reference, x);
        }
        check(std::string("double: ") + entry.name + " encloses MPFR at boundary arguments", true, enclosed);
    }

    // Neighbours of reduction, pole, domain, overflow and subnormal thresholds
    // exercise regions that a random exponent distribution can easily miss.
    for (const auto& entry : cases) {
        bool enclosed = true;
        for (double base : {-1024.0, -745.0, -744.0, -710.0, -708.0, -32.0,
                            -1.0, -0.5, -tiny, 0.0, tiny, 0.5, 1.0,
                            0x1.921fb54442d18p+0, 0x1.921fb54442d18p+1,
                            32.0, 708.0, 709.0, 710.0, 1024.0, 0x1p20}) {
            for (double x : {std::nextafter(base, -HUGE_VAL), base,
                             std::nextafter(base, HUGE_VAL)}) {
                mpfr_set_d(ref.a, x, MPFR_RNDN);
                entry.reference(ref.lo, ref.a, MPFR_RNDD);
                if (mpfr_nan_p(ref.lo)) continue;
                enclosed = enclosed && checkUnary(algebra, ref, entry.method, entry.reference, x);
            }
        }
        check(std::string("native ") + entry.name + " encloses threshold neighbours", true, enclosed);
    }

    // Direct kernel checks cover modulo (including enormous quotients), sqrt,
    // exactness and exceptional scales independently of interval image helpers.
    using namespace itv::detail;
    for (const auto& entry : {
             std::pair<BinaryOp, BinaryReference>{BinaryOp::Mul, mpfr_mul},
             {BinaryOp::Div, mpfr_div}, {BinaryOp::Fmod, mpfr_fmod},
             {BinaryOp::Remainder, mpfr_remainder}, {BinaryOp::Atan2, mpfr_atan2}}) {
        bool enclosed = true;
        for (int i = 0; i < 4000 && enclosed; ++i) {
            const double x = finiteSample(random), y = finiteSample(random);
            if (y == 0) continue;
            mpfr_set_d(ref.a, x, MPFR_RNDN); mpfr_set_d(ref.b, y, MPFR_RNDN);
            entry.second(ref.lo, ref.a, ref.b, MPFR_RNDD);
            entry.second(ref.hi, ref.a, ref.b, MPFR_RNDU);
            enclosed = ref.enclosed(interval(directedBinary(entry.first, x, y, Direction::Down),
                                             directedBinary(entry.first, x, y, Direction::Up)));
            if (!enclosed) std::cerr << "kernel failure op=" << int(entry.first)
                                    << " x=" << std::hexfloat << x << " y=" << y << '\n';
        }
        check("native binary kernel encloses 4000 MPFR references", true, enclosed);
    }
    for (const auto& entry : cases) {
        bool enclosed = true;
        for (int i = 0; i < 1000 && enclosed; ++i) {
            // Mix all binary exponents with ordinary arguments. The latter
            // exercise the series and reduction, rather than only broad hulls.
            const double x = i % 2 ? finiteSample(random)
                                   : std::uniform_real_distribution<double>(-32, 32)(random);
            mpfr_set_d(ref.a, x, MPFR_RNDN);
            entry.reference(ref.lo, ref.a, MPFR_RNDD);
            if (mpfr_nan_p(ref.lo)) continue;
            enclosed = checkUnary(algebra, ref, entry.method, entry.reference, x);
            if (!enclosed) std::cerr << "unary failure " << entry.name
                                    << " x=" << std::hexfloat << x << '\n';
        }
        check(std::string("native ") + entry.name + " encloses random MPFR references", true, enclosed);
    }
    bool scalesEnclosed = true, powersEnclosed = true, piEnclosed = true;
    for (int i = 0; i < 4000 && scalesEnclosed; ++i) {
        const double x = finiteSample(random);
        const int exponent = std::uniform_int_distribution<int>(-1500, 1500)(random);
        mpfr_set_d(ref.a, x, MPFR_RNDN);
        mpfr_mul_2si(ref.lo, ref.a, exponent, MPFR_RNDD);
        mpfr_mul_2si(ref.hi, ref.a, exponent, MPFR_RNDU);
        scalesEnclosed = ref.enclosed(interval(directedScale(x, exponent, Direction::Down),
                                               directedScale(x, exponent, Direction::Up)));
    }
    check("native power-of-two scaling encloses 4000 MPFR references", true, scalesEnclosed);
    for (int i = 0; i < 1000 && powersEnclosed; ++i) {
        const double x = std::abs(finiteSample(random));
        const double y = i % 2 ? std::uniform_real_distribution<double>(-16, 16)(random)
                               : std::uniform_int_distribution<int>(-31, 31)(random);
        powersEnclosed = checkBinary(algebra, ref, &interval_algebra::PowBounds, mpfr_pow, x, y);
        if (!powersEnclosed) std::cerr << "power failure x=" << std::hexfloat << x << " y=" << y << '\n';
    }
    check("native powers enclose random MPFR references", true, powersEnclosed);
    mpfr_const_pi(ref.lo, MPFR_RNDD); mpfr_const_pi(ref.hi, MPFR_RNDU);
    piEnclosed = ref.enclosed(interval(directedPi(Direction::Down), directedPi(Direction::Up)));
    check("native pi constants enclose MPFR pi", true, piEnclosed);
    mpfr_set_ui(ref.a, 2, MPFR_RNDN);
    mpfr_log(ref.lo, ref.a, MPFR_RNDD); mpfr_log(ref.hi, ref.a, MPFR_RNDU);
    check("native ln2 constants enclose MPFR log(2)", true,
          ref.enclosed(interval(0x1.62e42fefa39efp-1, 0x1.62e42fefa39f0p-1)));
    mpfr_set_ui(ref.a, 10, MPFR_RNDN);
    mpfr_log(ref.lo, ref.a, MPFR_RNDD); mpfr_log(ref.hi, ref.a, MPFR_RNDU);
    check("native ln10 constants enclose MPFR log(10)", true,
          ref.enclosed(interval(0x1.26bb1bbb55515p+1, 0x1.26bb1bbb55516p+1)));


    // Check transfer functions on real corridors as well as points, including
    // interior trigonometric extrema. Sampling validates implementation paths;
    // the analytic argument remains necessary for the inclusion guarantee.
    for (const auto& entry : cases) {
        bool enclosed = true;
        for (int i = 0; i < 200 && enclosed; ++i) {
            double lo = std::uniform_real_distribution<double>(-16, 16)(random);
            double hi = std::uniform_real_distribution<double>(-16, 16)(random);
            if (lo > hi) std::swap(lo, hi);
            const interval result = (algebra.*entry.method)(interval(lo, hi, -24));
            for (int sample = 0; sample <= 8 && enclosed; ++sample) {
                // Rounded interpolation can step outside the input at an endpoint.
                const double x = std::clamp(lo + (hi - lo) * (sample / 8.0), lo, hi);
                mpfr_set_d(ref.a, x, MPFR_RNDN);
                entry.reference(ref.lo, ref.a, MPFR_RNDD);
                entry.reference(ref.hi, ref.a, MPFR_RNDU);
                if (mpfr_nan_p(ref.lo)) continue;
                enclosed = ref.enclosed(result);
                if (!enclosed) std::cerr << "corridor failure " << entry.name
                                        << " lo=" << lo << " hi=" << hi << " sample=" << x << '\n';
            }
        }
        check(std::string("native ") + entry.name + " image encloses MPFR samples", true, enclosed);
    }

    // Periodic range decisions must use pi enclosures, not rounded fmod phases.
    const double pi = std::acos(-1.0), halfPi = pi / 2;
    check("double: sin includes an interior maximum near pi/2", true,
          algebra.SinBounds(interval(std::nextafter(halfPi, 0.0),
                                     std::nextafter(halfPi, HUGE_VAL))).hi() == 1);
    check("double: cos includes an interior minimum near pi", true,
          algebra.CosBounds(interval(std::nextafter(pi, 0.0),
                                     std::nextafter(pi, HUGE_VAL))).lo() == -1);
    check("double: tan includes a pole near pi/2", true,
          algebra.TanBounds(interval(std::nextafter(halfPi, 0.0),
                                     std::nextafter(halfPi, HUGE_VAL))).isUnbounded());
    check("double: periodic extrema survive huge argument uncertainty", true,
          algebra.SinBounds(interval(1e100, std::nextafter(1e100, HUGE_VAL))) == interval(-1, 1));
    check("double: atan2 encloses both sides of its negative-axis cut", true,
          algebra.Atan2Bounds(interval(-tiny, tiny), interval(-2, -1)).lo() < -pi &&
          algebra.Atan2Bounds(interval(-tiny, tiny), interval(-2, -1)).hi() > pi);
    check("double: negative power at zero covers infinity", true,
          algebra.PowBounds(interval(0, 1), interval(-0.5)).has(HUGE_VAL));
    check("double: infinite negative bases keep exceptional pow values", true,
          algebra.PowBounds(interval(-HUGE_VAL), interval(0.5)).has(HUGE_VAL) &&
          algebra.PowBounds(interval(-HUGE_VAL), interval(-0.5)).has(0));
    check("double: a negative reciprocal corridor keeps the included +zero", true,
          algebra.Inv(interval(-1, 0)).has(HUGE_VAL));
    check("double: negative-base negative integer powers keep their sign", true,
          algebra.PowBounds(interval(-2, -1), interval(-3, -1, -24)).has(-1) &&
          algebra.PowBounds(interval(-2, -1), interval(-3, -1, -24)).has(1));
    check("double: fmod of integer-looking operands has floating semantics", true,
          algebra.Fmod(interval(3), interval(2)).is(1));
    check("double: fmod cannot lose a tiny result in a huge quotient", true,
          algebra.Fmod(interval(1e100, 1e101), interval(tiny, 2 * tiny)).has(tiny));
    check("double: remainder half-divisor underflow rounds outward", true,
          algebra.Remainder(interval(-tiny, tiny), interval(tiny, tiny, -24)).has(tiny));

    // Chord coefficients and their evaluation can themselves round inward.
    // Check an exact rational t/T against both directed affine bounds.
    double lower0, lower1, upper0, upper1;
    achord(0, 1, 10, lower0, lower1);
    achord(0, 1, 10, upper0, upper1, true);
    AffItv chord{lower0, lower1, upper0, upper1, -24};
    bool chordContains = true;
    for (int t = 0; t <= 10; ++t) {
        mpfr_set_si(ref.a, t, MPFR_RNDN);
        mpfr_div_ui(ref.lo, ref.a, 10, MPFR_RNDD);
        mpfr_div_ui(ref.hi, ref.a, 10, MPFR_RNDU);
        chordContains = chordContains && ref.enclosed(interval(chord.lo(t), chord.hi(t)));
    }
    check("double: affine chord encloses exact rational interpolation", true, chordContains);
    affine_algebra affine(10);
    const AffItv affineSum = affine.Add(fromItv(interval(1, 1, -24)),
                                        fromItv(interval(0x1p-54)));
    check("double: affine addition does not collapse a half-ulp", true,
          affineSum.a0 == 1 && affineSum.b0 == nextOne);
    const AffItv infiniteJoin = ajoin({-HUGE_VAL, 0, 0, 0, -24},
                                      {0, 0, HUGE_VAL, 0, -24}, 10);
    check("double: affine join keeps both unbounded sides", true,
          infiniteJoin.lo(5) == -HUGE_VAL && infiniteJoin.hi(5) == HUGE_VAL);
    check("double: affine inclusion does not certify an absorbed rate", false,
          aleq({1, 0x1p-54, 1, 0x1p-54, -24}, {1, 0, 1, 0, -24}, 1));
    check("double: directed affine inclusion remains reflexive", true, aleq(chord, chord, 10));
    check("double: indeterminate affine corners keep numeric paths", true,
          toItv(affine.Add(fromItv(interval(HUGE_VAL)),
                          fromItv(interval(-HUGE_VAL, 0))), 10).has(HUGE_VAL));

    // Changing the host rounding mode must neither invalidate the native bounds
    // nor be an observable side effect of an interval operation.
    const int saved = std::fegetround();
    for (int rounding : {FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO, FE_TONEAREST}) {
        std::fesetround(rounding);
        check("double: inclusion and host rounding mode are preserved", true,
              checkBinary(algebra, ref, &interval_algebra::Div, mpfr_div, 1, 10) &&
              std::fegetround() == rounding);
        bool enclosed = true;
        for (int i = 0; i < 100 && enclosed; ++i) {
            const double x = finiteSample(random), y = finiteSample(random);
            enclosed = checkBinary(algebra, ref, &interval_algebra::Add, mpfr_add, x, y) &&
                       checkBinary(algebra, ref, &interval_algebra::Sub, mpfr_sub, x, y) &&
                       checkBinary(algebra, ref, &interval_algebra::Mul, mpfr_mul, x, y) &&
                       (y == 0 || checkBinary(algebra, ref, &interval_algebra::Div, mpfr_div, x, y)) &&
                       checkUnary(algebra, ref, &interval_algebra::Sqrt, mpfr_sqrt, std::abs(x));
        }
        check("native elementary bounds enclose MPFR in each CPU rounding mode", true, enclosed);

    }
    std::fesetround(saved);

    libmCompensation() = true;
    return reportCheckResults();
}
