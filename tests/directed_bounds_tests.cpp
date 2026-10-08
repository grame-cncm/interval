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

using namespace itv;

namespace {

// A 256-bit directed oracle is independent of the implementation's 53-bit
// evaluation. Checking both sides of its enclosure avoids trusting one rounded
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

    // Changing the host rounding mode must neither invalidate the MPFR bounds
    // nor be an observable side effect of an interval operation.
    const int saved = std::fegetround();
    for (int rounding : {FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO, FE_TONEAREST}) {
        std::fesetround(rounding);
        check("double: inclusion and host rounding mode are preserved", true,
              checkBinary(algebra, ref, &interval_algebra::Div, mpfr_div, 1, 10) &&
              std::fegetround() == rounding);
    }
    std::fesetround(saved);

    const mpfr_exp_t savedMin = mpfr_get_emin(), savedMax = mpfr_get_emax();
    mpfr_set_emin(-50);
    mpfr_set_emax(50);
    const interval restricted = algebra.Div(interval(tiny, tiny, -24), interval(2, 2, -24));
    check("double: a restricted caller MPFR range gives a safe fallback", true,
          restricted.lo() == -HUGE_VAL && restricted.hi() == HUGE_VAL &&
          mpfr_get_emin() == -50 && mpfr_get_emax() == 50);
    mpfr_set_emin(savedMin);
    mpfr_set_emax(savedMax);
    libmCompensation() = true;
    return reportCheckResults();
}
