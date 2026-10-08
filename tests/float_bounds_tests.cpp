#include <algorithm>
#include <array>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>

#ifdef INTERVAL_FLOAT_MPFR_ORACLE
#include <mpfr.h>
#endif

#include "interval/check.hh"
#include "interval/interval_algebra.hh"
#include "interval/affine_ops.hh"

using namespace itv;

namespace {

using BinaryMethod = interval (interval_algebra::*)(const interval&, const interval&) const;
using UnaryMethod = interval (interval_algebra::*)(const interval&) const;

// Negative LSB is the library's floating-nature marker. These domains contain
// exact target values; one in four samples also checks a nonpoint image.
interval domain(float x, bool wider = false)
{
    return {double(x), wider ? double(std::nextafter(x, INFINITY)) : double(x), -149};
}

// Sample signs, significands and every binary32 exponent, including subnormals.
// Uniform real distributions almost never reach underflow or huge cancellation.
float sample(std::mt19937_64& random)
{
    uint32_t bits;
    do { bits = uint32_t(random()); } while ((bits & 0x7f800000u) == 0x7f800000u);
    return std::bit_cast<float>(bits);
}

struct BinaryCase {
    const char* name;
    BinaryMethod method;
    float (*execute)(float, float);
};

const BinaryCase binaryCases[] = {
    {"add", &interval_algebra::Add, [](float a, float b) { return a + b; }},
    {"sub", &interval_algebra::Sub, [](float a, float b) { return a - b; }},
    {"mul", &interval_algebra::Mul, [](float a, float b) { return a * b; }},
    {"div", &interval_algebra::Div, [](float a, float b) { return a / b; }},
    {"fmod", &interval_algebra::Fmod, [](float a, float b) { return std::fmod(a, b); }},
    {"remainder", &interval_algebra::Remainder, [](float a, float b) { return std::remainder(a, b); }},
    {"pow", &interval_algebra::Pow, [](float a, float b) { return std::pow(a, b); }},
    {"atan2", &interval_algebra::Atan2, [](float a, float b) { return std::atan2(a, b); }}
};

struct UnaryCase {
    const char* name;
    UnaryMethod referenceBounds;
    UnaryMethod targetBounds;
    float (*execute)(float);
};

const UnaryCase unaryCases[] = {
    {"sqrt", &interval_algebra::Sqrt, &interval_algebra::Sqrt, [](float x) { return std::sqrt(x); }},
    {"acos", &interval_algebra::AcosBounds, &interval_algebra::Acos, [](float x) { return std::acos(x); }},
    {"acosh", &interval_algebra::AcoshBounds, &interval_algebra::Acosh, [](float x) { return std::acosh(x); }},
    {"asin", &interval_algebra::AsinBounds, &interval_algebra::Asin, [](float x) { return std::asin(x); }},
    {"asinh", &interval_algebra::AsinhBounds, &interval_algebra::Asinh, [](float x) { return std::asinh(x); }},
    {"atan", &interval_algebra::AtanBounds, &interval_algebra::Atan, [](float x) { return std::atan(x); }},
    {"atanh", &interval_algebra::AtanhBounds, &interval_algebra::Atanh, [](float x) { return std::atanh(x); }},
    {"cos", &interval_algebra::CosBounds, &interval_algebra::Cos, [](float x) { return std::cos(x); }},
    {"cosh", &interval_algebra::CoshBounds, &interval_algebra::Cosh, [](float x) { return std::cosh(x); }},
    {"exp", &interval_algebra::ExpBounds, &interval_algebra::Exp, [](float x) { return std::exp(x); }},
    {"log", &interval_algebra::LogBounds, &interval_algebra::Log, [](float x) { return std::log(x); }},
    {"log10", &interval_algebra::Log10Bounds, &interval_algebra::Log10, [](float x) { return std::log10(x); }},
    {"sin", &interval_algebra::SinBounds, &interval_algebra::Sin, [](float x) { return std::sin(x); }},
    {"sinh", &interval_algebra::SinhBounds, &interval_algebra::Sinh, [](float x) { return std::sinh(x); }},
    {"tan", &interval_algebra::TanBounds, &interval_algebra::Tan, [](float x) { return std::tan(x); }},
    {"tanh", &interval_algebra::TanhBounds, &interval_algebra::Tanh, [](float x) { return std::tanh(x); }}
};

#ifdef INTERVAL_FLOAT_MPFR_ORACLE
using BinaryReference = int (*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
using UnaryReference = int (*)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
const BinaryReference binaryReferences[] = {
    mpfr_add, mpfr_sub, mpfr_mul, mpfr_div, mpfr_fmod, mpfr_remainder, mpfr_pow, mpfr_atan2};
const UnaryReference unaryReferences[] = {
    mpfr_sqrt, mpfr_acos, mpfr_acosh, mpfr_asin, mpfr_asinh, mpfr_atan, mpfr_atanh,
    mpfr_cos, mpfr_cosh, mpfr_exp, mpfr_log, mpfr_log10, mpfr_sin, mpfr_sinh, mpfr_tan, mpfr_tanh};

// The reference enclosure and target rounding are different checks. Precision
// 256 with unrestricted exponents checks the real image; precision 24 with the
// binary32 exponent range and subnormalize checks gradual underflow and overflow.
class Oracle {
  public:
    mpfr_t a, b, c, lo, hi, rounded;
    Oracle()
    {
        mpfr_inits2(256, a, b, c, lo, hi, (mpfr_ptr)nullptr);
        mpfr_init2(rounded, 24);
    }
    ~Oracle() { mpfr_clears(a, b, c, lo, hi, rounded, (mpfr_ptr)nullptr); }
    Oracle(const Oracle&) = delete;
    Oracle& operator=(const Oracle&) = delete;

    bool enclosed(const interval& result) const
    {
        return mpfr_nan_p(lo) || (!result.isEmpty() &&
            mpfr_cmp_d(lo, result.lo()) >= 0 && mpfr_cmp_d(hi, result.hi()) <= 0);
    }
    // Inputs are exact binary32. Install limits only around the target operation,
    // and restore them before doing the wide reference or another library call.
    template <typename Operation> double target(Operation operation)
    {
        const mpfr_exp_t emin = mpfr_get_emin(), emax = mpfr_get_emax();
        mpfr_set_emin(-148);
        mpfr_set_emax(128);
        const int inexact = operation();
        mpfr_subnormalize(rounded, inexact, MPFR_RNDN);
        const double value = mpfr_get_d(rounded, MPFR_RNDN);
        mpfr_set_emin(emin);
        mpfr_set_emax(emax);
        return value;
    }
};
#endif

// Print the first failure's exact operands, not thousands of successful samples.
bool contains(const char* name, const interval& result, double value, float x, float y = 0)
{
    if (std::isnan(value) || result.has(value)) return true;
    std::cerr << name << " x=" << std::hexfloat << x << " y=" << y
              << " value=" << value << " interval=" << result << '\n';
    return false;
}

} // namespace

int main()
{
    static_assert(std::numeric_limits<float>::is_iec559 && std::numeric_limits<float>::digits == 24);
    programPrecision() = 1;
    libmCompensation() = true;
    const interval_algebra algebra;
    const double tiny = std::numeric_limits<float>::denorm_min();
    const double largest = std::numeric_limits<float>::max();

    check("float: exact literals remain points", true, algebra.FloatNum(16777217).is(16777216));
    check("float: constants are rounded once, domains outward", true,
          algebra.FloatNum(0.1).is(0.1f) && interval(0.1, 0.1, -24).has(0.1));
    check("float: exact operations retain point propagation", true,
          algebra.Mul(algebra.FloatNum(0.5), algebra.FloatNum(0.25)).is(0.125));
    check("float: integer injection remains exact", true, algebra.IntNum(16777217).is(16777217));
    check("float: positive underflow encloses real value and rounded zero", true,
          interval(tiny / 4, tiny / 4, -149).lo() == 0 &&
          interval(tiny / 4, tiny / 4, -149).hi() == tiny);
    check("float: negative underflow encloses real value and rounded zero", true,
          interval(-tiny / 4, -tiny / 4, -149).lo() == -tiny &&
          interval(-tiny / 4, -tiny / 4, -149).hi() == 0);
    check("float: overflow brackets max-finite and infinity", true,
          interval(2 * largest, 2 * largest, -24).lo() == largest &&
          interval(2 * largest, 2 * largest, -24).hi() == HUGE_VAL);
    check("float: literals below the overflow midpoint round to max-finite", true,
          algebra.FloatNum(std::nextafter(largest, HUGE_VAL)).is(largest) &&
          algebra.FloatNum(0x1.ffffffp127).is(HUGE_VAL));
    check("float: a positive estimated LSB cannot turn multiplication into int32", true,
          algebra.Mul(interval(0.5, 0.5, -1), interval(0x1p30, 0x1p30, 30)).lsb() < 0);
    const interval zero = algebra.Sub(algebra.FloatNum(1), algebra.FloatNum(1));
    check("float: computed zero retains floating nature for mixed conversion", true,
          zero.lsb() < 0 && algebra.Add(zero, algebra.IntNum(16777217)).is(16777216));
    check("float: integral-valued libm results retain floating nature", true,
          algebra.Add(algebra.Floor(algebra.FloatNum(1e9)),
                      algebra.Floor(algebra.FloatNum(2e9))).has(3e9f));
    check("float: exp10 with an integer exponent does not wrap as int32", true,
          algebra.Exp10(algebra.IntNum(10)).has(1e10f) &&
          algebra.Exp10(algebra.IntNum(-1)).has(0.1f));
    check("float: explicit zero remains an exact float point", true,
          algebra.FloatNum(0).isZero() && algebra.FloatNum(0).lsb() < 0);
    check("float: merged signed zero includes both reciprocal infinities", true,
          algebra.Inv(domain(-0.0f)).has(-HUGE_VAL) && algebra.Inv(domain(0.0f)).has(HUGE_VAL));
    check("float: negative power of signed zero includes negative infinity", true,
          algebra.Pow(domain(-0.0f), algebra.FloatNum(-1)).has(-HUGE_VAL));
    const interval slider = algebra.HSlider(interval(0), algebra.FloatNum(0),
        algebra.FloatNum(0), algebra.FloatNum(16777216), algebra.IntNum(1));
    check("float: a slider with integer step still has floating nature", true, slider.lsb() < 0);
    check("precision estimate: coincident function evaluations use a defined fallback", true,
          exactPrecisionUnary(std::sin, 1e20L, 0x1p-24L) == INT_MIN);

    // The affine bridge must use the same float rules. Endpoint chords through
    // rounded values alone cannot prove that intermediate staircase values fit.
    const affine_algebra affine(128);
    // The initial condition of a delay is a zero of the signal's nature : an integer
    // stays an integer (and keeps its int32 wrap) once delayed.
    const AffItv delayedInt = affine.Delay(fromItv(interval(0, 1000, 0)), affine.IntNum(1));
    check("float: an integer delayed in the affine domain stays an integer", true, delayedInt.lsb >= 0);
    const AffItv moving{4097, 0x1p-7, 4097, 0x1p-7, -24};
    const AffItv product = affine.Mul(moving, affine.FloatNum(4097));
    bool affineIncluded = true;
    for (int t = 0; t <= 128; ++t) {
        volatile float a = float(4097 + t * 0x1p-7);
        volatile float b = 4097;
        affineIncluded = affineIncluded && toItv(product, 128).has(a * b);
    }
    check("float: affine product encloses intermediate rounded staircase values", true, affineIncluded);
    const AffItv delay = affine.IntCast(affine.Sub(
        affine.Sub(product, affine.FloatNum(16785408)), affine.FloatNum(1)));
    check("float: affine delay also includes the negative counterexample", true, toItv(delay, 128).has(-1));
    check("float: affine constant casts round integer values", true,
          toItv(affine.FloatCast(affine.IntNum(16777217)), 128).is(16777216));

#ifndef __EMSCRIPTEN__
    // Directed narrowing corrects either adjacent conversion result. Changing
    // the host rounding mode must not make a mathematical endpoint disappear.
    const int savedMode = std::fegetround();
    bool directedModes = true;
    for (int mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO}) {
        if (std::fesetround(mode) != 0) continue;
        for (double x : {1 + 0x1p-24, -1 - 0x1p-24, tiny / 4, -tiny / 4, 2 * largest, -2 * largest})
            directedModes = directedModes && interval(x, x, -24).has(x);
        directedModes = directedModes && algebra.FloatNum(16777217).is(16777216)
                        && algebra.FloatNum(0.1).is(0.1f);
        directedModes = directedModes && std::fegetround() == mode;
    }
    std::fesetround(savedMode);
    check("float: directed narrowing encloses bounds and preserves host rounding mode", true, directedModes);
#endif

    // Integer-to-float conversion must occur before a mixed operation. Test both
    // signs and ties, rather than only a FloatCast used in isolation.
    for (int x : {16777217, -16777217, INT_MAX, INT_MIN + 1}) {
        for (const auto& entry : binaryCases) {
            const interval result = (algebra.*entry.method)(algebra.IntNum(x), algebra.FloatNum(-16777216));
            check(std::string("float: mixed conversion before ") + entry.name, true,
                  contains(entry.name, result, entry.execute(float(x), -16777216.0f), float(x)));
        }
    }

    std::mt19937_64 random(0xf10a7b0u);
#ifdef INTERVAL_FLOAT_MPFR_ORACLE
    Oracle oracle;
#endif
    for (size_t index = 0; index < std::size(binaryCases); ++index) {
        const auto& entry = binaryCases[index];
        bool included = true;
        const int count = index < 6 ? 4000 : 600;
        for (int i = 0; i < count && included; ++i) {
            const float x = sample(random), y = sample(random);
            const interval result = (algebra.*entry.method)(domain(x, i % 4 == 0), domain(y, i % 4 == 0));
            included = contains(entry.name, result, entry.execute(x, y), x, y);
#ifdef INTERVAL_FLOAT_MPFR_ORACLE
            mpfr_set_d(oracle.a, x, MPFR_RNDN); mpfr_set_d(oracle.b, y, MPFR_RNDN);
            const auto operation = binaryReferences[index];
            operation(oracle.lo, oracle.a, oracle.b, MPFR_RNDD);
            operation(oracle.hi, oracle.a, oracle.b, MPFR_RNDU);
            const double target = oracle.target([&] { return operation(oracle.rounded, oracle.a, oracle.b, MPFR_RNDN); });
            included = included && oracle.enclosed(result) && contains(entry.name, result, target, x, y);
#endif
        }
        check(std::string("float: sampled point/nonpoint ") + entry.name + " inclusion", true, included);
    }

    // These thresholds supplement the exponent distribution: libm domain edges,
    // argument reduction, poles, and binary32 exp/sinh overflow/underflow.
    const float thresholds[] = {-INFINITY, -float(largest), -1e20f, -104, -103, -89, -88,
        -100, -1, -0.5, -float(tiny), 0, float(tiny), 0.5, 1, 2, 88, 89, 100, 1e20f, float(largest), INFINITY};
    for (size_t index = 0; index < std::size(unaryCases); ++index) {
        const auto& entry = unaryCases[index];
        bool included = true;
        auto test = [&](float x) {
            const interval input = domain(x);
            const interval result = (algebra.*entry.targetBounds)(input);
            if (!contains(entry.name, result, entry.execute(x), x)) return false;
#ifdef INTERVAL_FLOAT_MPFR_ORACLE
            mpfr_set_d(oracle.a, x, MPFR_RNDN);
            const auto operation = unaryReferences[index];
            operation(oracle.lo, oracle.a, MPFR_RNDD);
            operation(oracle.hi, oracle.a, MPFR_RNDU);
            const interval reference = (algebra.*entry.referenceBounds)(input);
            const double target = oracle.target([&] { return operation(oracle.rounded, oracle.a, MPFR_RNDN); });
            if (!oracle.enclosed(reference) || !contains(entry.name, reference, target, x)) return false;
#endif
            return true;
        };
        for (float base : thresholds) {
            for (float x : {std::nextafter(base, -INFINITY), base, std::nextafter(base, INFINITY)})
                included = included && test(x);
        }
        for (int i = 0; i < 400 && included; ++i) {
            const float x = i % 2 ? sample(random) : std::uniform_real_distribution<float>(-32, 32)(random);
            included = test(x);
        }
        check(std::string("float: ") + entry.name + " threshold/random inclusion", true, included);
    }

    // Contraction replaces two roundings by one. Outward composition must cover
    // that one rounding, including finite results after a product would overflow.
    bool contractions = true;
    for (int i = 0; i < 4000 && contractions; ++i) {
        const float a = sample(random), b = sample(random), c = sample(random);
        const interval result = algebra.Add(algebra.Mul(domain(a), domain(b)), domain(c));
        contractions = contains("fma", result, std::fma(a, b, c), a, b);
#ifdef INTERVAL_FLOAT_MPFR_ORACLE
        mpfr_set_d(oracle.a, a, MPFR_RNDN); mpfr_set_d(oracle.b, b, MPFR_RNDN); mpfr_set_d(oracle.c, c, MPFR_RNDN);
        const double target = oracle.target([&] {
            return mpfr_fma(oracle.rounded, oracle.a, oracle.b, oracle.c, MPFR_RNDN);
        });
        contractions = contractions && contains("mpfr fma", result, target, a, b);
#endif
    }
    check("float: 4000 contracted product/addition inclusions", true, contractions);
    return reportCheckResults();
}
