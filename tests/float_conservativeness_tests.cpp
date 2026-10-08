#include <array>
#include <cfenv>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string_view>

#ifdef INTERVAL_FLOAT_MPFR_ORACLE
#include <mpfr.h>
#endif

#include "interval/interval_algebra.hh"

namespace {

using itv::interval;
using itv::interval_algebra;

// Domains already contain exact binary32 values. Negative LSB explicitly marks
// floating inputs, so a missing value cannot be blamed on the integer overload.
interval floatDomain(float lo, float hi)
{
    return {double(lo), double(hi), -24};
}

// Each witness checks inclusion rather than freezing historical erroneous
// endpoints. Premises concern the inputs/target execution only; a wider sound
// result satisfies the regression just as a tighter sound result does.
struct Witness {
    std::string_view name;
    interval predicted;
    double observed;
    bool validInputs;
};

// The sum lies halfway between adjacent floats above 2^24. Nearest-even
// execution differs from retaining the exactly computed integer-valued bound.
Witness addition(const interval_algebra& algebra)
{
    volatile float x = 16777216.0f, y = 1.0f;
    volatile float result = x + y;
    const interval a = floatDomain(x, x), b = floatDomain(y, y);
    return {"addition", algebra.Add(a, b), double(result),
            a.has(x) && b.has(y) && result == 16777216.0f};
}

// A nonpoint domain prevents singleton handling from explaining the defect:
// the rounded lower-end product falls below the real-valued lower endpoint.
Witness multiplication(const interval_algebra& algebra)
{
    volatile float a = 4097.0f, b = 4097.0f;
    volatile float result = a * b;
    const interval domain = floatDomain(4097.0f, 4098.0f);
    return {"multiplication", algebra.Mul(domain, domain), double(result),
            domain.has(a) && domain.has(b) && result == 16785408.0f};
}

// Amplify that missed lower endpoint into a false nonnegative-delay proof.
// Only defined separate float operations and an in-range int cast are needed.
Witness delay(const interval_algebra& algebra)
{
    volatile float a = 4097.0f, b = 4097.0f;
    volatile float product = a * b;
    volatile float difference = product - 16785408.0f;
    volatile float result = difference - 1.0f;
    const interval domain = floatDomain(4097.0f, 4098.0f);
    const interval predicted = algebra.IntCast(algebra.Sub(
        algebra.Sub(algebra.Mul(domain, domain), floatDomain(16785408.0f, 16785408.0f)),
        floatDomain(1.0f, 1.0f)));
    // The three volatile stages enforce separate binary32 roundings. All values
    // are finite and the final -1 is representable as int: this cast is defined.
    return {"delay", predicted, double(int(result)),
            domain.has(a) && domain.has(b) && result == -1.0f};
}

// Explicit FloatCast must describe conversion, even when its input is an exact
// integer singleton; preserving that integer is not a valid float result.
Witness conversion(const interval_algebra& algebra)
{
    volatile int source = 16777217;
    volatile float result = float(source);
    const interval integer = algebra.IntNum(source);
    return {"conversion", algebra.FloatCast(integer), double(result),
            integer.has(source) && result == 16777216.0f};
}

// Preserve the small product residual with FMA, then scale it to an integer:
// an analyzer that follows only separate roundings loses this negative value.
Witness contraction(const interval_algebra& algebra)
{
    volatile float a = 1.0f + 0x1p-13f, b = 1.0f - 0x1p-13f;
    volatile float product = a * b;
    volatile float separate = product - 1.0f;
    // Explicit fma represents an allowed contraction of a*b-1, without asking
    // the compiler to reassociate anything or enabling fast-math for the tests.
    volatile float fused = std::fma(float(a), float(b), -1.0f);
    volatile float scaled = fused * 0x1p26f;
    const interval difference = algebra.Sub(
        algebra.Mul(floatDomain(a, a), floatDomain(b, b)), floatDomain(1.0f, 1.0f));
    const interval predicted = algebra.IntCast(
        algebra.Mul(difference, floatDomain(0x1p26f, 0x1p26f)));
    return {"contraction", predicted, double(int(scaled)),
            separate == 0 && fused == -0x1p-26f && scaled == -1 &&
            floatDomain(a, a).has(a) && floatDomain(b, b).has(b)};
}

// Keep both a modest argument and the historically broken huge argument: the
// former tests tight reduction and the latter tests the conservative fallback.
Witness sineAt(const interval_algebra& algebra, float input, std::string_view name)
{
    volatile float x = input;
    volatile float result = std::sin(float(x));
    const interval domain = floatDomain(x, x);
    // Test the public compensated function, with compensation enabled. A point
    // input is a runtime parameter domain, not evidence of constant folding.
    return {name, algebra.Sin(domain), double(result),
            domain.has(x) && std::isfinite(result) && std::abs(result) <= 1};
}

Witness sine(const interval_algebra& algebra) { return sineAt(algebra, 100.0f, "sine"); }

// The legacy path cast log2(0) to int at this magnitude. The native reference
// avoids that metadata path and returns a sound trig range for huge arguments.
Witness largeSine(const interval_algebra& algebra)
{
    return sineAt(algebra, 1.0e20f, "sine-large");
}

#ifdef INTERVAL_FLOAT_MPFR_ORACLE
// These witnesses have only normal or zero finite results. Precision 24 with nearest
// rounding therefore reproduces binary32 at every stage, without needing MPFR
// exponent-range/subnormal emulation. MPFR is confined to this optional target.
class FloatOracle {
  public:
    mpfr_t a, b, result;
    FloatOracle() { mpfr_inits2(24, a, b, result, (mpfr_ptr)nullptr); }
    ~FloatOracle() { mpfr_clears(a, b, result, (mpfr_ptr)nullptr); }
    FloatOracle(const FloatOracle&) = delete;
    FloatOracle& operator=(const FloatOracle&) = delete;
};

// Reference the target operation sequence, including each binary32 rounding;
// do not compute the entire expression at high precision and round only once.
double oracleValue(std::string_view name)
{
    FloatOracle ref;
    if (name == "addition") {
        mpfr_set_ui(ref.a, 16777216, MPFR_RNDN);
        mpfr_add_ui(ref.result, ref.a, 1, MPFR_RNDN);
    } else if (name == "multiplication" || name == "delay") {
        mpfr_set_ui(ref.a, 4097, MPFR_RNDN);
        mpfr_mul(ref.result, ref.a, ref.a, MPFR_RNDN);
        if (name == "delay") {
            mpfr_sub_ui(ref.result, ref.result, 16785408, MPFR_RNDN);
            mpfr_sub_ui(ref.result, ref.result, 1, MPFR_RNDN);
        }
    } else if (name == "conversion") {
        mpfr_set_ui(ref.result, 16777217, MPFR_RNDN);
    } else if (name == "contraction") {
        mpfr_set_d(ref.a, double(1.0f + 0x1p-13f), MPFR_RNDN);
        mpfr_set_d(ref.b, double(1.0f - 0x1p-13f), MPFR_RNDN);
        mpfr_set_si(ref.result, -1, MPFR_RNDN);
        mpfr_fma(ref.result, ref.a, ref.b, ref.result, MPFR_RNDN);
        mpfr_mul_2ui(ref.result, ref.result, 26, MPFR_RNDN);
    } else {
        mpfr_set_d(ref.a, name == "sine-large" ? double(1.0e20f) : 100.0, MPFR_RNDN);
        mpfr_sin(ref.result, ref.a, MPFR_RNDN);
    }
    return mpfr_get_d(ref.result, MPFR_RNDN);
}
#endif

}  // namespace

// Test driver: require inclusion for the original counterexamples and the large
// sine diagnostic. Every excluded value or invalid witness now fails CTest too.
int main(int argc, char** argv)
{
    std::string_view selected;
    for (int i = 1; i < argc; ++i) {
        const std::string_view argument = argv[i];
        if (argument == "--require-inclusion") {}
        else if (argument == "--case" && i + 1 < argc) selected = argv[++i];
        else {
            std::cerr << "Usage: " << argv[0]
                      << " [--require-inclusion] [--case NAME]\n";
            return 2;
        }
    }
    static_assert(std::numeric_limits<float>::is_iec559 && std::numeric_limits<float>::digits == 24,
                  "Float witnesses require IEEE binary32.");
#ifndef __EMSCRIPTEN__
    if (std::fegetround() != FE_TONEAREST) {
        std::cerr << "These witnesses require nearest-even rounding.\n";
        return 2;
    }
#endif
    itv::programPrecision() = 1;
    itv::libmCompensation() = true;
    const interval_algebra algebra;
    // Build witnesses lazily so --case isolates a rule. Large sine now belongs
    // to the default suite: the conservative path must avoid precision-metadata UB.
    struct Case {
        std::string_view name;
        Witness (*build)(const interval_algebra&);
    };
    const std::array<Case, 7> cases{{
        {"addition", addition}, {"multiplication", multiplication},
        {"delay", delay}, {"conversion", conversion},
        {"contraction", contraction}, {"sine", sine}, {"sine-large", largeSine}}};
    std::cout << std::setprecision(17);
    int tested = 0, failures = 0;
    for (const auto& test : cases) {
        if (!selected.empty() && selected != test.name) continue;
        const Witness witness = test.build(algebra);
        ++tested;
        if (!witness.validInputs || !std::isfinite(witness.observed)) {
            std::cerr << "INVALID WITNESS: " << witness.name << '\n';
            return 2;
        }
        bool included = witness.predicted.has(witness.observed);
#ifdef INTERVAL_FLOAT_MPFR_ORACLE
        const double reference = oracleValue(witness.name);
        // Elementary operations/FMA have a known correct-rounding contract.
        // The host sine may have its own error, so compare its observed value
        // and the independent correctly-rounded reference separately.
        if (witness.name != "sine" && witness.name != "sine-large" && reference != witness.observed) {
            std::cerr << "INVALID TARGET ROUNDING: " << witness.name << '\n';
            return 2;
        }
        included = included && witness.predicted.has(reference);
#endif
        failures += !included;
        std::cout << (included ? "PASS inclusion: " : "FAIL inclusion: ")
                  << witness.name << "\n  interval=[" << witness.predicted.lo()
                  << ", " << witness.predicted.hi() << "] observed=" << witness.observed
                  << " observed_included=" << witness.predicted.has(witness.observed);
#ifdef INTERVAL_FLOAT_MPFR_ORACLE
        std::cout << " mpfr24=" << reference
                  << " mpfr_included=" << witness.predicted.has(reference);
#endif
        std::cout << '\n';
    }
    if (!tested) {
        std::cerr << "Unknown witness: " << selected << '\n';
        return 2;
    }
    std::cout << tested << " checks, " << failures << " failures (required float inclusion).\n";
    return failures ? 1 : 0;
}
