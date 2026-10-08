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

// Each witness checks inclusion rather than freezing the library's erroneous
// endpoints. Premises concern the inputs/target execution only; a future wider
// sound result will satisfy the ordinary check and fail the known-gap mode.
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

// The small argument isolates enclosure failure without triggering the separate
// precision-metadata overflow exposed by the opt-in large-argument diagnostic.
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

// At this magnitude exactPrecisionUnary subtracts indistinguishable evaluations,
// then casts log2(0) to int. UBSan must expose that bug, rather than suppress it.
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

// Test driver: normal execution demands inclusion and currently fails. The
// explicit known-gap mode keeps CTest useful while documenting these defects;
// invalid inputs remain errors, and a repaired rule must leave the known-gap set.
int main(int argc, char** argv)
{
    bool expectGaps = false;
    std::string_view selected;
    for (int i = 1; i < argc; ++i) {
        const std::string_view argument = argv[i];
        if (argument == "--expect-known-gaps") expectGaps = true;
        else if (argument == "--require-inclusion") expectGaps = false;
        else if (argument == "--case" && i + 1 < argc) selected = argv[++i];
        else {
            std::cerr << "Usage: " << argv[0]
                      << " [--require-inclusion|--expect-known-gaps] [--case NAME]\n";
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
    // Build witnesses lazily so --case really isolates a rule. The large sine
    // diagnostic is opt-in: its analyzer UB must not taint the six inclusion tests.
    struct Case {
        std::string_view name;
        Witness (*build)(const interval_algebra&);
        bool optIn;
    };
    const std::array<Case, 7> cases{{
        {"addition", addition, false}, {"multiplication", multiplication, false},
        {"delay", delay, false}, {"conversion", conversion, false},
        {"contraction", contraction, false}, {"sine", sine, false},
        {"sine-large", largeSine, true}}};
    std::cout << std::setprecision(17);
    int tested = 0, failures = 0;
    for (const auto& test : cases) {
        if ((!selected.empty() && selected != test.name) || (selected.empty() && test.optIn)) continue;
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
        const bool success = expectGaps ? !included : included;
        failures += !success;
        std::cout << (expectGaps ? (success ? "KNOWN GAP: " : "GAP NOW COVERED: ")
                                : (success ? "PASS inclusion: " : "FAIL inclusion: "))
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
    std::cout << tested << " checks, " << failures << " failures ("
              << (expectGaps ? "expected known gaps; this does not certify inclusion"
                             : "required float inclusion") << ").\n";
    return failures ? 1 : 0;
}
