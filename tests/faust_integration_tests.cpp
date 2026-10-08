#include <array>
#include <cstdint>
#include <limits>
#include <random>
#include <string>

#include "interval/affine_ops.hh"
#include "interval/check.hh"
#include "interval/cxx_compat.hh"

using namespace itv;

namespace {

using Binary = interval (interval_algebra::*)(const interval&, const interval&) const;
using Unary = interval (interval_algebra::*)(const interval&) const;

// Runtime witnesses wrap in unsigned arithmetic, then copy the bit pattern to
// int32. They do not rely on the analyzer's implementation or signed overflow.
int32_t wrap(uint32_t bits) { return compat::bit_cast<int32_t>(bits); }
int32_t noiseStep(int32_t state)
{
    return wrap(UINT32_C(1103515245) * (UINT32_C(12345) + uint32_t(state)));
}

void precisionTests(int precision)
{
    programPrecision() = precision;
    const std::string prefix = precision == 1 ? "float: " : "double: ";
    const interval_algebra a;
    const affine_algebra affine(2147483648.0);
    const interval ints(-3, 7, 0), floats(-3, 7, -24);

    // Integral-valued floats must keep their nature through zeros, large values,
    // casts and rounding; otherwise a following multiply silently becomes int32.
    check(prefix + "floating zero keeps its nature", true,
          a.FloatNum(0).lsb() < 0 && interval(0, 0, -24).lsb() < 0 &&
          a.Sub(a.FloatNum(1), a.FloatNum(1)).lsb() < 0);
    check(prefix + "integer zero retains its requested grid", true,
          interval(0, 0, 8).lsb() == 8);
    const Unary floatRules[] = {&interval_algebra::FloatCast, &interval_algebra::Round,
        &interval_algebra::Rint, &interval_algebra::Floor, &interval_algebra::Ceil,
        &interval_algebra::Sqrt, &interval_algebra::SinBounds, &interval_algebra::CosBounds,
        &interval_algebra::ExpBounds, &interval_algebra::LogBounds};
    for (size_t i = 0; i < std::size(floatRules); ++i) {
        for (const interval& x : {a.IntNum(0), a.IntNum(1), a.FloatNum(0), a.FloatNum(1)}) {
            const interval result = (a.*floatRules[i])(x);
            check(prefix + "floating unary rule " + std::to_string(i), true,
                  result.isEmpty() || result.lsb() < 0);
        }
    }
    check(prefix + "rounded floating values do not wrap", true,
          a.Mul(a.Round(a.FloatNum(65536)), a.Round(a.FloatNum(65536))).has(4294967296.0));
    check(prefix + "mixed product keeps floating nature on a coarse grid", true,
          a.Mul(a.FloatNum(0.5), interval(0x1p30, 0x1p30, 30)).lsb() < 0);

    const Binary preserving[] = {&interval_algebra::Add, &interval_algebra::Sub,
        &interval_algebra::Mul, &interval_algebra::Mod, &interval_algebra::Min,
        &interval_algebra::Max};
    for (size_t i = 0; i < std::size(preserving); ++i) {
        const Binary op = preserving[i];
        check(prefix + "integer binary nature " + std::to_string(i), true,
              (a.*op)(ints, a.IntNum(2)).lsb() >= 0);
        check(prefix + "mixed binary nature " + std::to_string(i), true,
              (a.*op)(floats, a.IntNum(2)).lsb() < 0 &&
              (a.*op)(ints, a.FloatNum(2)).lsb() < 0);
    }
    check(prefix + "integer power stays integer", true,
          a.Pow(ints, a.IntNum(2)).lsb() >= 0);
    check(prefix + "floating power stays floating", true,
          a.Pow(a.FloatNum(0), a.IntNum(2)).lsb() < 0 &&
          a.Pow(a.FloatNum(2), a.IntNum(0)).lsb() < 0);
    const Binary integerResults[] = {&interval_algebra::And, &interval_algebra::Or,
        &interval_algebra::Xor, &interval_algebra::Eq, &interval_algebra::Ne,
        &interval_algebra::Lt, &interval_algebra::Le, &interval_algebra::Gt,
        &interval_algebra::Ge};
    for (size_t i = 0; i < std::size(integerResults); ++i)
        check(prefix + "integer result after floating inputs " + std::to_string(i), true,
              (a.*integerResults[i])(floats, a.FloatNum(2)).lsb() >= 0);
    check(prefix + "right shift keeps integer nature and truncated bounds", true,
          a.ARsh(a.IntNum(3), a.IntNum(1)).is(1) &&
          a.ARsh(a.IntNum(3), a.IntNum(1)).lsb() >= 0);

    // Infinite widening and affine horizon overflow both mean unknown int32,
    // including negative wrapped states. Check binary and conversion paths.
    const interval widened[] = {interval(0, HUGE_VAL, 0), interval(0, 1.5e22, 0),
        interval(-HUGE_VAL, -1, 0), interval(-HUGE_VAL, HUGE_VAL, 0)};
    for (const interval& x : widened) {
        check(prefix + "widened integer negation", true,
              a.Neg(x).has(INT_MIN) && a.Neg(x).has(INT_MAX) && a.Neg(x).lsb() >= 0);
        check(prefix + "widened integer absolute value keeps wrapped INT_MIN", true,
              a.Abs(x).has(INT_MIN) && a.Abs(x).lsb() >= 0);
        check(prefix + "widened integer conversion includes negative states", true,
              a.FloatCast(x).has(programBound(INT_MIN)) && a.FloatCast(x).lsb() < 0);
        check(prefix + "mixed widened sum includes negative states", true,
              a.Add(x, a.FloatNum(0.5)).has(programBound(INT_MIN) + 0.5));
        check(prefix + "affine mixed widened sum includes negative states", true,
              toItv(affine.Add(fromItv(x), affine.FloatNum(0.5)), affine.horizon())
                  .has(programBound(INT_MIN) + 0.5));
        check(prefix + "affine mixed widened product includes negative states", true,
              toItv(affine.Mul(fromItv(x), affine.FloatNum(0.5)), affine.horizon())
                  .has(programBound(INT_MIN) * 0.5));
        check(prefix + "mixed widened comparison is undecided", true,
              a.Lt(x, a.FloatNum(0)).has(0) && a.Lt(x, a.FloatNum(0)).has(1));
        check(prefix + "widened integer bitwise and keeps negative states", true,
              a.And(x, a.IntNum(-1)).has(-1));
    }

    // Point and corridor witnesses exercise all legal shift counts and both
    // signed orders, including INT_MIN. Unsigned left shift is the wrap oracle.
    std::mt19937 random(0x7261696eu);
    bool shiftIncluded = true;
    for (int i = 0; i < 128; ++i) {
        const int32_t value = wrap(random());
        const int32_t other = wrap(random());
        const interval corridor(std::min(value, other), std::max(value, other), 0);
        for (int count = 0; count < 32; ++count) {
            const interval k = a.IntNum(count);
            const int32_t left = wrap(uint32_t(value) << count);
            const int32_t arithmetic = value >> count;
            const int32_t logical = wrap(uint32_t(value) >> count);
            for (const interval& x : {a.IntNum(value), corridor})
                shiftIncluded = shiftIncluded && a.Lsh(x, k).has(left) &&
                    a.ARsh(x, k).has(arithmetic) && a.LRsh(x, k).has(logical) &&
                    a.Lsh(x, k).lsb() >= 0 && a.ARsh(x, k).lsb() >= 0;
        }
        shiftIncluded = shiftIncluded && a.Not(corridor).has(wrap(~uint32_t(value)));
    }
    check(prefix + "all int32 shift counts and complement include runtime witnesses", true,
          shiftIncluded);
    check(prefix + "full int32 complement terminates and keeps both ends", true,
          a.Not(interval(INT_MIN, INT_MAX, 0)) == interval(INT_MIN, INT_MAX, 0));
    check(prefix + "pure integer comparisons do not round to binary32", true,
          a.Gt(a.IntNum(16777217), a.IntNum(16777216)).is(1));
    check(prefix + "mixed comparisons use target conversion", true,
          a.Eq(a.IntNum(16777217), a.FloatNum(16777216)).is(precision == 1 ? 1 : 0));
    const interval cast = toItv(affine.FloatCast(affine.IntNum(16777217)), affine.horizon());
    check(prefix + "affine constant float cast retains floating nature", true,
          cast.lsb() < 0 && cast.is(programBound(16777217)));
    for (const interval& x : {ints, floats, a.IntNum(0), a.FloatNum(0)}) {
        const AffItv delayed = affine.Delay(fromItv(x), affine.IntNum(1));
        check(prefix + "affine delay preserves nature", true,
              (delayed.lsb >= 0) == (x.lsb() >= 0));
    }
    const AffItv zero = affine.IntNum(0);
    check(prefix + "soundfile metadata stays integer", true,
          affine.SoundFile(zero).lsb >= 0 && affine.SoundFileRate(zero, zero).lsb >= 0 &&
          affine.SoundFileLength(zero, zero).lsb >= 0);

    // Replay the LCG recursion and the second intermediate noise from rain.dsp.
    // A post-widening transfer check certifies the invariant; runtime unsigned
    // witnesses check both noise stages, not just the final recursive state.
    AffItv state = affine.IntNum(0);
    auto step = [&](const AffItv& x) {
        return affine.Mul(affine.IntNum(1103515245),
                          affine.Add(affine.IntNum(12345), affine.Mem(x)));
    };
    for (int i = 0; i < 16; ++i) state = awiden(state, step(state), affine.horizon());
    check(prefix + "LCG fixed point is inductive", true,
          toItv(step(state), affine.horizon()) <= toItv(state, affine.horizon()));
    const AffItv second = step(step(state));
    const interval noise = toItv(affine.Mul(affine.FloatCast(second),
                                           affine.FloatNum(0x1p-31)), affine.horizon());
    const interval predicate = a.Lt(a.Abs(noise), a.FloatNum(1.5));
    bool included = state.lsb >= 0 && second.lsb >= 0 && predicate.is(1);
    int32_t runtime = 0;
    for (int i = 0; i < 4096; ++i) {
        runtime = noiseStep(runtime);
        const int32_t intermediate = noiseStep(runtime);
        const double value = precision == 1 ? double(float(intermediate) * 0x1p-31f)
                                             : double(intermediate) * 0x1p-31;
        included = included && toItv(state, affine.horizon()).has(runtime) && noise.has(value);
    }
    check(prefix + "rain second noise remains included and its predicate true", true, included);
}

} // namespace

int main()
{
    // Exercise the C++17 fallbacks with endpoints which bit algorithms often miss.
    check("C++17: countl_zero includes zero and both end bits", true,
          compat::countl_zero(0) == 64 && compat::countl_zero(1) == 63 &&
          compat::countl_zero(UINT64_C(1) << 63) == 0);
    check("C++17: bit_cast preserves signed zero and subnormal bits", true,
          compat::bit_cast<uint32_t>(-0.0f) == UINT32_C(0x80000000) &&
          compat::bit_cast<float>(UINT32_C(1)) == std::numeric_limits<float>::denorm_min());
    precisionTests(1);
    precisionTests(2);
    return reportCheckResults();
}
