// Copyright 2026 Yann Orlarey, Stéphane Letz. Apache-2.0 (see LICENSE).
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <sstream>
#include <string>

#include "interval/affine_ops.hh"
#include "interval/check.hh"

#ifdef INTERVAL_INVALIDITY_MPFR_ORACLE
#include <mpfr.h>
#endif

using namespace itv;

namespace {
using Unary = interval (interval_algebra::*)(const interval&) const;
using Binary = interval (interval_algebra::*)(const interval&, const interval&) const;
using AUnary = AffItv (affine_algebra::*)(const AffItv&) const;
using ABinary = AffItv (affine_algebra::*)(const AffItv&, const AffItv&) const;

// Public transfer matrices check the observable propagation contract, including
// early-empty exits, rather than reproducing the implementation's predicates.
struct UnaryCase { const char* name; Unary ordinary; AUnary affine; };
const UnaryCase unaryCases[] = {
#define U(name) {#name, &interval_algebra::name, &affine_algebra::name}
    U(Abs), U(Neg), U(Inv), U(FloatCast), U(IntCast), U(Not), U(BitCast),
    U(Acos), U(Acosh), U(Asin), U(Asinh), U(Atan), U(Atanh), U(Cos), U(Cosh),
    U(Exp), U(Exp10), U(Log), U(Log10), U(Sin), U(Sinh), U(Sqrt), U(Tan), U(Tanh),
    U(Ceil), U(Floor), U(Rint), U(Round), U(Mem), U(Gen), U(Highest), U(Lowest),
    U(Input), U(Button), U(Checkbox), U(SoundFile)
#undef U
};
struct BinaryCase { const char* name; Binary ordinary; ABinary affine; };
const BinaryCase binaryCases[] = {
#define B(name) {#name, &interval_algebra::name, &affine_algebra::name}
    B(Add), B(Sub), B(Mul), B(Div), B(Mod), B(Fmod), B(Remainder), B(Pow), B(Atan2),
    B(Min), B(Max), B(And), B(Or), B(Xor), B(Lsh), B(ARsh), B(LRsh),
    B(Eq), B(Ne), B(Gt), B(Ge), B(Lt), B(Le), B(Delay), B(Prefix),
    B(FixPointUpdate), B(Attach), B(Enable), B(Control), B(RDTbl), B(Output)
#undef B
};

void representationTests()
{
    const interval valid(0, 2), invalid = valid.withInvalid(), only = empty().withInvalid();
    check("ordinary: true flag is sticky", true, invalid.withInvalid(false).mayBeInvalid());
    check("ordinary: empty numeric part is not an invalid execution", false, empty().mayBeInvalid());
    check("ordinary: NaN literal is invalid-only", true,
          interval(NAN).isEmpty() && interval(NAN).mayBeInvalid());
    check("ordinary: NaN endpoints are invalid-only", true,
          interval(NAN, NAN).isEmpty() && interval(NAN, NAN).mayBeInvalid());
    check("ordinary: a nonempty interval with an alert is not valid", false, invalid.isValid());
    check("ordinary: an alerted point is not a foldable constant", false,
          interval(1).withInvalid().isconst());
    check("ordinary: singleton shortcuts cannot erase an alert", true,
          !interval(0).withInvalid().isZero() && !interval(1).withInvalid().is(1));
    check("ordinary: equality observes flag-only changes", false, valid == invalid);
    check("ordinary: distinct empty states", false, empty() == only);
    check("ordinary: empty metadata does not perform an invalid cast", true,
          empty().msb() == 0 && only.msb() == 0 && !only.ispowerof2() && !only.isbitmask());
    check("ordinary: inclusion permits adding an alert", true, valid <= invalid);
    check("ordinary: inclusion forbids removing an alert", false, invalid <= valid);
    check("ordinary: invalid-only is included in an alerted numeric corridor", true, only <= invalid);
    check("ordinary: invalid-only is not numeric bottom", false, only <= empty());
    checkExact("ordinary: join does not discard invalid-only", reunion(valid, only), invalid);
    checkExact("ordinary: join does not discard invalid-only on the left", reunion(only, valid), invalid);
    check("ordinary: refinement cannot erase invalidity", true,
          intersection(invalid, interval(0, 1)).mayBeInvalid());
    check("ordinary: disjoint refinement preserves alert", true,
          intersection(invalid, interval(3, 4)).isEmpty() &&
          intersection(invalid, interval(3, 4)).mayBeInvalid());
    check("ordinary: diagnostic string includes invalidity", true,
          only.to_string().find("invalid") != std::string::npos);
    std::ostringstream stream; stream << invalid;
    check("ordinary: printed interval includes invalidity", true,
          stream.str().find("invalid") != std::string::npos);
    check("ordinary: libm compensation retains an alert", true,
          libmBounds(invalid, -HUGE_VAL, HUGE_VAL).mayBeInvalid());

    const AffItv a = fromItv(valid), bad = fromItv(invalid), none = fromItv(only);
    check("affine: invalid-only survives both bridges", true,
          none.isEmpty() && none.mayBeInvalid && toItv(none, 32).mayBeInvalid());
    checkExact("affine: flagged corridor round-trip", toItv(bad, 32), invalid);
    check("affine: inclusion observes flag-only change", false, aleq(bad, a, 32));
    check("affine: inclusion permits alert", true, aleq(a, bad, 32));
    check("affine: invalid-only is not bottom", false, aleq(none, aempty(), 32));
    check("affine: join with invalid-only", true, ajoin(a, none, 32).mayBeInvalid);
    check("affine: widen keeps an old alert", true, awiden(bad, a, 32).mayBeInvalid);
    check("affine: widen detects a newly raised alert", true, awiden(a, bad, 32).mayBeInvalid);
    check("affine: empty widening retains an old alert", true, awiden(none, a, 32).mayBeInvalid);
    check("affine: wrapping hull retains an alert", true,
          toItv(AffItv{0, 0, 0, 1e9, 0, true}, 32).mayBeInvalid());
}

void propagationTests(int precision)
{
    programPrecision() = precision;
    const std::string prefix = precision == 1 ? "float: " : "double: ";
    const interval_algebra a;
    const affine_algebra affine(32);
    const interval inputs[] = {a.FloatNum(0.5).withInvalid(), empty().withInvalid()};
    for (const auto& method : unaryCases) {
        bool ordinary = true, temporal = true;
        for (const interval& x : inputs) {
            ordinary = ordinary && (a.*method.ordinary)(x).mayBeInvalid();
            temporal = temporal && (affine.*method.affine)(fromItv(x)).mayBeInvalid;
        }
        check(prefix + method.name + " preserves alerts in ordinary and empty inputs", true, ordinary);
        check(prefix + method.name + " preserves alerts in affine and empty inputs", true, temporal);
        const interval value = a.FloatNum(0.5);
        const interval baseline = (a.*method.ordinary)(value);
        const interval alerted = (a.*method.ordinary)(value.withInvalid());
        check(prefix + method.name + " adding an alert retains the numeric image", true,
              baseline.isEmpty() || (!alerted.isEmpty() && alerted.has(baseline.lo()) && alerted.has(baseline.hi())));
        const interval affineBaseline = toItv((affine.*method.affine)(fromItv(value)), 32);
        const interval affineAlerted = toItv((affine.*method.affine)(fromItv(value.withInvalid())), 32);
        check(prefix + method.name + " adding an affine alert retains the numeric image", true,
              affineBaseline.isEmpty() || (!affineAlerted.isEmpty() &&
                  affineAlerted.has(affineBaseline.lo()) && affineAlerted.has(affineBaseline.hi())));
    }
    const interval good = a.IntNum(1);
    for (const auto& method : binaryCases) {
        bool ordinary = true, temporal = true;
        for (const interval& x : {good.withInvalid(), empty().withInvalid()}) {
            ordinary = ordinary && (a.*method.ordinary)(x, good).mayBeInvalid() &&
                (a.*method.ordinary)(good, x).mayBeInvalid();
            temporal = temporal && (affine.*method.affine)(fromItv(x), fromItv(good)).mayBeInvalid &&
                (affine.*method.affine)(fromItv(good), fromItv(x)).mayBeInvalid;
        }
        check(prefix + method.name + " propagates each ordinary operand", true, ordinary);
        check(prefix + method.name + " propagates each affine operand", true, temporal);
        // An alert adds possible invalid executions without removing valid ones.
        // Exercise both numeric natures because singleton shortcuts differ.
        bool retainsOrdinary = true, retainsAffine = true;
        for (const interval& value : {good, a.FloatNum(0.5)}) {
            const interval baseline = (a.*method.ordinary)(value, good);
            const interval affineBaseline = toItv((affine.*method.affine)(fromItv(value), fromItv(good)), 32);
            for (int position = 0; position < 2; ++position) {
                const interval left = position == 0 ? value.withInvalid() : value;
                const interval right = position == 1 ? good.withInvalid() : good;
                const interval alerted = (a.*method.ordinary)(left, right);
                const interval affineAlerted = toItv((affine.*method.affine)(fromItv(left), fromItv(right)), 32);
                retainsOrdinary = retainsOrdinary && (baseline.isEmpty() ||
                    (!alerted.isEmpty() && alerted.has(baseline.lo()) && alerted.has(baseline.hi())));
                retainsAffine = retainsAffine && (affineBaseline.isEmpty() ||
                    (!affineAlerted.isEmpty() && affineAlerted.has(affineBaseline.lo()) &&
                     affineAlerted.has(affineBaseline.hi())));
            }
        }
        check(prefix + method.name + " adding an alert retains the ordinary numeric image", true, retainsOrdinary);
        check(prefix + method.name + " adding an alert retains the affine numeric image", true, retainsAffine);
    }
    const Unary bounds[] = {&interval_algebra::AcosBounds, &interval_algebra::AcoshBounds,
        &interval_algebra::AsinBounds, &interval_algebra::AsinhBounds, &interval_algebra::AtanBounds,
        &interval_algebra::AtanhBounds, &interval_algebra::CosBounds, &interval_algebra::CoshBounds,
        &interval_algebra::ExpBounds, &interval_algebra::LogBounds, &interval_algebra::Log10Bounds,
        &interval_algebra::SinBounds, &interval_algebra::SinhBounds, &interval_algebra::TanBounds,
        &interval_algebra::TanhBounds};
    bool references = true;
    for (auto method : bounds)
        for (const interval& x : inputs) references = references && (a.*method)(x).mayBeInvalid();
    check(prefix + "reference Bounds entry points also propagate alerts", true, references);
    check(prefix + "reference binary Bounds also propagate alerts", true,
          a.PowBounds(inputs[0], good).mayBeInvalid() &&
          a.Atan2Bounds(good, inputs[1]).mayBeInvalid());
}

void domainTests(int precision)
{
    programPrecision() = precision;
    const interval_algebra a;
    const affine_algebra affine(16);
    const std::string p = precision == 1 ? "float domains: " : "double domains: ";
    const interval zero = a.FloatNum(0), one = a.FloatNum(1), inf = a.FloatNum(HUGE_VAL);
    const interval negInf = a.FloatNum(-HUGE_VAL);
    check(p + "sqrt partially invalid keeps valid image and alert", true,
          a.Sqrt(interval(-1, 4)).has(0) && a.Sqrt(interval(-1, 4)).has(2) &&
          a.Sqrt(interval(-1, 4)).mayBeInvalid());
    check(p + "sqrt invalid-only remains distinct from bottom", true,
          a.Sqrt(interval(-4, -1)).isEmpty() && a.Sqrt(interval(-4, -1)).mayBeInvalid());
    check(p + "sqrt valid domain has no alert", false, a.Sqrt(interval(0, 4)).mayBeInvalid());
    check(p + "sqrt infinity is valid", false, a.Sqrt(inf).mayBeInvalid());
    check(p + "negative log input is invalid", true, a.Log(interval(-1, 2)).mayBeInvalid());
    check(p + "log zero is numeric -infinity, not NaN", true,
          a.Log(zero).has(-HUGE_VAL) && !a.Log(zero).mayBeInvalid());
    check(p + "log10 invalid-only domain", true, a.Log10(negInf).mayBeInvalid());
    for (Unary method : {&interval_algebra::Asin, &interval_algebra::Acos, &interval_algebra::Atanh}) {
        check(p + "unit domain accepts endpoints", false, (a.*method)(interval(-1, 1)).mayBeInvalid());
        check(p + "unit domain rejects excluded values", true, (a.*method)(interval(-2, 2)).mayBeInvalid());
    }
    check(p + "acosh domain lower boundary", false, a.Acosh(one).mayBeInvalid());
    check(p + "acosh domain violation", true, a.Acosh(interval(0, 2)).mayBeInvalid());
    for (Unary method : {&interval_algebra::Sin, &interval_algebra::Cos, &interval_algebra::Tan}) {
        check(p + "trigonometric infinity is invalid-only", true,
              (a.*method)(inf).isEmpty() && (a.*method)(inf).mayBeInvalid());
        check(p + "trigonometric unbounded input keeps finite results and alert", true,
              !(a.*method)(interval(-HUGE_VAL, HUGE_VAL)).isEmpty() &&
              (a.*method)(interval(-HUGE_VAL, HUGE_VAL)).mayBeInvalid());
    }
    check(p + "atan infinity is valid", false, a.Atan(inf).mayBeInvalid());
    check(p + "exp infinity is valid", false, a.Exp(inf).mayBeInvalid());
    check(p + "exp negative infinity is exact zero", true,
          a.Exp(negInf).has(0) && !a.Exp(negInf).mayBeInvalid());
    check(p + "opposite infinities cannot be added", true, a.Add(inf, negInf).mayBeInvalid());
    check(p + "same infinities can be added", false, a.Add(inf, inf).mayBeInvalid());
    check(p + "same infinities cannot be subtracted", true, a.Sub(inf, inf).mayBeInvalid());
    check(p + "zero times infinity is invalid", true, a.Mul(zero, inf).mayBeInvalid());
    check(p + "nonzero times infinity is valid", false, a.Mul(one, inf).mayBeInvalid());
    check(p + "zero over zero is invalid", true, a.Div(zero, zero).mayBeInvalid());
    check(p + "infinity over infinity is invalid", true, a.Div(inf, inf).mayBeInvalid());
    check(p + "nonzero over zero is infinite but not NaN", false, a.Div(one, zero).mayBeInvalid());
    check(p + "reciprocal zero is infinite but not NaN", false, a.Inv(zero).mayBeInvalid());
    check(p + "merged signed zero includes both reciprocal infinities", true,
          a.Inv(zero).has(-HUGE_VAL) && a.Inv(zero).has(HUGE_VAL));
    for (Binary method : {static_cast<Binary>(&interval_algebra::Mod), &interval_algebra::Fmod,
                          &interval_algebra::Remainder}) {
        check(p + "remainder with zero divisor is invalid", true, (a.*method)(one, zero).mayBeInvalid());
        check(p + "infinite dividend is invalid", true, (a.*method)(inf, one).mayBeInvalid());
        check(p + "infinite divisor and finite dividend is valid", false, (a.*method)(one, inf).mayBeInvalid());
    }
    check(p + "integer modulo with possible zero divisor", true,
          a.Mod(a.IntNum(5), interval(0, 3, 0)).mayBeInvalid());
    check(p + "scalar modulo overload preserves a domain violation", true,
          a.Mod(a.IntNum(5), 0.0).mayBeInvalid());
    check(p + "INT_MIN modulo -1 is undefined", true,
          a.Mod(a.IntNum(INT_MIN), a.IntNum(-1)).mayBeInvalid());
    check(p + "wrapping addition is valid", true,
          a.Add(a.IntNum(INT_MAX), a.IntNum(1)).has(INT_MIN) &&
          !a.Add(a.IntNum(INT_MAX), a.IntNum(1)).mayBeInvalid());
    check(p + "wrapping left shift is valid under the contract", false,
          a.Lsh(a.IntNum(INT_MIN), a.IntNum(1)).mayBeInvalid());
    for (Binary method : {&interval_algebra::Lsh, &interval_algebra::ARsh, &interval_algebra::LRsh}) {
        check(p + "negative shift count is invalid", true, (a.*method)(one, a.IntNum(-1)).mayBeInvalid());
        check(p + "shift count 32 is invalid", true, (a.*method)(one, a.IntNum(32)).mayBeInvalid());
        check(p + "legal shift count has no alert", false, (a.*method)(one, a.IntNum(31)).mayBeInvalid());
        check(p + "invalid conversion before a shift is retained", true, (a.*method)(inf, one).mayBeInvalid());
    }
    check(p + "negative base fractional exponent is invalid", true,
          a.Pow(a.FloatNum(-2), a.FloatNum(0.5)).mayBeInvalid());
    check(p + "negative base integer exponent is valid", false,
          a.Pow(a.FloatNum(-2), a.FloatNum(3)).mayBeInvalid());
    check(p + "integer power cannot silently treat a negative exponent as zero", true,
          a.Pow(a.IntNum(2), a.IntNum(-1)).mayBeInvalid());
    check(p + "floating negative powers remain valid", false,
          a.Pow(a.FloatNum(2), a.FloatNum(-1)).mayBeInvalid());
    check(p + "zero to a negative power is infinite but not NaN", false,
          a.Pow(zero, a.FloatNum(-1)).mayBeInvalid());
    check(p + "merged signed zero includes negative odd-power infinity", true,
          a.Pow(zero, a.FloatNum(-1)).has(-HUGE_VAL));
    check(p + "negative delay amount is invalid", true,
          a.Delay(one, a.IntNum(-1)).mayBeInvalid());
    check(p + "invalid delay conversion is invalid", true, a.Delay(one, inf).mayBeInvalid());
    check(p + "assertion cannot hide excluded inputs", true,
          a.AssertBounds(a.FloatNum(-1), one, interval(-2, 2)).mayBeInvalid());
    check(p + "uncertain assertion bounds check all argument tuples", true,
          a.AssertBounds(interval(0, 2), a.FloatNum(5), interval(1, 4)).mayBeInvalid());
    check(p + "uncertain assertion upper bound check", true,
          a.AssertBounds(zero, interval(3, 5), interval(1, 4)).mayBeInvalid());
    check(p + "overflow then subtraction detects an emergent NaN", true,
          a.Sub(a.Mul(a.FloatNum(std::numeric_limits<double>::max()), a.FloatNum(2)),
                a.Mul(a.FloatNum(std::numeric_limits<double>::max()), a.FloatNum(2))).mayBeInvalid());
    check(p + "affine moving conversion crossing int32 is invalid", true,
          affine.IntCast(AffItv{0, 0, 0, 1e9, -24}).mayBeInvalid);
    check(p + "affine opposite infinities are invalid", true,
          affine.Add(fromItv(inf), fromItv(negInf)).mayBeInvalid);
    check(p + "affine negative delay is invalid", true,
          affine.Delay(fromItv(one), affine.IntNum(-1)).mayBeInvalid);
    check(p + "affine integer wrapping is not invalid", false,
          affine.IntCast(AffItv{0, 0, 0, 1e9, 0}).mayBeInvalid);
}

void conversionTests(int precision)
{
    programPrecision() = precision;
    const interval_algebra a;
    const std::string p = precision == 1 ? "float casts: " : "double casts: ";
    for (double value : {HUGE_VAL, -HUGE_VAL, 1e30, -1e30, double(NAN)}) {
        const interval result = a.IntCast(a.FloatNum(value));
        check(p + "invalid-only cast has no valid integer result", true,
              result.isEmpty() && result.mayBeInvalid() && result.lsb() >= 0);
    }
    check(p + "partial cast keeps valid truncated image and alert", true,
          a.IntCast(interval(0, 1e30)).has(0) && a.IntCast(interval(0, 1e30)).has(INT_MAX) &&
          a.IntCast(interval(0, 1e30)).mayBeInvalid());
    check(p + "ordinary fractional truncation is valid", true,
          a.IntCast(interval(-3.8, 4.9)).has(-3) && a.IntCast(interval(-3.8, 4.9)).has(4) &&
          !a.IntCast(interval(-3.8, 4.9)).mayBeInvalid());
    check(p + "NaN flag persists through a float cast", true,
          a.FloatCast(a.FloatNum(NAN)).mayBeInvalid());
    check(p + "widened integer remains defined int32", false,
          a.IntCast(interval(-HUGE_VAL, HUGE_VAL, 0)).mayBeInvalid());
    check(p + "integer-to-float rounding can make a subsequent cast invalid", precision == 1,
          a.IntCast(a.FloatCast(a.IntNum(INT_MAX))).mayBeInvalid());
    if (precision == 2) {
        for (double value : {-2147483648.75, 2147483647.75,
                             std::nextafter(-2147483649.0, HUGE_VAL),
                             std::nextafter(2147483648.0, -HUGE_VAL)})
            check(p + "fraction just beyond int32 can still truncate validly", false,
                  a.IntCast(a.FloatNum(value)).mayBeInvalid());
    } else {
        check(p + "float just below upper excluded boundary is valid", false,
              a.IntCast(a.FloatNum(std::nextafter(2147483648.0f, -INFINITY))).mayBeInvalid());
        check(p + "float INT_MIN is valid", false, a.IntCast(a.FloatNum(INT_MIN)).mayBeInvalid());
    }
    for (double value : {-2147483649.0, 2147483648.0})
        if (precision == 2 || value > 0)
            check(p + "excluded truncation boundary is invalid", true,
                  a.IntCast(a.FloatNum(value)).mayBeInvalid());
    if (precision == 1) {
        check(p + "rounded lower boundary becomes valid INT_MIN", false,
              a.IntCast(a.FloatNum(-2147483649.0)).mayBeInvalid());
        check(p + "first float below INT_MIN is invalid", true,
              a.IntCast(a.FloatNum(std::nextafter(float(INT_MIN), -INFINITY))).mayBeInvalid());
    }
}

void structuralTests(int precision)
{
    programPrecision() = precision;
    const interval_algebra a;
    const affine_algebra affine(32);
    const std::string p = precision == 1 ? "float structural: " : "double structural: ";
    const interval label = a.Label("test"), lo = a.FloatNum(0), hi = a.FloatNum(2);
    const interval init = a.FloatNum(1), step = a.FloatNum(0.5), alert = init.withInvalid();
    check(p + "label is numeric bottom, not NaN", false, label.mayBeInvalid());
    check(p + "literal NaN injection alerts", true, a.FloatNum(NAN).mayBeInvalid());
    check(p + "UI alerts propagate", true,
          a.VSlider(label, alert, lo, hi, step).mayBeInvalid() &&
          a.HSlider(label, init, lo, hi, step.withInvalid()).mayBeInvalid() &&
          a.NumEntry(label, init, lo.withInvalid(), hi, step).mayBeInvalid());
    using Widget = interval (interval_algebra::*)(const interval&, const interval&,
        const interval&, const interval&, const interval&) const;
    using AWidget = AffItv (affine_algebra::*)(const AffItv&, const AffItv&,
        const AffItv&, const AffItv&, const AffItv&) const;
    const Widget widgets[] = {&interval_algebra::VSlider, &interval_algebra::HSlider,
                             &interval_algebra::NumEntry};
    const AWidget awidgets[] = {&affine_algebra::VSlider, &affine_algebra::HSlider,
                               &affine_algebra::NumEntry};
    bool alerts = true;
    const affine_algebra defaults(32, true);
    for (size_t k = 0; k < 3; ++k) {
        for (size_t input = 0; input < 5; ++input) {
            std::array<interval, 5> args{label, init, lo, hi, step};
            args[input] = args[input].withInvalid();
            alerts = alerts && (a.*widgets[k])(args[0], args[1], args[2], args[3], args[4]).mayBeInvalid();
            for (const affine_algebra* mode : {&affine, &defaults})
                alerts = alerts && (mode->*awidgets[k])(fromItv(args[0]), fromItv(args[1]),
                    fromItv(args[2]), fromItv(args[3]), fromItv(args[4])).mayBeInvalid;
        }
    }
    check(p + "every widget operand retains alerts, including default-parameter mode", true, alerts);
    check(p + "a widget default outside its declared range is invalid", true,
          a.HSlider(label, a.FloatNum(3), lo, hi, step).mayBeInvalid());
    check(p + "an infinite widget step is invalid without an analyzer cast overflow", true,
          a.HSlider(label, init, lo, hi, a.FloatNum(HUGE_VAL)).mayBeInvalid());
    check(p + "uncertain widget bounds are a joint-domain obligation", true,
          a.HSlider(label, init, interval(0, 2), a.FloatNum(5), step).mayBeInvalid());
    check(p + "selector invalidity cannot be dropped", true,
          a.Select2(alert, lo, hi).mayBeInvalid() &&
          affine.Select2(fromItv(alert), fromItv(lo), fromItv(hi)).mayBeInvalid);
    check(p + "waveform alert in any sample survives", true,
          a.Waveform({lo, empty().withInvalid(), hi}).mayBeInvalid() &&
          affine.Waveform({fromItv(lo), fromItv(alert)}).mayBeInvalid);
    check(p + "display-bound alerts propagate", true,
          a.HBargraph(label, alert, hi, lo).mayBeInvalid() &&
          affine.VBargraph(fromItv(label), fromItv(lo), fromItv(alert), fromItv(hi)).mayBeInvalid);
    check(p + "valid table write has no alert", false,
          a.WRTbl(a.IntNum(100), lo, a.IntNum(99), hi).mayBeInvalid());
    check(p + "negative table write index is invalid", true,
          a.WRTbl(a.IntNum(100), lo, interval(-50, 149, 0), hi).mayBeInvalid() &&
          affine.WRTbl(affine.IntNum(100), fromItv(lo), fromItv(interval(-50, 149, 0)),
                        fromItv(hi)).mayBeInvalid);
    check(p + "unknown table extent cannot certify a read", true,
          a.RDTbl(hi, a.IntNum(0)).mayBeInvalid());
    check(p + "unknown soundfile index extents are not silently certified", true,
          a.SoundFileBuffer(label, lo, lo, lo).mayBeInvalid() &&
          affine.SoundFileRate(fromItv(label), fromItv(lo)).mayBeInvalid);
    const interval wrappingState(INT_MIN, INT_MAX, 0);
    const interval wrappingIndex = a.Sub(a.Mod(wrappingState, a.IntNum(200)), a.IntNum(50));
    checkExact(p + "modulo outside recurrence retains the negative wrapping remainder",
               wrappingIndex, interval(-249, 149, 0));
    check(p + "wrapping arithmetic is valid but the subsequent table access is not", true,
          !wrappingIndex.mayBeInvalid() &&
          a.WRTbl(a.IntNum(100), lo, wrappingIndex, hi).mayBeInvalid());
    check(p + "unknown foreign function is not certified", true,
          a.ForeignFunction(0, {}).mayBeInvalid() && affine.ForeignFunction(0, {}).mayBeInvalid);
    check(p + "unknown foreign float can contain NaN", true,
          a.ForeignVar(1, label, label).mayBeInvalid() &&
          affine.ForeignConst(1, fromItv(label), fromItv(label)).mayBeInvalid);
    check(p + "foreign integer storage cannot contain NaN", false,
          a.ForeignConst(0, label, label).mayBeInvalid());
}

// Independent runtime witnesses use libm solely to observe NaN creation; the
// test does not assume that libm gives correctly rounded finite results.
void nanWitnessTests(int precision)
{
    programPrecision() = precision;
    const interval_algebra a;
    const std::string p = precision == 1 ? "float NaN witnesses: " : "double NaN witnesses: ";
    auto point = [&](double x) { return a.FloatNum(x); };
    struct Case { const char* name; Unary method; double (*function)(double); };
    const Case cases[] = {
        {"sqrt", &interval_algebra::Sqrt, static_cast<double(*)(double)>(std::sqrt)},
        {"acos", &interval_algebra::Acos, static_cast<double(*)(double)>(std::acos)},
        {"asin", &interval_algebra::Asin, static_cast<double(*)(double)>(std::asin)},
        {"acosh", &interval_algebra::Acosh, static_cast<double(*)(double)>(std::acosh)},
        {"atanh", &interval_algebra::Atanh, static_cast<double(*)(double)>(std::atanh)},
        {"log", &interval_algebra::Log, static_cast<double(*)(double)>(std::log)},
        {"log10", &interval_algebra::Log10, static_cast<double(*)(double)>(std::log10)},
        {"sin", &interval_algebra::Sin, static_cast<double(*)(double)>(std::sin)},
        {"cos", &interval_algebra::Cos, static_cast<double(*)(double)>(std::cos)},
        {"tan", &interval_algebra::Tan, static_cast<double(*)(double)>(std::tan)}
    };
    const double samples[] = {-HUGE_VAL, -4, -2, -1.25, -1, -0.5, -0.0, 0, 0.5, 1, 1.25, 2, 4, HUGE_VAL};
    for (const auto& method : cases) {
        bool included = true;
        for (double value : samples) {
            const interval input = point(value);
            // Values are exactly representable in both target formats here.
            if (std::isnan(method.function(value)))
                included = included && (a.*method.method)(input).mayBeInvalid();
        }
        check(p + method.name + " flags every observed NaN", true, included);
    }
    struct PairCase { const char* name; Binary method; double (*function)(double, double); };
    const PairCase pairs[] = {
        {"add", &interval_algebra::Add, [](double x, double y) { return x + y; }},
        {"sub", &interval_algebra::Sub, [](double x, double y) { return x - y; }},
        {"mul", &interval_algebra::Mul, [](double x, double y) { return x * y; }},
        {"div", &interval_algebra::Div, [](double x, double y) { return x / y; }},
        {"fmod", &interval_algebra::Fmod, static_cast<double(*)(double, double)>(std::fmod)},
        {"remainder", &interval_algebra::Remainder,
                      static_cast<double(*)(double, double)>(std::remainder)},
        {"pow", &interval_algebra::Pow, static_cast<double(*)(double, double)>(std::pow)},
        {"atan2", &interval_algebra::Atan2, static_cast<double(*)(double, double)>(std::atan2)}
    };
    for (const auto& method : pairs) {
        bool included = true;
        for (double x : samples) {
            for (double y : samples) {
                if (std::isnan(method.function(x, y)))
                    included = included && (a.*method.method)(point(x), point(y)).mayBeInvalid();
            }
        }
        check(p + method.name + " flags every NaN in the argument-pair matrix", true, included);
    }
    std::mt19937 random(0x1a11du);
    bool casts = true;
    for (int i = 0; i < 4096; ++i) {
        const uint32_t bits = random();
        const float value = compat::bit_cast<float>(bits);
        const bool invalid = !std::isfinite(value) || std::trunc(double(value)) < INT_MIN ||
                             std::trunc(double(value)) > INT_MAX;
        const interval result = a.IntCast(point(value));
        casts = casts && (!invalid || result.mayBeInvalid());
        if (!invalid) casts = casts && result.has(std::trunc(double(value)));
        // Never execute an invalid C++ conversion in the oracle itself.
    }
    check(p + "4096 independent float bit patterns certify only defined casts", true, casts);
}

#ifdef INTERVAL_INVALIDITY_MPFR_ORACLE
// MPFR is an optional independent domain oracle, not a production dependency.
// Its NaN classification is checked in addition to the host-libm witnesses.
void mpfrDomainTests(int precision)
{
    programPrecision() = precision;
    const interval_algebra a;
    const std::string p = precision == 1 ? "MPFR float: " : "MPFR double: ";
    mpfr_t x, y, result;
    mpfr_inits2(256, x, y, result, nullptr);
    using MPUnary = int (*)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
    using MPBinary = int (*)(mpfr_ptr, mpfr_srcptr, mpfr_srcptr, mpfr_rnd_t);
    struct UCase { const char* name; Unary method; MPUnary oracle; };
    const UCase unary[] = {
        {"sqrt", &interval_algebra::Sqrt, mpfr_sqrt},
        {"acos", &interval_algebra::Acos, mpfr_acos},
        {"acosh", &interval_algebra::Acosh, mpfr_acosh},
        {"asin", &interval_algebra::Asin, mpfr_asin},
        {"asinh", &interval_algebra::Asinh, mpfr_asinh},
        {"atan", &interval_algebra::Atan, mpfr_atan},
        {"atanh", &interval_algebra::Atanh, mpfr_atanh},
        {"sin", &interval_algebra::Sin, mpfr_sin},
        {"cos", &interval_algebra::Cos, mpfr_cos},
        {"tan", &interval_algebra::Tan, mpfr_tan},
        {"log", &interval_algebra::Log, mpfr_log},
        {"log10", &interval_algebra::Log10, mpfr_log10},
        {"exp", &interval_algebra::Exp, mpfr_exp}
    };
    const double samples[] = {-HUGE_VAL, -4, -2, -1.25, -1, -0.5, -0.0, 0, 0.5, 1, 1.25, 2, 4, HUGE_VAL};
    for (const auto& test : unary) {
        bool sound = true;
        for (double value : samples) {
            mpfr_set_d(x, value, MPFR_RNDN);
            test.oracle(result, x, MPFR_RNDN);
            if (mpfr_nan_p(result)) sound = sound && (a.*test.method)(a.FloatNum(value)).mayBeInvalid();
        }
        check(p + test.name + " never misses an MPFR NaN", true, sound);
    }
    struct BCase { const char* name; Binary method; MPBinary oracle; };
    const BCase binary[] = {
        {"add", &interval_algebra::Add, mpfr_add},
        {"sub", &interval_algebra::Sub, mpfr_sub},
        {"mul", &interval_algebra::Mul, mpfr_mul},
        {"div", &interval_algebra::Div, mpfr_div},
        {"pow", &interval_algebra::Pow, mpfr_pow},
        {"atan2", &interval_algebra::Atan2, mpfr_atan2},
        {"fmod", &interval_algebra::Fmod, mpfr_fmod},
        {"remainder", &interval_algebra::Remainder, mpfr_remainder}
    };
    for (const auto& test : binary) {
        bool sound = true;
        for (double vx : samples) {
            for (double vy : samples) {
                mpfr_set_d(x, vx, MPFR_RNDN); mpfr_set_d(y, vy, MPFR_RNDN);
                test.oracle(result, x, y, MPFR_RNDN);
                if (mpfr_nan_p(result))
                    sound = sound && (a.*test.method)(a.FloatNum(vx), a.FloatNum(vy)).mayBeInvalid();
            }
        }
        check(p + test.name + " never misses an MPFR NaN pair", true, sound);
    }
    std::mt19937 random(0xca570u);
    bool sound = true;
    for (int i = 0; i < 4096; ++i) {
        const double value = compat::bit_cast<float>(uint32_t(random()));
        mpfr_set_d(x, value, MPFR_RNDN);
        const bool exceptional = mpfr_nan_p(x) || mpfr_inf_p(x);
        if (!exceptional) mpfr_trunc(result, x);
        const bool outside = exceptional || mpfr_cmp_si(result, INT_MIN) < 0 ||
                             mpfr_cmp_si(result, INT_MAX) > 0;
        const interval enclosed = a.IntCast(a.FloatNum(value));
        sound = sound && (!outside || enclosed.mayBeInvalid());
        if (!outside) sound = sound && enclosed.has(mpfr_get_d(result, MPFR_RNDN));
    }
    check(p + "int32 conversion agrees with independent MPFR truncation", true, sound);
    mpfr_clears(x, y, result, nullptr);
}
#endif
}  // namespace

int main()
{
    representationTests();
    for (int precision : {1, 2}) {
        propagationTests(precision);
        domainTests(precision);
        conversionTests(precision);
        structuralTests(precision);
        nanWitnessTests(precision);
#ifdef INTERVAL_INVALIDITY_MPFR_ORACLE
        mpfrDomainTests(precision);
#endif
    }
    return reportCheckResults();
}
