#pragma once

#include "interval_def.hh"

#include "FaustAlgebra.hh"

namespace itv {
class interval_algebra : public FaustAlgebra<interval> {
   private:
    interval iPow(const interval& x, const interval& y) const;  // integer power, when x can be negative
    interval fPow(const interval& x, const interval& y) const;  // float power, when x is positive

   public:
    // Injections of external values
    interval IntNum(int x) const override;
    interval Int64Num(int64_t x) const override;
    // Public API: inject a known literal as a point, rounded once to binary32 in float
    // mode. Float/double retain floating nature even for integral-looking literals;
    // quad/fixed keep legacy injection. Explicit literals need no analyzer-error widening.
    // NaN keeps the library's historical empty representation.
    interval FloatNum(double x) const override;
    interval Label(const std::string& x) const override;

    // Missing operations
    interval FixPointUpdate(const interval& x, const interval& y) const override;
    interval Input(const interval& c) const override;
    interval Output(const interval& c, const interval& y) const override;
    interval HBargraph(const interval& name, const interval& lo, const interval& hi,
                       const interval& signal) const override;
    interval VBargraph(const interval& name, const interval& lo, const interval& hi,
                       const interval& signal) const override;
    interval Gen(const interval& x) const override;
    interval Attach(const interval& x, const interval& y) const override;
    interval Enable(const interval& x, const interval& control) const override;
    interval Control(const interval& x, const interval& control) const override;
    interval AssertBounds(const interval& lo, const interval& hi, const interval& x) const override;
    interval Highest(const interval& x) const override;
    interval Lowest(const interval& x) const override;
    interval BitCast(const interval& x) const override;
    interval Select2(const interval& x, const interval& y, const interval& z) const override;
    interval Prefix(const interval& x, const interval& y) const override;
    interval RDTbl(const interval& wtbl, const interval& ri) const override;
    interval WRTbl(const interval& n, const interval& g, const interval& wi,
                   const interval& ws) const override;
    interval SoundFile(const interval& label) const override;
    interval SoundFileRate(const interval& sf, const interval& x) const override;
    interval SoundFileLength(const interval& sf, const interval& x) const override;
    interval SoundFileBuffer(const interval& sf, const interval& x, const interval& y,
                             const interval& z) const override;
    interval Waveform(const std::vector<interval>& w) const override;

    // Foreign functions
    interval ForeignFunction(int resultNature, const std::vector<interval>& args) const override;
    interval ForeignVar(int declaredNature, const interval& name,
                        const interval& file) const override;
    interval ForeignConst(int declaredNature, const interval& name,
                          const interval& file) const override;

    // User interface elements
    interval Button(const interval& name) const override;
    interval Checkbox(const interval& name) const override;
    interval VSlider(const interval& name, const interval& init, const interval& lo,
                     const interval& hi, const interval& step) const override;
    interval HSlider(const interval& name, const interval& init, const interval& lo,
                     const interval& hi, const interval& step) const override;
    interval NumEntry(const interval& name, const interval& init, const interval& lo,
                      const interval& hi, const interval& step) const override;

    // Public API: enclose wrapping int32 absolute value or floating absolute value.
    // Empty stays empty. Normalize widened integer corridors before endpoint work;
    // the result retains the input nature, including integral-valued floating zero.
    interval Abs(const interval& x) const override;
    void     testAbs();
    //
    // Public API: enclose addition; float/double floating bounds round outward at each
    // endpoint, while the existing int32 wrapping and other precision paths remain.
    interval Add(const interval& x, const interval& y) const override;
    void     testAdd();
    //
    // Public API: enclose subtraction; float/double floating bounds round outward at each
    // endpoint, while the existing int32 wrapping and other precision paths remain.
    interval Sub(const interval& x, const interval& y) const override;
    void     testSub();
    //
    // Public API: enclose multiplication; float/double floating corners round outward.
    // Zero times an unbounded endpoint keeps the historical numeric-hull convention;
    // possible NaN values are not represented separately by this interval type.
    interval Mul(const interval& x, const interval& y) const override;
    void     testMul();
    //
    // Public API: enclose floating division with directly divided endpoints. Float/double
    // bounds round outward; quad/fixed retain their previous evaluation.
    // Empty operands yield empty; zero or an indeterminate infinite corner gives
    // [-inf, +inf]. The LSB estimate remains separate from numeric bound inclusion.
    interval Div(const interval& x, const interval& y) const override;
    void     testDiv();
    //
    // Public API: reciprocal bounds and their LSB estimate; empty stays empty.
    // Zero gives both signed infinities in float (signed zero is not tracked); double
    // retains its historical +inf convention. Float/double reciprocal endpoints
    // round outward; LSB remains an estimate.
    interval Inv(const interval& x) const override;
    void     testInv();
    //
    // Public API: enclose wrapping int32 negation or floating sign reversal.
    // Empty stays empty. Normalize widened integer corridors before endpoint work;
    // the result retains the input nature, including integral-valued floating zero.
    interval Neg(const interval& x) const override;
    void     testNeg();
    //
    interval Mod(const interval& x, double m) const;
    // Public API: integer C modulo for integer operands, otherwise numeric fmod.
    // Float/double floating bounds avoid rounded quotient tests and unsafe int casts;
    // invalid-only domains yield empty and NaN is not represented separately.
    interval Mod(const interval& x, const interval& y) const override;
    // Public API: floating modulo; in float/double, even integer-looking inputs
    // follow fmod rather than C integer modulo. NaN is not tracked separately.
    interval Fmod(const interval& x, const interval& y) const override;
    void     testMod();
    //

    interval Acos(const interval& x) const override;
    // Public API: numeric acos image on its valid domain. In float/double,
    // the native kernel encloses endpoints and extrema; quad/fixed retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval AcosBounds(const interval& x) const;
    void     testAcos();
    //
    interval Acosh(const interval& x) const override;
    // Public API: numeric acosh image on its valid domain. In float/double,
    // the native kernel encloses endpoints and extrema; quad/fixed retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval AcoshBounds(const interval& x) const;
    void     testAcosh();
    //
    // Public API: enclose int32 bitwise operations after integer conversion.
    // Normalize widened integers before bit analysis; empty stays empty.
    // Floating inputs are truncated and nonempty results remain integer.
    interval And(const interval& x, const interval& y) const override;
    void     testAnd();
    //
    interval Asin(const interval& x) const override;
    // Public API: numeric asin image on its valid domain. In float/double,
    // the native kernel encloses endpoints and extrema; quad/fixed retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval AsinBounds(const interval& x) const;
    void     testAsin();
    //
    interval Asinh(const interval& x) const override;
    // Public API: numeric asinh image on its valid domain. In float/double,
    // the native kernel encloses endpoints and extrema; quad/fixed retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval AsinhBounds(const interval& x) const;
    void     testAsinh();
    //
    interval Atan(const interval& x) const override;
    // Public API: numeric atan image on its valid domain. In float/double,
    // the native kernel encloses endpoints and extrema; quad/fixed retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval AtanBounds(const interval& x) const;
    void     testAtan();
    //
    interval Atan2(const interval& x, const interval& y) const override;
    // Public API: numeric atan2(y, x) image. Float/double endpoints round outward and a
    // possible negative-axis cut covers both signs of pi. Quad/fixed retain
    // their historical rule; LSB is an estimate and NaN is not tracked separately.
    interval Atan2Bounds(const interval& y, const interval& x) const;
    void     testAtan2();
    //
    interval Atanh(const interval& x) const override;
    // Public API: numeric atanh image on its valid domain. In float/double,
    // the native kernel encloses endpoints and extrema; quad/fixed retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval AtanhBounds(const interval& x) const;
    void     testAtanh();
    //
    // Public API: numeric ceil image; empty stays empty, and integer inputs
    // first convert to the target float. Results keep floating nature in every precision,
    // so subsequent operations cannot incorrectly enter the int32 wrapping branch.
    interval Ceil(const interval& x) const override;
    void     testCeil();
    interval Cos(const interval& x) const override;
    // Public API: numeric cos image on its valid domain. In float/double,
    // the native kernel encloses endpoints and extrema; quad/fixed retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval CosBounds(const interval& x) const;
    void     testCos();
    interval Cosh(const interval& x) const override;
    // Public API: numeric cosh image on its valid domain. In float/double,
    // the native kernel encloses endpoints and extrema; quad/fixed retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval CoshBounds(const interval& x) const;
    void     testCosh();
    interval Delay(const interval& x, const interval& y) const override;
    void     testDelay();
    // Public API: enclose comparison after target operand conversions.
    // Normalize widened integers before deciding a predicate.
    // Empty stays empty; a nonempty result is integer zero, one or both.
    interval Eq(const interval& x, const interval& y) const override;
    void     testEq();
    interval Exp(const interval& x) const override;
    // Public API: numeric exp image on its valid domain. In float/double,
    // the native kernel encloses endpoints and extrema; quad/fixed retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval ExpBounds(const interval& x) const;
    void     testExp();
    interval Exp10(const interval& x) const override;
    // Public API: enclose conversion to the program's floating type. In float mode,
    // convert the source endpoints once; IEEE narrowing is monotone. This is a known
    // conversion, not a bound calculation with an unknown analyzer rounding residual.
    // Normalize widened integer sources before conversion, retaining a negative
    // floating LSB in every precision. NaN keeps the historical empty convention.
    interval FloatCast(const interval& x) const override;
    void     testFloatCast();
    // Public API: numeric floor image; empty stays empty, and integer inputs
    // first convert to the target float. Results keep floating nature in every precision,
    // so subsequent operations cannot incorrectly enter the int32 wrapping branch.
    interval Floor(const interval& x) const override;
    void     testFloor();
    // Public API: enclose comparison after target operand conversions.
    // Normalize widened integers before deciding a predicate.
    // Empty stays empty; a nonempty result is integer zero, one or both.
    interval Ge(const interval& x, const interval& y) const override;
    void     testGe();
    // Public API: enclose comparison after target operand conversions.
    // Normalize widened integers before deciding a predicate.
    // Empty stays empty; a nonempty result is integer zero, one or both.
    interval Gt(const interval& x, const interval& y) const override;
    void     testGt();
    interval IntCast(const interval& x) const override;
    void     testIntCast();
    // Public API: enclose comparison after target operand conversions.
    // Normalize widened integers before deciding a predicate.
    // Empty stays empty; a nonempty result is integer zero, one or both.
    interval Le(const interval& x, const interval& y) const override;
    void     testLe();
    interval Log(const interval& x) const override;
    // Public API: numeric log image on its valid domain. In float/double,
    // the native kernel encloses endpoints and extrema; quad/fixed retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval LogBounds(const interval& x) const;
    void     testLog();
    interval Log10(const interval& x) const override;
    // Public API: numeric log10 image on its valid domain. In float/double,
    // the native kernel encloses endpoints and extrema; quad/fixed retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval Log10Bounds(const interval& x) const;
    void     testLog10();
    // Public API: enclose wrapping int32 left shift for counts in [0,31].
    // Convert operands to integers first. Invalid counts have no portable execution
    // contract and conservatively return full int32; this does not define such shifts.
    interval Lsh(const interval& x, const interval& y) const override;
    void     testLsh();
    // Public API: enclose comparison after target operand conversions.
    // Normalize widened integers before deciding a predicate.
    // Empty stays empty; a nonempty result is integer zero, one or both.
    interval Lt(const interval& x, const interval& y) const override;
    void     testLt();
    // Public API: enclose maximum after target operand conversions.
    // Empty stays empty; widened integers recover their signed hull.
    // The result is floating exactly when either operand is floating.
    interval Max(const interval& x, const interval& y) const override;
    void     testMax();
    interval Mem(const interval& x) const override;
    void     testMem();
    // Public API: enclose minimum after target operand conversions.
    // Empty stays empty; widened integers recover their signed hull.
    // The result is floating exactly when either operand is floating.
    interval Min(const interval& x, const interval& y) const override;
    void     testMin();
    // Public API: enclose comparison after target operand conversions.
    // Normalize widened integers before deciding a predicate.
    // Empty stays empty; a nonempty result is integer zero, one or both.
    interval Ne(const interval& x, const interval& y) const override;
    void     testNe();
    // Public API: enclose int32 bitwise complement after integer conversion.
    // ~x = -1-x reverses signed order, so endpoints suffice even for full int32;
    // enumerating that range would not terminate. The low bit is always significant.
    interval Not(const interval& x) const override;
    void     testNot();
    // Public API: enclose int32 bitwise operations after integer conversion.
    // Normalize widened integers before bit analysis; empty stays empty.
    // Floating inputs are truncated and nonempty results remain integer.
    interval Or(const interval& x, const interval& y) const override;
    void     testOr();
    interval Pow(const interval& x, const interval& y) const override;  // for all cases
    // Public API: numeric power image. Float/double floating powers round outward; negative
    // bases use integer exponents by parity. The existing nonnegative integer-power
    // path retains wrapping and LSB estimates. NaN is not tracked separately.
    interval PowBounds(const interval& x, const interval& y) const;
    void     testPow();
    // Public API: numeric IEEE remainder image. Float/double half-divisor bounds round
    // outward, including subnormals; invalid-only domains yield empty. NaN is not
    // tracked separately and LSB remains an estimate.
    interval Remainder(const interval& x, const interval& y) const override;
    void     testRemainder();
    // Public API: numeric rint image; empty stays empty, and integer inputs
    // first convert to the target float. Results keep floating nature in every precision,
    // so subsequent operations cannot incorrectly enter the int32 wrapping branch.
    interval Rint(const interval& x) const override;
    void     testRint();
    // Public API: numeric round image; empty stays empty, and integer inputs
    // first convert to the target float. Results keep floating nature in every precision,
    // so subsequent operations cannot incorrectly enter the int32 wrapping branch.
    interval Round(const interval& x) const override;
    void     testRound();
    // Public API: enclose arithmetic int32 right shift for counts in [0,31].
    // Arithmetic shifts round negative quotients downward; real scaling alone can
    // miss the actual integer result. Invalid counts conservatively return full int32.
    interval ARsh(const interval& x, const interval& y) const override;
    // Public API: enclose logical int32 right shift for counts in [0,31].
    // A positive shift sees the uint32 bit pattern; shift zero retains the original
    // signed value. Split at zero to handle the unsigned ordering discontinuity.
    interval LRsh(const interval& x, const interval& y) const override;
    void     testRsh();
    interval Sin(const interval& x) const override;
    // Public API: numeric sin image on its valid domain. In float/double,
    // the native kernel encloses endpoints and extrema; quad/fixed retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval SinBounds(const interval& x) const;
    void     testSin();
    interval Sinh(const interval& x) const override;
    // Public API: numeric sinh image on its valid domain. In float/double,
    // the native kernel encloses endpoints and extrema; quad/fixed retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval SinhBounds(const interval& x) const;
    void     testSinh();
    // Public API: square-root bounds over x intersected with [0, +inf], with an LSB
    // estimate. An empty domain yields empty; a zero-only domain yields exact zero.
    // Float/double bounds round outward independently of the host libm and rounding mode.
    interval Sqrt(const interval& x) const override;
    void     testSqrt();
    interval Tan(const interval& x) const override;
    // Public API: numeric tan image on its valid domain. In float/double,
    // the native kernel encloses endpoints and extrema; quad/fixed retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval TanBounds(const interval& x) const;
    void     testTan();
    interval Tanh(const interval& x) const override;
    // Public API: numeric tanh image on its valid domain. In float/double,
    // the native kernel encloses endpoints and extrema; quad/fixed retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval TanhBounds(const interval& x) const;
    void     testTanh();
    // Public API: enclose int32 bitwise operations after integer conversion.
    // Normalize widened integers before bit analysis; empty stays empty.
    // Floating inputs are truncated and nonempty results remain integer.
    interval Xor(const interval& x, const interval& y) const override;
    void     testXor();

    void testAll();
};
}  // namespace itv
