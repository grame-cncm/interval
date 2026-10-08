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
    // Public API: inject the known binary64 literal as a point. In double mode retain
    // floating nature even for integral-looking literals; other modes keep legacy
    // injection. A known literal needs no widening for analyzer rounding error.
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

    interval Abs(const interval& x) const override;
    void     testAbs();
    //
    // Public API: enclose addition; double floating bounds round outward at each
    // endpoint, while the existing int32 wrapping and other precision paths remain.
    interval Add(const interval& x, const interval& y) const override;
    void     testAdd();
    //
    // Public API: enclose subtraction; double floating bounds round outward at each
    // endpoint, while the existing int32 wrapping and other precision paths remain.
    interval Sub(const interval& x, const interval& y) const override;
    void     testSub();
    //
    // Public API: enclose multiplication; double floating corners round outward.
    // Zero times an unbounded endpoint keeps the historical numeric-hull convention;
    // possible NaN values are not represented separately by this interval type.
    interval Mul(const interval& x, const interval& y) const override;
    void     testMul();
    //
    // Public API: enclose floating division with directly divided endpoints. Double
    // bounds round outward; other precisions retain their previous evaluation.
    // Empty operands yield empty; zero or an indeterminate infinite corner gives
    // [-inf, +inf]. The LSB estimate remains separate from numeric bound inclusion.
    interval Div(const interval& x, const interval& y) const override;
    void     testDiv();
    //
    // Public API: reciprocal bounds and their LSB estimate; empty stays empty.
    // The exact zero point yields +inf with default floating LSB, preserving the
    // historical unsigned-zero convention without estimating a precision at zero.
    // Double reciprocal endpoints round outward; LSB remains an estimate.
    interval Inv(const interval& x) const override;
    void     testInv();
    //
    interval Neg(const interval& x) const override;
    void     testNeg();
    //
    interval Mod(const interval& x, double m) const;
    // Public API: integer C modulo for integer operands, otherwise numeric fmod.
    // Double floating bounds avoid rounded quotient tests and unsafe int casts;
    // invalid-only domains yield empty and NaN is not represented separately.
    interval Mod(const interval& x, const interval& y) const override;
    // Public API: floating modulo; in double precision, even integer-looking inputs
    // follow fmod rather than C integer modulo. NaN is not tracked separately.
    interval Fmod(const interval& x, const interval& y) const override;
    void     testMod();
    //

    interval Acos(const interval& x) const override;
    // Public API: numeric acos image on its valid domain. In double precision,
    // the native kernel encloses endpoints and extrema; other precisions retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval AcosBounds(const interval& x) const;
    void     testAcos();
    //
    interval Acosh(const interval& x) const override;
    // Public API: numeric acosh image on its valid domain. In double precision,
    // the native kernel encloses endpoints and extrema; other precisions retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval AcoshBounds(const interval& x) const;
    void     testAcosh();
    //
    interval And(const interval& x, const interval& y) const override;
    void     testAnd();
    //
    interval Asin(const interval& x) const override;
    // Public API: numeric asin image on its valid domain. In double precision,
    // the native kernel encloses endpoints and extrema; other precisions retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval AsinBounds(const interval& x) const;
    void     testAsin();
    //
    interval Asinh(const interval& x) const override;
    // Public API: numeric asinh image on its valid domain. In double precision,
    // the native kernel encloses endpoints and extrema; other precisions retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval AsinhBounds(const interval& x) const;
    void     testAsinh();
    //
    interval Atan(const interval& x) const override;
    // Public API: numeric atan image on its valid domain. In double precision,
    // the native kernel encloses endpoints and extrema; other precisions retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval AtanBounds(const interval& x) const;
    void     testAtan();
    //
    interval Atan2(const interval& x, const interval& y) const override;
    // Public API: numeric atan2(y, x) image. Double endpoints round outward and a
    // possible negative-axis cut covers both signs of pi. Other precisions retain
    // their historical rule; LSB is an estimate and NaN is not tracked separately.
    interval Atan2Bounds(const interval& y, const interval& x) const;
    void     testAtan2();
    //
    interval Atanh(const interval& x) const override;
    // Public API: numeric atanh image on its valid domain. In double precision,
    // the native kernel encloses endpoints and extrema; other precisions retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval AtanhBounds(const interval& x) const;
    void     testAtanh();
    //
    interval Ceil(const interval& x) const override;
    void     testCeil();
    interval Cos(const interval& x) const override;
    // Public API: numeric cos image on its valid domain. In double precision,
    // the native kernel encloses endpoints and extrema; other precisions retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval CosBounds(const interval& x) const;
    void     testCos();
    interval Cosh(const interval& x) const override;
    // Public API: numeric cosh image on its valid domain. In double precision,
    // the native kernel encloses endpoints and extrema; other precisions retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval CoshBounds(const interval& x) const;
    void     testCosh();
    interval Delay(const interval& x, const interval& y) const override;
    void     testDelay();
    interval Eq(const interval& x, const interval& y) const override;
    void     testEq();
    interval Exp(const interval& x) const override;
    // Public API: numeric exp image on its valid domain. In double precision,
    // the native kernel encloses endpoints and extrema; other precisions retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval ExpBounds(const interval& x) const;
    void     testExp();
    interval Exp10(const interval& x) const override;
    interval FloatCast(const interval& x) const override;
    void     testFloatCast();
    interval Floor(const interval& x) const override;
    void     testFloor();
    interval Ge(const interval& x, const interval& y) const override;
    void     testGe();
    interval Gt(const interval& x, const interval& y) const override;
    void     testGt();
    interval IntCast(const interval& x) const override;
    void     testIntCast();
    interval Le(const interval& x, const interval& y) const override;
    void     testLe();
    interval Log(const interval& x) const override;
    // Public API: numeric log image on its valid domain. In double precision,
    // the native kernel encloses endpoints and extrema; other precisions retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval LogBounds(const interval& x) const;
    void     testLog();
    interval Log10(const interval& x) const override;
    // Public API: numeric log10 image on its valid domain. In double precision,
    // the native kernel encloses endpoints and extrema; other precisions retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval Log10Bounds(const interval& x) const;
    void     testLog10();
    interval Lsh(const interval& x, const interval& y) const override;
    void     testLsh();
    interval Lt(const interval& x, const interval& y) const override;
    void     testLt();
    interval Max(const interval& x, const interval& y) const override;
    void     testMax();
    interval Mem(const interval& x) const override;
    void     testMem();
    interval Min(const interval& x, const interval& y) const override;
    void     testMin();
    interval Ne(const interval& x, const interval& y) const override;
    void     testNe();
    interval Not(const interval& x) const override;
    void     testNot();
    interval Or(const interval& x, const interval& y) const override;
    void     testOr();
    interval Pow(const interval& x, const interval& y) const override;  // for all cases
    // Public API: numeric power image. Double floating powers round outward; negative
    // bases use integer exponents by parity. The existing nonnegative integer-power
    // path retains wrapping and LSB estimates. NaN is not tracked separately.
    interval PowBounds(const interval& x, const interval& y) const;
    void     testPow();
    // Public API: numeric IEEE remainder image. Double half-divisor bounds round
    // outward, including subnormals; invalid-only domains yield empty. NaN is not
    // tracked separately and LSB remains an estimate.
    interval Remainder(const interval& x, const interval& y) const override;
    void     testRemainder();
    interval Rint(const interval& x) const override;
    void     testRint();
    interval Round(const interval& x) const override;
    void     testRound();
    interval ARsh(const interval& x, const interval& y) const override;
    interval LRsh(const interval& x, const interval& y) const override;
    void     testRsh();
    interval Sin(const interval& x) const override;
    // Public API: numeric sin image on its valid domain. In double precision,
    // the native kernel encloses endpoints and extrema; other precisions retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval SinBounds(const interval& x) const;
    void     testSin();
    interval Sinh(const interval& x) const override;
    // Public API: numeric sinh image on its valid domain. In double precision,
    // the native kernel encloses endpoints and extrema; other precisions retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval SinhBounds(const interval& x) const;
    void     testSinh();
    // Public API: square-root bounds over x intersected with [0, +inf], with an LSB
    // estimate. An empty domain yields empty; a zero-only domain yields exact zero.
    // Double bounds round outward independently of the host libm and rounding mode.
    interval Sqrt(const interval& x) const override;
    void     testSqrt();
    interval Tan(const interval& x) const override;
    // Public API: numeric tan image on its valid domain. In double precision,
    // the native kernel encloses endpoints and extrema; other precisions retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval TanBounds(const interval& x) const;
    void     testTan();
    interval Tanh(const interval& x) const override;
    // Public API: numeric tanh image on its valid domain. In double precision,
    // the native kernel encloses endpoints and extrema; other precisions retain
    // their historical rule. LSB is an estimate; NaN is not tracked separately.
    interval TanhBounds(const interval& x) const;
    void     testTanh();
    interval Xor(const interval& x, const interval& y) const override;
    void     testXor();

    void testAll();
};
}  // namespace itv
