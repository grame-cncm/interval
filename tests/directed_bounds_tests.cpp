#include <cfenv>
#include <cmath>
#include <limits>

#include "interval/check.hh"
#include "interval/interval_algebra.hh"
#include "interval/affine_ops.hh"
#include "interval/directed_rounding.hh"

using namespace itv;

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
    // An integer recursion widened to +inf is any int32 under the wrapping semantics :
    // the integer rules must not take the floating path (no wrap) on it.
    const interval widenedInt(0, HUGE_VAL, 0);
    check("double: integer + widened integer keeps the int32 wrap", true,
          algebra.Add(algebra.IntNum(12345), widenedInt).has(-1));
    check("double: integer * widened integer keeps the int32 wrap", true,
          algebra.Mul(algebra.IntNum(1103515245), widenedInt).has(-1));
    check("double: integer - widened integer stays within int32", true,
          algebra.Sub(algebra.IntNum(0), widenedInt) == interval(INT_MIN, INT_MAX, 0));
    check("double: widened integer % 100 is an integer remainder", true,
          algebra.Mod(widenedInt, algebra.IntNum(100)) == interval(-99, 99, 0));
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

    // A floating product residual may underflow to zero even when the product
    // is inexact. Exact significand comparisons must retain the neighbouring bit.
    check("double: positive inexact subnormal product is bracketed", true,
          algebra.Mul(interval(tiny, tiny, -1074), interval(0.75)) == interval(0, tiny));
    check("double: negative inexact subnormal product is bracketed", true,
          algebra.Mul(interval(-tiny, -tiny, -1074), interval(0.75)) == interval(-tiny, 0));
    check("double: an exact subnormal product stays a singleton", true,
          algebra.Mul(interval(2 * tiny, 2 * tiny, -1074), interval(0.5)).is(tiny));
    check("double: an exact subnormal quotient stays a singleton", true,
          algebra.Div(interval(2 * tiny, 2 * tiny, -1074), interval(2, 2, -24)).is(tiny));
    check("double: negative overflow keeps max-finite above -inf", true,
          algebra.Mul(interval(-largest, -largest, -24), interval(2, 2, -24))
              == interval(-HUGE_VAL, -largest));
    check("double: exact transcendental identities stay points", true,
          algebra.SinBounds(interval(0)).is(0) && algebra.CosBounds(interval(0)).is(1) &&
          algebra.ExpBounds(interval(0)).is(1) && algebra.LogBounds(interval(1)).is(0));
    check("double: large singleton trig inputs use an explicit conservative fallback", true,
          algebra.SinBounds(interval(1e100)) == interval(-1, 1) &&
          algebra.CosBounds(interval(1e100)) == interval(-1, 1));
    check("double: fmod with a huge quotient is exactly representable", true,
          algebra.Fmod(interval(largest), interval(3 * tiny)).is(2 * tiny));
    using namespace itv::detail;
    check("double: scalar IEEE remainder uses ties-to-even", true,
          directedBinary(BinaryOp::Remainder, 3, 2, Direction::Down) == -1 &&
          directedBinary(BinaryOp::Remainder, 5, 2, Direction::Up) == 1 &&
          directedBinary(BinaryOp::Remainder, -3, 2, Direction::Up) == 1);

    check("double: known positive images survive tiny reference calculations", true,
          algebra.AtanBounds(interval(0, tiny)).lo() >= 0 &&
          algebra.AsinBounds(interval(0, tiny)).lo() >= 0 &&
          algebra.LogBounds(interval(1, nextOne)).lo() >= 0);

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

    // A rational chord must contain the correctly-rounded t/10 values even
    // when coefficient evaluation introduces a second rounding.
    double lower0, lower1, upper0, upper1;
    achord(0, 1, 10, lower0, lower1);
    achord(0, 1, 10, upper0, upper1, true);
    AffItv chord{lower0, lower1, upper0, upper1, -24};
    bool chordContains = true;
    for (int t = 0; t <= 10; ++t)
        chordContains = chordContains && interval(chord.lo(t), chord.hi(t)).has(double(t) / 10);
    check("double: affine chord encloses rational interpolation", true, chordContains);
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

#ifndef __EMSCRIPTEN__
    // The kernel never changes a caller's rounding mode. The quotient endpoints
    // still enclose exact 1/10 with any supported IEEE rounding direction.
    const int saved = std::fegetround();
    for (int rounding : {FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO, FE_TONEAREST}) {
        std::fesetround(rounding);
        const interval result = algebra.Div(interval(1, 1, -24), interval(10, 10, -24));
        check("double: inclusion and host rounding mode are preserved", true,
              result.lo() <= std::nextafter(0.1, 0.0) && result.hi() >= 0.1 &&
              std::fegetround() == rounding);
    }
    std::fesetround(saved);
#endif
    libmCompensation() = true;
    return reportCheckResults();
}
