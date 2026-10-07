#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <sstream>

#include "interval/check.hh"
#include "interval/interval_algebra.hh"
#include "interval/affine_ops.hh"
#include "interval/affint.hh"
#include "interval/interval_def.hh"

using namespace itv;

int main()
{
    // the defaults, read before any test changes them
    check("the program is in double precision by default", true, programPrecision() == 2);
    check("the libm is compensated by default", true, libmCompensation());

    const interval all = fullFinite();
    const interval nil = empty();

    checkExact("default interval equals fullFinite", interval(), all);
    check("fullFinite is valid", true, all.isValid());
    check("fullFinite is not empty", false, all.isEmpty());
    check("fullFinite contains zero", true, all.hasZero());
    check("fullFinite contains the largest finite double", true,
          all.has(std::numeric_limits<double>::max()));
    check("fullFinite excludes positive infinity", false, all.has(HUGE_VAL));

    check("empty is empty", true, nil.isEmpty());
    check("empty is invalid", false, nil.isValid());
    checkExact("empty intervals compare equal", empty(), empty());
    check("empty inclusion is reflexive", true, nil <= nil);
    check("empty is included in every interval", true, nil <= all);
    check("nonempty is not included in empty", false, all <= nil);
    check("empty is strictly below nonempty", true, nil < all);

    std::ostringstream printedEmpty;
    printedEmpty << nil;
    check("empty()", nil);
    check("empty has an unambiguous representation", true, printedEmpty.str() == "empty()");

    const interval x(-2, 5, -8);
    checkExact("intersection with fullFinite", intersection(x, all), interval(-2, 5, -24));
    checkExact("union with empty", reunion(x, nil), x);
    checkExact("intersection with empty", intersection(x, nil), nil);
    checkExact("union with fullFinite", reunion(x, all), all);
    checkExact("disjoint intersection", intersection(interval(-2, -1), interval(1, 2)), empty());
    checkExact("overlapping intersection", intersection(interval(-2, 3), interval(1, 5)),
               interval(1, 3));
    checkExact("overlapping union", reunion(interval(-2, 3), interval(1, 5)), interval(-2, 5));

    interval_algebra algebra;
    // Direct division must enclose the machine quotient, not a product with an
    // independently rounded reciprocal. Exercise points, signs and interior values.
    check("double division contains 3/10", true,
          algebra.Div(interval(3.0), interval(10.0)).has(3.0 / 10.0));
    checkExact("division with an empty operand", algebra.Div(nil, x), nil);
    checkExact("division by an empty operand", algebra.Div(x, nil), nil);
    for (int precision : {1, 2}) {
        programPrecision() = precision;
        bool contains      = true;
        for (double sign : {-1.0, 1.0}) {
            const interval numerator(-3.7, 12.1), denominator(sign * 0.3, sign * 9.7);
            const interval quotient = algebra.Div(numerator, denominator);
            for (int i = 0; i <= 32; ++i) {
                for (int j = 0; j <= 32; ++j) {
                    // Clamp interpolation's rounding so the samples belong to
                    // the tested input intervals, including at their endpoints.
                    const double a     = std::clamp(numerator.lo() + numerator.size() * i / 32,
                                                    numerator.lo(), numerator.hi());
                    const double b     = std::clamp(denominator.lo() + denominator.size() * j / 32,
                                                    denominator.lo(), denominator.hi());
                    const double value = precision == 1 ? double(float(a) / float(b)) : a / b;
                    contains           = contains && quotient.has(value);
                }
            }
        }
        check(
            "division contains sampled machine quotients at precision " + std::to_string(precision),
            true, contains);
    }
    programPrecision() = 2;
    check("division across zero covers both infinities", true,
          algebra.Div(interval(1), interval(-1, 1)).has(-HUGE_VAL) &&
              algebra.Div(interval(1), interval(-1, 1)).has(HUGE_VAL));
    check("division at an indeterminate infinite corner stays conservative", true,
          algebra.Div(interval(1, HUGE_VAL), interval(1, HUGE_VAL)).has(1));
    checkExact("addition", algebra.Add(interval(-2, 3), interval(4, 5)), interval(2, 8));
    checkExact("subtraction", algebra.Sub(interval(-2, 3), interval(4, 5)), interval(-7, -1));
    checkExact("multiplication", algebra.Mul(interval(-2, 3), interval(4, 5)),
               interval(-10, 15, -48));
    checkExact("negation", algebra.Neg(interval(-2, 3)), interval(-3, 2));
    checkExact("absolute value", algebra.Abs(interval(-2, 3)), interval(0, 3));
    checkExact("minimum", algebra.Min(interval(-2, 3), interval(1, 5)), interval(-2, 3));
    checkExact("maximum", algebra.Max(interval(-2, 3), interval(1, 5)), interval(1, 5));
    checkExact("integer power without wrapping", algebra.Pow(interval(2, 3, 0), interval(2, 2, 0)),
               interval(4, 9, 2));
    checkExact("integer power with possible wrapping",
               algebra.Pow(interval(60, 62, 0), interval(8, 106, 0)),
               interval(static_cast<double>(INT32_MIN), static_cast<double>(INT32_MAX), 0));
    // the raw bounds of the algebra : the compensation of the libm is turned off here
    libmCompensation() = false;
    checkExact("fixed integer exponent preserves sufficient precision",
               algebra.Pow(interval(-0.375, -0.34375, -5), interval(2, 2, 0)),
               interval(0.1181640625, 0.140625, -6));
    libmCompensation() = true;


    // ---- affine-in-time intervals -----------------------------------------------------
    {
        const double    T = 1000;
        affine_algebra  aa(T);

        // the rate-0 subdomain embeds the ordinary intervals
        const AffItv c = fromItv(interval(-2, 3, -8));
        checkExact("affine: rate-0 round-trip", toItv(c, T), interval(-2, 3, -8));

        // the counter x = x@1 + 1 is STATIONARY at rate 1: 1 + delay(x) == x
        const AffItv x  = {0, 0, 9, 1, 0};  // lo(t) = 0, hi(t) = 9 + t
        const AffItv x2 = aa.Add(aa.IntNum(1), aa.Mem(x));
        check("affine: counter stationary (hi intercept)", true, x2.b0 == x.b0 && x2.b1 == x.b1);
        check("affine: counter stationary (lo)", true, x2.a0 >= 0 && x2.a1 == 0);

        // mod of a climbing corridor is a HORIZONTAL band -- the hull of the sawtooth
        const AffItv m = aa.Mod(x, fromItv(interval(2000, 2000, 0)));
        check("affine: mod kills the rate", true, m.isConst());
        check("affine: mod band is [0, 2000)", true, m.a0 >= 0 && m.b0 < 2000);

        // the INTEGER path of Mod (C semantics : sign of x, magnitude < |y|)
        itv::interval_algebra ia;
        checkExact("integer mod is [0, n-1]", ia.Mod(interval(1, 192000, 0), interval(2, 2, 0)),
                   interval(0, 1, 0));
        checkExact("integer mod follows x's sign",
                   ia.Mod(interval(-10, 10, 0), interval(3, 3, 0)), interval(-2, 2, 0));
        checkExact("integer mod keeps a small x",
                   ia.Mod(interval(3, 5, 0), interval(7, 7, 0)), interval(3, 5, 0));

        // integer Neg and Sub wrap around int32 like Add (computational semantics)
        const double IMIN = -2147483648.0, IMAX = 2147483647.0;
        checkExact("integer neg of INT_MIN wraps to INT_MIN",
                   ia.Neg(interval(IMIN, IMIN, 0)), interval(IMIN, IMIN, 0));
        checkExact("integer neg straddling the boundary widens to full int32",
                   ia.Neg(interval(IMIN, IMIN + 5, 0)), interval(IMIN, IMAX, 0));
        checkExact("integer neg without wrap stays exact",
                   ia.Neg(interval(-5, 10, 0)), interval(-10, 5, 0));
        checkExact("integer sub reaching -INT_MIN widens to full int32",
                   ia.Sub(interval(0, 0, 0), interval(IMIN, IMAX, 0)),
                   interval(IMIN, IMAX, 0));
        checkExact("integer sub without wrap stays exact",
                   ia.Sub(interval(0, 0, 0), interval(-10, 10, 0)), interval(-10, 10, 0));
        checkExact("integer sub wholly past the boundary wraps coherently",
                   ia.Sub(interval(IMIN, IMIN, 0), interval(1, 1, 0)),
                   interval(IMAX, IMAX, 0));

        // the general conformity pass : Abs, Lsh, LRsh, Remainder
        checkExact("integer abs containing INT_MIN reaches both ends",
                   ia.Abs(interval(IMIN, 5, 0)), interval(IMIN, IMAX, 0));
        checkExact("integer lsh overflowing widens to full int32",
                   ia.Lsh(interval(1, 1, 0), interval(31, 31, 0)),
                   interval(IMIN, IMAX, 0 + 31));
        checkExact("logical rsh of a negative by k>=1 is bounded",
                   ia.LRsh(interval(-8, -8, 0), interval(1, 1, 0)),
                   interval(0, IMAX, 0));
        checkExact("logical rsh by a possible 0 keeps the negatives reachable",
                   ia.LRsh(interval(-8, -8, 0), interval(0, 1, 0)),
                   interval(IMIN, IMAX, 0));
        checkExact("remainder band", ia.Remainder(interval(0, 100, 0), interval(2, 2, 0)),
                   interval(-1, 1, 0));
        checkExact("remainder identity below half the divisor",
                   ia.Remainder(interval(3, 3, 0), interval(8, 8, 0)),
                   interval(3, 3, 0));
        checkExact("int cast of a possibly out-of-range float is the full int32 range",
                   ia.IntCast(interval(1e10, 2e10, -24)), interval(IMIN, IMAX, 0));
        checkExact("int cast of an in-range float truncates exactly",
                   ia.IntCast(interval(-3.8, 4.9, -24)), interval(-3, 4, 0));

        // widening proposes the observed per-round rate, then escalates
        const AffItv w1 = awiden({0, 0, 8, 0, 0}, {0, 0, 9, 0, 0}, T);
        check("affine: widening proposes the rate", true, w1.b1 == 1 && w1.b0 == 9);
        const AffItv w2 = awiden({0, 0, 9, 1, 0}, {0, 0, 20, 1, 0}, T);
        check("affine: super-linear escalates to the int top", true,
              w2.b0 == 2147483647.0 && w2.b1 == 0);

        // the collapse caps an integer chain past its wrap date
        checkExact("affine: int32 cap at the bridge", toItv({0, 0, 10, 3e6, 0}, T),
                   interval(-2147483648.0, 2147483647.0, 0));

        // foreign entities are fullFinite (sound near-top), not empty
        check("affine: foreign is not neutral", false,
              aa.ForeignConst(0, aempty(), aempty()).isEmpty());

        // A reciprocal is curved, so endpoint interpolation excludes interior
        // values. Collapsing rated denominators also catches an interior pole.
        affine_algebra shortHorizon(10);
        const AffItv   numerator = fromItv(interval(1.0));
        for (double sign : {-1.0, 1.0}) {
            const AffItv denominator = {sign, sign, sign, sign, -24};
            const AffItv quotient    = shortHorizon.Div(numerator, denominator);
            bool         contains    = quotient.isConst();
            for (int t = 0; t <= 10; ++t) {
                const double value = 1.0 / (sign * (1 + t));
                contains           = contains && quotient.lo(t) <= value && value <= quotient.hi(t);
            }
            check("affine: rated denominator contains interior reciprocals", true, contains);
        }
        const AffItv pole = shortHorizon.Div(numerator, {-1, 0.2, -1, 0.2, -24});
        check("affine: interior zero crossing covers both infinities", true,
              pole.lo(5) == -HUGE_VAL && pole.hi(5) == HUGE_VAL);
        const AffItv linear = shortHorizon.Div({2, 2, 2, 2, -24}, fromItv(interval(2.0)));
        check("affine: constant denominator preserves the numerator rate", true,
              linear.a0 == 1 && linear.b0 == 1 && linear.a1 == 1 && linear.b1 == 1);
    }

    // ---- the precision of the program (programPrecision, programBound, ulpMargin) ----
    {
        interval_algebra algebra;
        const int        saved = programPrecision();

        // double precision : the bounds stay doubles, nothing is rounded
        check("double: 0.1 stays the double 0.1", true, interval(0.1).lo() == 0.1);
        check("double: a bound between two floats is kept", true,
              interval(0, 99.99999999).hi() == 99.99999999);

        // single precision : the bounds of a float-carried value are floats
        programPrecision() = 1;
        check("single: 0.1 is the float 0.1f", true, interval(0.1).lo() == double(0.1f));
        check("single: a constant stays a point", true,
              algebra.Add(interval(0.1), interval(0.2)).isconst());
        check("single: 0.1 + 0.2 is the float sum the program computes", true,
              algebra.Add(interval(0.1), interval(0.2)).lo() == double(0.1f + 0.2f));
        check("single: a bound between two floats rounds to the nearest float", true,
              interval(0, 99.99999999).hi() == 100.0);
        check("single: an integer interval (lsb >= 0) is not rounded", true,
              interval(0, 16777217.0, 0).hi() == 16777217.0);
        check("single: an integer bound beyond 2^24 is kept", true,
              interval(0, 16777217.0).hi() == 16777217.0);
        check("single: a tiny positive bound underflows to zero", true,
              interval(1e-300, 1).lo() == 0);
        checkExact("single: underflow normalizes a zero point", interval(1e-300, 1e-300, -149),
                   interval(0));
        // Consumers of an underflowed bound must handle zero before estimating
        // precision; derivatives at these singular points give NaN or infinity.
        const interval underflowed(1e-300);
        checkExact("single: log of underflowed zero", algebra.Log(underflowed),
                   interval(-HUGE_VAL, -HUGE_VAL, -24));
        checkExact("single: log10 of underflowed zero", algebra.Log10(underflowed),
                   interval(-HUGE_VAL, -HUGE_VAL, -24));
        checkExact("single: reciprocal of underflowed zero", algebra.Inv(underflowed),
                   interval(HUGE_VAL, HUGE_VAL, -24));
        checkExact("single: sqrt of underflowed zero", algebra.Sqrt(underflowed), interval(0));
        // Check both signs and the transition from a subnormal result to zero.
        // Volatile inputs ensure the reference executes at float precision.
        volatile float smallestNormal = std::numeric_limits<float>::min();
        volatile float factor         = 0.1f;
        for (double sign : {-1.0, 1.0}) {
            const interval product = algebra.Mul(
                interval(sign * double(smallestNormal), sign * double(smallestNormal), -149),
                interval(double(factor)));
            check("single: multiplication contains its subnormal result", true,
                  product.has(float(sign) * smallestNormal * factor));
        }
        volatile float smallestSubnormal = std::numeric_limits<float>::denorm_min();
        check("single: arithmetic underflow contains the machine zero", true,
              algebra.Mul(interval(double(smallestSubnormal)), interval(0.25))
                  .has(smallestSubnormal * 0.25f));
        check("single: division contains a subnormal quotient", true,
              algebra.Div(interval(double(smallestNormal)), interval(10.0))
                  .has(smallestNormal / 10.0f));
        check("single: division by a subnormal contains float overflow", true,
              algebra.Div(interval(1.0), interval(double(smallestSubnormal)))
                  .has(1.0f / smallestSubnormal));
        check("single: Pow keeps zero out of its positive domain", true,
              algebra.Pow(interval(0, 1), interval(0.5)).has(0.5));
        check("single: ulpMargin is k float ulps at the magnitude", true,
              ulpMargin(0, 100, 4) == 4 * 0x1p-23 * 100);

        // the float sums and products of the program stay in their intervals (fixed seed)
        {
            interval X(-3.7, 12.1), Y(0.003, 99.99);
            interval S = algebra.Add(X, Y), P = algebra.Mul(X, Y);
            std::mt19937                          gen(42);
            std::uniform_real_distribution<float> dx(float(X.lo()), float(X.hi()));
            std::uniform_real_distribution<float> dy(float(Y.lo()), float(Y.hi()));
            bool ok = S.has(float(X.hi()) + float(Y.hi())) && P.has(float(X.hi()) * float(Y.hi()));
            for (int i = 0; i < 100000 && ok; i++) {
                float a = dx(gen), b = dy(gen);
                ok      = S.has(a + b) && P.has(a * b);
            }
            check("single: the float sums and products stay in their intervals", true, ok);
        }

        // double precision : nothing to round, the bounds are already doubles
        programPrecision() = 2;
        check("double: 0.1 stays the double 0.1", true, interval(0.1).lo() == 0.1);
        check("double: ulpMargin is k double ulps at the magnitude", true,
              ulpMargin(0, 100, 4) == 4 * 0x1p-52 * 100);

        // the compensation of the libm : 2 ulps outward, within the image of the function,
        // never for a point, an exact 0 or an integer result
        {
            programPrecision() = 1;
            interval x(0, 1.00000596);
            interval s = algebra.Sin(x);
            libmCompensation() = false;
            interval raw = algebra.Sin(x);
            libmCompensation() = true;
            float  up2 = std::nextafter(std::nextafter(float(raw.hi()), 2.0f), 2.0f);
            check("libm: sin widens its upper bound by 2 float ulps", true, s.hi() == double(up2));
            check("libm: sin keeps its exact 0 bound", true, s.lo() == 0);
            check("libm: the image of sin caps the widening", true,
                  algebra.Sin(interval(0, 7)).hi() == 1 && algebra.Sin(interval(0, 7)).lo() == -1);
            check("libm: a constant stays a point", true, algebra.Sin(interval(0.5)).isconst());
            check("libm: exp stays >= 0", true, algebra.Exp(interval(-200, 0)).lo() >= 0);
            // the witness : 118.83905 * sinf(1.00000596) is 100.0 in float, one ulp above the
            // bound computed with sin in double ; the compensated bound covers it
            // volatile : a sinf of a constant would be folded by the C++ compiler, with a
            // correct rounding, instead of calling the libm
            volatile float vx   = 1.00000596f;
            float          prog = 118.83905f * std::sin(float(vx));
            interval p = algebra.Mul(interval(118.83905), algebra.Sin(interval(0, 1.00000596)));
            check("libm: the float program value stays in its interval", true, p.has(prog));
            programPrecision() = 2;
            libmCompensation() = false;
            check("libm: without compensation, the raw bounds", true,
                  algebra.Sin(interval(0, 1.00000596)).hi() == std::sin(1.00000596));
            libmCompensation() = true;
        }

        programPrecision() = saved;
    }

    return reportCheckResults();
}
