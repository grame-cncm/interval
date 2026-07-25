#include <cmath>
#include <limits>
#include <sstream>

#include "interval/check.hh"
#include "interval/interval_algebra.hh"
#include "interval/affine_ops.hh"
#include "interval/affint.hh"
#include "interval/interval_def.hh"

using namespace itv;

int main()
{
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
    checkExact("fixed integer exponent preserves sufficient precision",
               algebra.Pow(interval(-0.375, -0.34375, -5), interval(2, 2, 0)),
               interval(0.1181640625, 0.140625, -6));


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
    }

    return reportCheckResults();
}
