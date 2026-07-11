#include <cmath>
#include <limits>
#include <sstream>

#include "interval/check.hh"
#include "interval/interval_algebra.hh"
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

    return reportCheckResults();
}
