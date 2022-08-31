#include <algorithm>
#include <functional>
#include <random>

#include "check.hh"
#include "interval.hh"
#include "interval_algebra.hh"

namespace itv {
//------------------------------------------------------------------------------------------
// Interval Not
// interval Not(const interval& x) const;
// void testNot() const;

interval interval_algebra::Not(const interval& x) const
{
    if (x.isEmpty()) return x;
    int x0 = x.lo();
    int x1 = x.hi();

    int z0 = INT32_MAX;
    int z1 = INT32_MIN;

    for (int i = x0; i <= x1; i++) {
        int z = ~i;
        if (z < z0) z0 = z;
        if (z > z1) z1 = z;
    }
    return interval(double(z0), double(z1));
}

static double myNot(double x)
{
    int a = int(x);
    int b = ~a;
    return double(b);
}

void interval_algebra::testNot() const
{
    analyzeUnaryMethod(10, 1000, "not", interval(-10, -1), myNot, &interval_algebra::Not);
    analyzeUnaryMethod(10, 1000, "not", interval(10, 12), myNot, &interval_algebra::Not);
}
}  // namespace itv
