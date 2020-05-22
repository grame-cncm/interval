#include <algorithm>
#include <functional>
#include <random>

#include "check.hh"
#include "interval.hh"
#include "interval_algebra.hh"

//------------------------------------------------------------------------------------------
// Interval Sqrt
// interval Sqrt(const interval& x) const;
// void testSqrt() const;

interval interval_algebra::Sqrt(const interval& x) const
{
    if (x.isEmpty()) return x;
    if (x.lo() < 0)
        return {};  // sqrt of negative numbers
    else
        return {sqrt(x.lo()), sqrt(x.hi())};
}

void interval_algebra::testSqrt() const
{
    analyzeumth(10, 1000, "sqrt", interval(0, 10), sqrt, &interval_algebra::Sqrt);
}
