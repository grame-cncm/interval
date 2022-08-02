#include <algorithm>
#include <functional>
#include <random>

#include "check.hh"
#include "interval.hh"
#include "interval_algebra.hh"

namespace itv {
//------------------------------------------------------------------------------------------
// Interval Acos
// interval Acos(const interval& x) const;
// void testAcos() const;

interval interval_algebra::Abs(const interval& x) const
{
    if (x.lo() >= 0) return x;
    if (x.hi() <= 0) return {-x.hi(), -x.lo()};
    return {0, std::max(std::abs(x.lo()), std::abs(x.hi()))};
}

void interval_algebra::testAbs() const
{
    analyzeUnaryMethod(10, 1000, "abs", interval(-10, 10), std::abs, &interval_algebra::Abs);
}
}  // namespace itv