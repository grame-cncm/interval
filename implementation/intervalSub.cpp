#include <algorithm>
#include <functional>
#include <random>

#include "check.hh"
#include "interval.hh"
#include "interval_algebra.hh"

namespace itv {
//------------------------------------------------------------------------------------------
// Interval substraction

interval interval_algebra::Sub(const interval& x, const interval& y) const
{
    if (x.isEmpty() || y.isEmpty()) {
        return {};
    } else {
        return {x.lo() - y.hi(), x.hi() - y.lo()};
    }
}

void interval_algebra::testSub() const
{
    check("test algebra Sub", Sub(interval(0, 100), interval(10, 500)), interval(-500, 90));
}
}  // namespace itv
