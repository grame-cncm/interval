#include <algorithm>
#include <functional>
#include <random>

#include "check.hh"
#include "interval.hh"
#include "interval_algebra.hh"

//------------------------------------------------------------------------------------------
// Interval Asinh
// interval Asinh(const interval& x) const;
// void testAsinh() const;
static interval domain(-HUGE_VAL, HUGE_VAL);

interval interval_algebra::Asinh(const interval& x) const
{
    return {asinh(x.lo()), asinh(x.hi())};
}

void interval_algebra::testAsinh() const
{
    analyzeumth(10, 1000, "asinh", interval(-10, 10), asinh, &interval_algebra::Asinh);
}
