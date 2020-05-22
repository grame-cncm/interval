#include <algorithm>
#include <functional>
#include <random>

#include "check.hh"
#include "interval.hh"
#include "interval_algebra.hh"

//------------------------------------------------------------------------------------------
// Interval Tanh
// interval Tanh(const interval& x) const;
// void testTanh() const;

interval interval_algebra::Tanh(const interval& x) const
{
    if (x.isEmpty()) return {};
    return {tanh(x.lo()), tanh(x.hi())};
}

void interval_algebra::testTanh() const
{
    analyzeumth(20, 2000, "tanh", interval(-10 * M_PI, 10 * M_PI), tanh, &interval_algebra::Tanh);
}
