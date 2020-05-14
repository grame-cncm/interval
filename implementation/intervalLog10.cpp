#include <algorithm>
#include <functional>
#include <random>

#include "check.hh"
#include "interval.hh"
#include "interval_algebra.hh"

//------------------------------------------------------------------------------------------
// Interval Log10
// interval Log10(const interval& x) const;
// void testLog10() const;

interval interval_algebra::Log10(const interval& x) const
{
    if (x.isEmpty()) {
        return {};
    } else {
        interval i = intersection(interval(0, HUGE_VAL), x);
        return {log10(i.lo()), log10(i.hi())};
    }
}

void interval_algebra::testLog10() const
{
    check("test algebra Log", Log10(interval(1, 10)), interval(log10(1), log10(10)));
    check("test algebra Log", Log10(interval(0, 10)), interval(log10(0), log10(10)));
    check("test algebra Log", Log10(interval(-10, 10)), interval(log10(0), log10(10)));
}
