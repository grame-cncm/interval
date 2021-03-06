#include <algorithm>
#include <functional>
#include <random>

#include "check.hh"
#include "interval.hh"
#include "interval_algebra.hh"

//------------------------------------------------------------------------------------------
// Interval And
// interval And(const interval& x, const interval& y) const;
// void testAnd() const;

//----------------------booleans&bits--------------------------------------
inline int bitmask(double x)
{
    int v = int(x);
    for (int i = 1; i < 32; i *= 2) {
        v |= v >> i;
    }
    return v;
}

// 00111
// 01000  +1
// 00000  &
bool isBitMask(int x)
{
    return ((x + 1) & x) == 0;
}

// A&B <= B et A&B <= A
// inline interval operator&(const interval& x, const interval& y)
// {
//     if (x.valid && y.valid) {
//         if (x.lo >= 0 && y.lo >= 0) {
//             return interval(0, bitmask(x.hi) & bitmask(y.hi));
//         } else if (y.lo >= 0) {
//             return interval(0, bitmask(y.hi));
//         } else if (x.lo >= 0) {
//             return interval(0, bitmask(y.hi));
//         } else {
//             return interval();
//         }
//     } else if (x.valid && x.lo >= 0) {
//         return interval(0, bitmask(x.hi));
//     } else if (y.valid && y.lo >= 0) {
//         return interval(0, bitmask(y.hi));
//     } else {
//         return interval();
//     }
// }
static double myAnd(double x, double y)
{
    int a = int(x);
    int b = int(y);
    int c = a & b;
    return double(c);
}

// we assume mask = 2^n-1
// periodicity 2^n
// resulting interval <= [0..mask]

// 3.9 .. 10.8 => 6.9 ==> 6
// 3 .. 10 => 7

interval bmAnd(const interval& x, int mask)
{
    int lo = int(x.lo());
    int hi = int(x.hi());
    int w  = hi - lo;
    int p  = mask + 1;

    if (w >= p) {
        return interval{0, double(mask)};
    } else {
        // shit x to be positive
        lo = lo % p;
        if (lo < 0) lo += p;
        hi = lo + w;

        if (hi < p) {
            return interval{double(lo), double(hi)};
        } else {
            return interval{0, double(mask)};
        }
    }
}
/*
interval interval_algebra::And(const interval& x, const interval& y) const
{
    if (x.isEmpty() || y.isEmpty()) return {};
    if (x.isconst()) {
        return bmAnd(y, bitmask(x.hi()));
    } else if (y.isconst()) {
        return bmAnd(x, bitmask(y.hi()));
    } else {
        if (y.lo() > 0) {
            return bmAnd(x, bitmask(y.hi()));
        } else if (x.lo() > 0) {
            return bmAnd(y, bitmask(x.hi()));
        } else {
            // we don't know !
            return {-13, -11};
        }
    }
}
*/

// BRUTE FORCE
interval interval_algebra::And(const interval& x, const interval& y) const
{
    if (x.isEmpty() || y.isEmpty()) return {};
    int x0 = x.lo();
    int x1 = x.hi();
    int y0 = y.lo();
    int y1 = y.hi();

    int z0 = INT32_MAX;
    int z1 = INT32_MIN;

    for (int i = x0; i <= x1; i++) {
        for (int j = y0; j <= y1; j++) {
            int z = i & j;
            if (z < z0) z0 = z;
            if (z > z1) z1 = z;
        }
    }
    return interval(double(z0), double(z1));
}

void interval_algebra::testAnd() const
{
    analyzeBinaryMethod(10, 2000, "And", interval(-1000, -800), interval(127), myAnd, &interval_algebra::And);
    analyzeBinaryMethod(10, 2000, "And", interval(-1000, -800), interval(123), myAnd, &interval_algebra::And);
    analyzeBinaryMethod(10, 2000, "And", interval(-128, 128), interval(127), myAnd, &interval_algebra::And);
    analyzeBinaryMethod(10, 2000, "And", interval(0, 1000), interval(63, 127), myAnd, &interval_algebra::And);
    analyzeBinaryMethod(10, 2000, "And", interval(-1000, 1000), interval(63, 127), myAnd, &interval_algebra::And);
}
