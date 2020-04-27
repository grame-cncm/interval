/************************************************************************
 ************************************************************************
    FAUST compiler
    Copyright (C) 2003-2018 GRAME, Centre National de Creation Musicale
    ---------------------------------------------------------------------
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 ************************************************************************
 ************************************************************************/

#pragma once

#include <cmath>
#include <iostream>

// ***************************************************************************
//
//     An Interval is a (possibly empty) set of numbers approximated by two boundaries
//
//****************************************************************************

class interval {
   private:
    bool   fEmpty{true};  ///< true when the interval is empty
    double fLo{NAN};      ///< minimal value
    double fHi{NAN};      ///< maximal value

   public:
    //-------------------------------------------------------------------------
    // constructors
    //-------------------------------------------------------------------------

    interval() = default;

    interval(double n, double m)
    {
        if (std::isnan(n) || std::isnan(m)) {
            fEmpty = true;
            fLo    = NAN;
            fHi    = NAN;
        } else {
            fEmpty = false;
            fLo    = std::min(n, m);
            fHi    = std::max(n, m);
        }
    }

    explicit interval(double n) : interval(n, n) {}

    // interval(const interval& r) : fEmpty(r.empty()), fLo(r.lo()), fHi(r.hi()) {}

    //-------------------------------------------------------------------------
    // basic properties
    //-------------------------------------------------------------------------

    bool   empty() const { return fEmpty; }
    double lo() const { return fLo; }
    double hi() const { return fHi; }
};

//-------------------------------------------------------------------------
// printing
//-------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& dst, const interval& i)
{
    if (i.empty()) {
        return dst << "interval()";
    } else if (i.lo() == i.hi()) {
        return dst << "interval(" << i.lo() << ")";
    } else {
        return dst << "interval(" << i.lo() << "," << i.hi() << ")";
    }
}

//-------------------------------------------------------------------------
// set operations
//-------------------------------------------------------------------------

inline interval intersection(const interval& i, const interval& j)
{
    if (i.empty()) {
        return i;
    } else if (j.empty()) {
        return j;
    } else {
        double l = std::max(i.lo(), j.lo());
        double h = std::min(i.hi(), j.hi());
        if (l > h) {
            return {};
        } else {
            return {l, h};
        }
    }
}

inline interval reunion(const interval& i, const interval& j)
{
    if (i.empty()) {
        return j;
    } else if (j.empty()) {
        return i;
    } else {
        double l = std::min(i.lo(), j.lo());
        double h = std::max(i.hi(), j.hi());
        return {l, h};
    }
}

//-------------------------------------------------------------------------
// predicates
//-------------------------------------------------------------------------

// basic predicates
inline bool operator==(const interval& i, const interval& j)
{
    return (i.empty() && j.empty()) || ((i.lo() == j.lo()) && (i.hi() == j.hi()));
}

inline bool operator<=(const interval& i, const interval& j)
{
    return reunion(i, j) == j;
}

// additional predicates
inline bool operator!=(const interval& i, const interval& j)
{
    return !(i == j);
}

inline bool operator<(const interval& i, const interval& j)
{
    return (i <= j) && (i != j);
}

inline bool operator>=(const interval& i, const interval& j)
{
    return j <= i;
}

inline bool operator>(const interval& i, const interval& j)
{
    return j < i;
}
