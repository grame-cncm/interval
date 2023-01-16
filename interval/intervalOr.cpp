/* Copyright 2023 Yann ORLAREY
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include <algorithm>
#include <functional>
#include <random>

#include "bitwiseOperations.hh"
#include "check.hh"
#include "interval_algebra.hh"
#include "interval_def.hh"

namespace itv {
//------------------------------------------------------------------------------------------
// Interval Or
// interval Or(const interval& x, const interval& y) const;
// void testOr() const;
static double myOr(double x, double y)
{
    int a = saturatedIntCast(x);
    int b = saturatedIntCast(y);
    int c = a | b;
    return double(c);
}

// BRUTE FORCE
interval interval_algebra::Or(const interval& x, const interval& y) const
{
    if (x.isEmpty() || y.isEmpty()) return {};
    int x0 = saturatedIntCast(x.lo());
    int x1 = saturatedIntCast(x.hi());
    int y0 = saturatedIntCast(y.lo());
    int y1 = saturatedIntCast(y.hi());

    SInterval z = bitwiseSignedOr({x0, x1}, {y0, y1});
    return {double(z.lo), double(z.hi)};
}

void interval_algebra::testOr() const
{
    analyzeBinaryMethod(10, 2000, "Or", interval(-1000, -800), interval(127), myOr, &interval_algebra::Or);
    analyzeBinaryMethod(10, 2000, "Or", interval(-1000, -800), interval(123), myOr, &interval_algebra::Or);
    analyzeBinaryMethod(10, 2000, "Or", interval(-128, 128), interval(127), myOr, &interval_algebra::Or);
    analyzeBinaryMethod(10, 2000, "Or", interval(0, 1000), interval(63, 127), myOr, &interval_algebra::Or);
    analyzeBinaryMethod(10, 2000, "Or", interval(-1000, 1000), interval(63, 127), myOr, &interval_algebra::Or);
}
}  // namespace itv
