/* Copyright 2020-2026 Yann Orlarey, Agathe Herrou, Stéphane Letz
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
#include <cmath>

#include "directed_rounding.hh"
#include "interval_algebra.hh"
#include "interval_def.hh"

namespace itv {

//------------------------------------------------------------------------------------------
// The compound operations Frac, Cfrac, Mix and Tpt : each one is its definition, written once as
// an operation so that its repeated operand appears once. Evaluated as two independent
// values, x - floor(x) loses that the result lies in [0, 1] (the dependency problem).
//
// The bounds describe the program, which computes x once and then the operation : the
// subtraction of floor(x) or rint(x) is exact (Sterbenz : the rounded value is within a
// factor 2 of x, or 0), so the result lies in [0, 1[ (frac of x >= 0), [0, 1] (any x) or
// [-1/2, 1/2] (cfrac). One error is left : a compiler may contract the product computing x
// with the subtraction (fma(a, b, -floor(a * b))), reading x exactly on one side and
// rounded on the other ; the result may then leave its range by half an ulp of x. The
// bounds keep that margin.
//------------------------------------------------------------------------------------------

namespace {

/// the unit roundoff of the program (float for the precisions without one)
double unitRoundoff()
{
    return programPrecision() == 2 ? std::ldexp(1.0, -53) : std::ldexp(1.0, -24);
}

/// the largest program value below 1
double belowOne()
{
    if (programPrecision() == 1) return double(std::nextafter(1.0f, 0.0f));
    if (programPrecision() == 2) return std::nextafter(1.0, 0.0);
    return 1;
}

/// half an ulp of the largest magnitude of x, plus the loss of a flush to zero
double contractionMargin(const interval& x)
{
    const double mag = std::max(std::abs(x.lo()), std::abs(x.hi()));
    return unitRoundoff() * mag * (1 + std::ldexp(1.0, -20)) + std::ldexp(1.0, -125);
}

/// x - round(x) for a rounding round (floor or rint) of range [lo, hi] when x spans
/// several cells ; [x.lo - r, x.hi - r] when x lies in the cell of the integer r
interval remainder(const interval& input, double (*round)(double), double lo, double hi)
{
    const interval x = detail::floatingOperand(input);
    if (x.isEmpty()) return empty();
    if (input.lsb() >= 0) return interval(0, 0, -1);  // an integer is its own rounding
    if (!std::isfinite(x.lo()) || !std::isfinite(x.hi())) return interval(lo, hi, -1);
    const double rl = round(x.lo()), rh = round(x.hi());
    double       a = lo, b = hi;
    if (rl == rh) {
        a = std::max(lo, detail::directedBinary(detail::BinaryOp::Sub, x.lo(), rl, detail::Direction::Down));
        b = std::min(hi, detail::directedBinary(detail::BinaryOp::Sub, x.hi(), rl, detail::Direction::Up));
    }
    const double d = contractionMargin(x);
    return interval(a - d, b + d, std::min(x.lsb(), -1));
}

double floorOf(double v)
{
    return std::floor(v);
}

double rintOf(double v)
{
    return std::rint(v);
}

}  // namespace

// The definitions, operation by operation
interval interval_algebra::mixDefinition(const interval& c, const interval& a, const interval& b) const
{
    const interval w = Sub(FloatNum(1.0), c);
    const interval p = Mul(w, a);
    const interval q = Mul(c, b);
    return Add(p, q);
}

interval interval_algebra::tptDefinition(const interval& x) const
{
    const interval d = Add(x, FloatNum(1.0));
    return Div(x, d);
}

// Numeric kernel : x - floor(x) lies in [0, 1[ for x >= 0, [0, 1] otherwise
interval interval_algebra::numericFrac(const interval& x) const
{
    const bool nonnegative = !x.isEmpty() && x.lo() >= 0;
    return remainder(x, floorOf, 0, nonnegative ? belowOne() : 1.0);
}

// Numeric kernel : x - rint(x) lies in [-1/2, 1/2]
interval interval_algebra::numericCfrac(const interval& x) const
{
    return remainder(x, rintOf, -0.5, 0.5);
}

// Numeric kernel : (1 - c) a + c b with c in [0, 1] lies between a and b. For a fixed c
// the mix is linear in c, so its extremes are reached at the two ends of c : computed in
// double (a 2^-50 relative slack covers their own rounding), then widened by 3 rounding
// errors of the program at any c in between (1 - c, the two products, the sum). The
// weights are nonnegative : a mix of nonnegative operands stays nonnegative whatever the
// rounding (and symmetrically). Outside of its premise (c beyond [0, 1], an unbounded
// operand), the mix is its definition.
interval interval_algebra::numericMix(const interval& c, const interval& a, const interval& b) const
{
    if (c.isEmpty() || a.isEmpty() || b.isEmpty()) return empty();
    const interval fc = detail::floatingOperand(c), fa = detail::floatingOperand(a),
                   fb = detail::floatingOperand(b);
    if (fc.lo() < 0 || fc.hi() > 1 || !fa.isBounded() || !fb.isBounded()) {
        return mixDefinition(c, a, b);
    }
    double lo = HUGE_VAL, hi = -HUGE_VAL;
    for (double e : {fc.lo(), fc.hi()}) {
        const double w = 1 - e;
        lo             = std::min(lo, w * fa.lo() + e * fb.lo());
        hi             = std::max(hi, w * fa.hi() + e * fb.hi());
    }
    const double mag = std::max({std::abs(fa.lo()), std::abs(fa.hi()), std::abs(fb.lo()), std::abs(fb.hi())});
    const double eps = (3 * unitRoundoff() + std::ldexp(1.0, -50)) * mag + std::ldexp(1.0, -122);
    double       rlo = lo - eps, rhi = hi + eps;
    if (fa.lo() >= 0 && fb.lo() >= 0) rlo = std::max(rlo, 0.0);
    if (fa.hi() <= 0 && fb.hi() <= 0) rhi = std::min(rhi, 0.0);
    return interval(rlo, rhi, std::min({fc.lsb(), fa.lsb(), fb.lsb()}));
}

// Numeric kernel : x / (x + 1) is increasing on ]-1, +inf[, so for a finite x in that
// domain its image is [lo / (lo + 1), hi / (hi + 1)]. The ends are computed in double
// (a 2^-50 relative slack covers their own rounding), then widened by the 2 rounding
// errors of the program (the sum, the division) relative to each end ; the result is
// also intersected with the bounds of the definition, both being sound. Outside of
// that domain (x reaching -1, or unbounded), Tpt is its definition.
interval interval_algebra::numericTpt(const interval& x) const
{
    if (x.isEmpty()) return empty();
    const interval def = tptDefinition(x);
    const interval fx  = detail::floatingOperand(x);
    if (!(fx.lo() > -1) || !std::isfinite(fx.hi()) || def.isEmpty()) return def;
    const double rel = 2 * unitRoundoff() + std::ldexp(1.0, -50);
    const double lo  = fx.lo() / (fx.lo() + 1);
    const double hi  = fx.hi() / (fx.hi() + 1);
    const double rlo = lo - std::abs(lo) * rel - std::ldexp(1.0, -122);
    const double rhi = hi + std::abs(hi) * rel + std::ldexp(1.0, -122);
    return interval(std::max(rlo, def.lo()), std::min(rhi, def.hi()), def.lsb());
}

// Public API : an x that may be infinite gives a possible NaN (inf - inf)
interval interval_algebra::Frac(const interval& x) const
{
    const bool infinite = !x.isEmpty() && (!std::isfinite(x.lo()) || !std::isfinite(x.hi()));
    return numericFrac(x).withInvalid(x.mayBeInvalid() || infinite);
}

interval interval_algebra::Cfrac(const interval& x) const
{
    const bool infinite = !x.isEmpty() && (!std::isfinite(x.lo()) || !std::isfinite(x.hi()));
    return numericCfrac(x).withInvalid(x.mayBeInvalid() || infinite);
}

interval interval_algebra::Mix(const interval& c, const interval& a, const interval& b) const
{
    return numericMix(c, a, b).withInvalid(c.mayBeInvalid() || a.mayBeInvalid() || b.mayBeInvalid());
}

// Public API : inside its domain Tpt is valid when x is ; outside, the definition says
interval interval_algebra::Tpt(const interval& x) const
{
    if (x.isEmpty()) return empty();
    const interval fx = detail::floatingOperand(x);
    if (!(fx.lo() > -1) || !std::isfinite(fx.hi())) return tptDefinition(x);
    return numericTpt(x).withInvalid(x.mayBeInvalid());
}

}  // namespace itv
