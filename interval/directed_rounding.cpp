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
#include "directed_rounding.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <mpfr.h>

namespace itv::detail {
namespace {

static_assert(std::numeric_limits<double>::is_iec559 &&
              std::numeric_limits<double>::digits == 53 &&
              std::numeric_limits<double>::max_exponent == 1024,
              "Directed double bounds require IEEE binary64 storage.");

// An embedding application can restrict MPFR's thread-local exponent range.
// Never load a double inexactly under that restriction; use an unrestricted
// numeric bound instead of changing the caller's MPFR configuration.
bool usableExponentRange()
{
    return mpfr_get_emin() <= -1073 && mpfr_get_emax() >= 1024;
}

// Local MPFR values own their storage. Binary64 inputs load exactly at 53 bits;
// a wider exponent range prevents premature underflow before directed get_d.
class Number {
   public:
    mpfr_t value;
    explicit Number(mpfr_prec_t precision = 53) { mpfr_init2(value, precision); }
    ~Number() { mpfr_clear(value); }
    Number(const Number&) = delete;
    Number& operator=(const Number&) = delete;
};

mpfr_rnd_t mode(Direction direction)
{
    return direction == Direction::Down ? MPFR_RNDD : MPFR_RNDU;
}

// Test whether a critical point pi*(offset + period*k) may occur in [lo, hi].
// Interval division by enclosing values of pi avoids false negatives from host
// argument reduction; at huge phases uncertainty deliberately includes extrema.
bool containsPiLattice(double lo, double hi, double offset, unsigned period)
{
    if (!usableExponentRange()) return true;
    if (!std::isfinite(lo) || !std::isfinite(hi)) return true;
    Number piLo(128), piHi(128), lower(128), upper(128), input(128);
    mpfr_const_pi(piLo.value, MPFR_RNDD);
    mpfr_const_pi(piHi.value, MPFR_RNDU);
    mpfr_set_d(input.value, lo, MPFR_RNDN);
    mpfr_div(lower.value, input.value, lo < 0 ? piLo.value : piHi.value, MPFR_RNDD);
    mpfr_sub_d(lower.value, lower.value, offset, MPFR_RNDD);
    mpfr_div_ui(lower.value, lower.value, period, MPFR_RNDD);
    mpfr_set_d(input.value, hi, MPFR_RNDN);
    mpfr_div(upper.value, input.value, hi < 0 ? piHi.value : piLo.value, MPFR_RNDU);
    mpfr_sub_d(upper.value, upper.value, offset, MPFR_RNDU);
    mpfr_div_ui(upper.value, upper.value, period, MPFR_RNDU);
    mpfr_ceil(lower.value, lower.value);
    mpfr_floor(upper.value, upper.value);
    return mpfr_cmp(lower.value, upper.value) <= 0;
}

interval positivePower(double lo, double hi, double elo, double ehi, int lsb)
{
    double lower = HUGE_VAL, upper = -HUGE_VAL;
    for (double base : {lo, hi}) {
        for (double exponent : {elo, ehi}) {
            const double a = directedBinary(BinaryOp::Pow, base, exponent, Direction::Down);
            const double b = directedBinary(BinaryOp::Pow, base, exponent, Direction::Up);
            // Indeterminate extended-real corners must not discard finite values
            // attained inside the rectangle (e.g. 1^inf or inf^0).
            if (std::isnan(a) || std::isnan(b)) return {0, HUGE_VAL, lsb};
            lower = std::min(lower, a);
            upper = std::max(upper, b);
        }
    }
    return {lower, upper, lsb};
}

}  // namespace

// Internal kernel: enclose the exact binary operation in the requested direction.
// Both MPFR evaluation and conversion are directed, so subnormal rounding cannot
// undo the inclusion established at MPFR's wider exponent range.
double directedBinary(BinaryOp op, double x, double y, Direction direction)
{
    if (!usableExponentRange()) return direction == Direction::Down ? -HUGE_VAL : HUGE_VAL;
    Number a, b, result;
    mpfr_set_d(a.value, x, MPFR_RNDN);
    mpfr_set_d(b.value, y, MPFR_RNDN);
    const auto rounding = mode(direction);
    switch (op) {
        case BinaryOp::Add: mpfr_add(result.value, a.value, b.value, rounding); break;
        case BinaryOp::Sub: mpfr_sub(result.value, a.value, b.value, rounding); break;
        case BinaryOp::Mul: mpfr_mul(result.value, a.value, b.value, rounding); break;
        case BinaryOp::Div: mpfr_div(result.value, a.value, b.value, rounding); break;
        case BinaryOp::Pow: mpfr_pow(result.value, a.value, b.value, rounding); break;
        case BinaryOp::Atan2: mpfr_atan2(result.value, a.value, b.value, rounding); break;
        case BinaryOp::Fmod: mpfr_fmod(result.value, a.value, b.value, rounding); break;
        case BinaryOp::Remainder: mpfr_remainder(result.value, a.value, b.value, rounding); break;
    }
    return mpfr_get_d(result.value, rounding);
}

// Internal kernel: enclose the exact unary operation; invalid arguments give NaN.
double directedUnary(UnaryOp op, double x, Direction direction)
{
    if (!usableExponentRange()) return direction == Direction::Down ? -HUGE_VAL : HUGE_VAL;
    Number input, result;
    mpfr_set_d(input.value, x, MPFR_RNDN);
    const auto rounding = mode(direction);
    switch (op) {
        case UnaryOp::Sqrt: mpfr_sqrt(result.value, input.value, rounding); break;
        case UnaryOp::Acos: mpfr_acos(result.value, input.value, rounding); break;
        case UnaryOp::Acosh: mpfr_acosh(result.value, input.value, rounding); break;
        case UnaryOp::Asin: mpfr_asin(result.value, input.value, rounding); break;
        case UnaryOp::Asinh: mpfr_asinh(result.value, input.value, rounding); break;
        case UnaryOp::Atan: mpfr_atan(result.value, input.value, rounding); break;
        case UnaryOp::Atanh: mpfr_atanh(result.value, input.value, rounding); break;
        case UnaryOp::Cos: mpfr_cos(result.value, input.value, rounding); break;
        case UnaryOp::Cosh: mpfr_cosh(result.value, input.value, rounding); break;
        case UnaryOp::Exp: mpfr_exp(result.value, input.value, rounding); break;
        case UnaryOp::Log: mpfr_log(result.value, input.value, rounding); break;
        case UnaryOp::Log10: mpfr_log10(result.value, input.value, rounding); break;
        case UnaryOp::Sin: mpfr_sin(result.value, input.value, rounding); break;
        case UnaryOp::Sinh: mpfr_sinh(result.value, input.value, rounding); break;
        case UnaryOp::Tan: mpfr_tan(result.value, input.value, rounding); break;
        case UnaryOp::Tanh: mpfr_tanh(result.value, input.value, rounding); break;
    }
    return mpfr_get_d(result.value, rounding);
}

// Internal kernel: bound pi itself rather than treating the nearest M_PI as exact.
double directedPi(Direction direction)
{
    if (!usableExponentRange()) return direction == Direction::Down ? -HUGE_VAL : HUGE_VAL;
    Number result;
    mpfr_const_pi(result.value, mode(direction));
    return mpfr_get_d(result.value, mode(direction));
}

// Numeric image on the valid real domain; endpoint evaluation is supplemented by
// all interior extrema/poles. LSB remains an estimate, not part of this proof.
interval doubleUnaryBounds(UnaryOp op, const interval& x)
{
    if (x.isEmpty()) return empty();
    double lo = x.lo(), hi = x.hi();
    const int lsb = std::min(x.lsb(), -24);
    switch (op) {
        case UnaryOp::Sqrt: case UnaryOp::Log: case UnaryOp::Log10:
            if (hi < 0) return empty();
            lo = std::max(lo, 0.0);
            break;
        case UnaryOp::Acosh:
            if (hi < 1) return empty();
            lo = std::max(lo, 1.0);
            break;
        case UnaryOp::Acos: case UnaryOp::Asin: case UnaryOp::Atanh:
            if (hi < -1 || lo > 1) return empty();
            lo = std::max(lo, -1.0);
            hi = std::min(hi, 1.0);
            break;
        default: break;
    }
    if (op == UnaryOp::Sin || op == UnaryOp::Cos || op == UnaryOp::Tan) {
        // Infinities are invalid trigonometric arguments. An unbounded interval
        // still contains finite arguments, whose image must remain represented.
        if (lo == hi && std::isinf(lo)) return empty();
        if (op == UnaryOp::Tan && containsPiLattice(lo, hi, 0.5, 1))
            return {-HUGE_VAL, HUGE_VAL, lsb};
        if (!std::isfinite(lo) || !std::isfinite(hi)) return {-1, 1, lsb};
        double lower = std::min(directedUnary(op, lo, Direction::Down),
                                directedUnary(op, hi, Direction::Down));
        double upper = std::max(directedUnary(op, lo, Direction::Up),
                                directedUnary(op, hi, Direction::Up));
        if (op == UnaryOp::Sin) {
            if (containsPiLattice(lo, hi, 0.5, 2)) upper = 1;
            if (containsPiLattice(lo, hi, -0.5, 2)) lower = -1;
        } else if (op == UnaryOp::Cos) {
            if (containsPiLattice(lo, hi, 0, 2)) upper = 1;
            if (containsPiLattice(lo, hi, 1, 2)) lower = -1;
        }
        return {lower, upper, lsb};
    }
    if (op == UnaryOp::Cosh) {
        return {x.hasZero() ? 1 : std::min(directedUnary(op, lo, Direction::Down),
                                          directedUnary(op, hi, Direction::Down)),
                std::max(directedUnary(op, lo, Direction::Up),
                         directedUnary(op, hi, Direction::Up)), lsb};
    }
    if (op == UnaryOp::Acos) std::swap(lo, hi); // decreasing on [-1, 1]
    return {directedUnary(op, lo, Direction::Down), directedUnary(op, hi, Direction::Up), lsb};
}

// Numeric atan2 image. Rectangular extrema are corners except at the negative
// x-axis cut or the origin; include both signs of pi whenever that cut is possible.
interval doubleAtan2Bounds(const interval& y, const interval& x)
{
    if (x.isEmpty() || y.isEmpty()) return empty();
    const int lsb = std::min({x.lsb(), y.lsb(), -24});
    if (x.lo() <= 0 && y.hasZero()) {
        const double pi = directedPi(Direction::Up);
        return {-pi, pi, lsb};
    }
    double lower = HUGE_VAL, upper = -HUGE_VAL;
    for (double a : {y.lo(), y.hi()}) {
        for (double b : {x.lo(), x.hi()}) {
            lower = std::min(lower, directedBinary(BinaryOp::Atan2, a, b, Direction::Down));
            upper = std::max(upper, directedBinary(BinaryOp::Atan2, a, b, Direction::Up));
        }
    }
    return {lower, upper, lsb};
}

// Numeric pow image. Positive bases use monotonicity in each coordinate (the
// direction changes at base 1/exponent 0, but their value 1 is covered by corners).
// Negative bases contribute only integer exponents, separated by parity.
interval doublePowBounds(const interval& x, const interval& y)
{
    if (x.isEmpty() || y.isEmpty()) return empty();
    const int lsb = std::min({x.lsb(), y.lsb(), -24});
    interval result = empty();
    if (x.hi() >= 0) result = positivePower(std::max(0.0, x.lo()), x.hi(), y.lo(), y.hi(), lsb);
    // MPFR/IEEE pow also defines some noninteger powers of -inf, although the
    // corresponding finite negative bases have no real image. Keep those numeric
    // exceptional results instead of dropping the infinite endpoint entirely.
    if (x.lo() == -HUGE_VAL) {
        if (y.hi() > 0) result = reunion(result, interval(-HUGE_VAL, HUGE_VAL, lsb));
        if (y.lo() < 0) result = reunion(result, interval(0));
        if (y.hasZero()) result = reunion(result, interval(1));
    }
    if (x.lo() < 0) {
        double first = std::ceil(y.lo()), last = std::floor(y.hi());
        if (first <= last) {
            // Beyond 2^53, binary64 cannot distinguish neighbouring integer
            // parities. An unrestricted hull is preferable to dropping a sign.
            if (!std::isfinite(first) || !std::isfinite(last) ||
                std::abs(first) >= 0x1p53 || std::abs(last) >= 0x1p53)
                return {-HUGE_VAL, HUGE_VAL, lsb};
            for (int parity : {0, 1}) {
                const double a = first + (std::fmod(std::abs(first), 2.0) == parity ? 0 : 1);
                const double b = last - (std::fmod(std::abs(last), 2.0) == parity ? 0 : 1);
                if (a > b) continue;
                const interval magnitude = positivePower(std::max(0.0, -x.hi()), -x.lo(), a, b, lsb);
                result = reunion(result, parity == 0 ? magnitude
                                                    : interval(-magnitude.hi(), -magnitude.lo(), lsb));
            }
        }
    }
    return result;
}

// Numeric fmod image. The sign follows x and |fmod| <= min(|x|, |y|); this
// deliberately avoids a rounded quotient and an unsafe conversion to int.
interval doubleFmodBounds(const interval& x, const interval& y)
{
    if (x.isEmpty() || y.isEmpty() || y.isZero()) return empty();
    const int lsb = std::min({x.lsb(), y.lsb(), -24});
    if (x.isconst() && y.isconst()) {
        return {directedBinary(BinaryOp::Fmod, x.lo(), y.lo(), Direction::Down),
                directedBinary(BinaryOp::Fmod, x.hi(), y.hi(), Direction::Up), lsb};
    }
    const double magnitude = std::min(std::max(std::abs(x.lo()), std::abs(x.hi())),
                                      std::max(std::abs(y.lo()), std::abs(y.hi())));
    return {x.lo() >= 0 ? 0 : -magnitude, x.hi() <= 0 ? 0 : magnitude, lsb};
}

}  // namespace itv::detail
