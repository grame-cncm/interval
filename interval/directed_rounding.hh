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
#pragma once

#include <utility>

#include "interval_def.hh"

namespace itv::detail {

// Internal numerical kernel. Inputs are exact binary64 values; directed results
// enclose the mathematical operation, including binary64 underflow and overflow.
// Exact dyadic comparisons and bounded reference series do not change host fenv.
enum class Direction { Down, Up };
enum class BinaryOp { Add, Sub, Mul, Div, Pow, Atan2, Fmod, Remainder };
enum class UnaryOp {
    Sqrt, Acos, Acosh, Asin, Asinh, Atan, Atanh, Cos, Cosh,
    Exp, Log, Log10, Sin, Sinh, Tan, Tanh
};

// Float reuses the binary64 proof, then interval construction narrows outward.
// Quad/fixed keep their old rules; the analyzer never narrows a reference step.
inline bool usesNativeBounds() { return programPrecision() == 1 || programPrecision() == 2; }

// Mixed floating/int expressions convert integer operands to the target float
// before evaluation. Enclosing that conversion avoids computing on a non-float
// integer endpoint (e.g. 16777217). Already-floating domains are unchanged.
inline interval floatingOperand(const interval& x)
{
    return programPrecision() == 1 && x.lsb() >= 0
        ? interval(programBound(x.lo()), programBound(x.hi()), -24) : x;
}

// A numerical floating result must not reenter int32 wrapping merely because
// its estimated grid has a nonnegative LSB. Explicit IntCast still marks integers.
inline int floatingLSB(int lsb)
{
    return programPrecision() == 1 ? std::min(lsb, -1) : lsb;
}

// LSB is the legacy integer marker, but an out-of-range or infinite bound must
// never reach a C++ int cast. Such corridors are handled numerically instead.
inline bool hasInt32Bounds(const interval& x)
{
    return x.lsb() >= 0 && x.lo() >= -2147483648.0 && x.hi() <= 2147483647.0;
}

// Under the wrapping semantics, an integer interval (nonnegative LSB) whose bounds
// leave int32 describes any int32 : an infinite bound from the widening of an integer
// recursion, or a finite one from the affine domain evaluated at its horizon (a
// counter reaches 2^31). Normalizing it keeps the integer path : the floating path
// ignores the wrap and certifies values the program never computes (12345 + [0, +inf]
// is not [12345, +inf] in int32). A float literal enters through FloatNum, whose
// negative LSB keeps it on the floating path. Applied when both operands are integers.
inline interval int32Hull(const interval& x)
{
    if (x.isEmpty() || x.lsb() < 0 || hasInt32Bounds(x)) return x;
    return {-2147483648.0, 2147483647.0, x.lsb()};
}

inline std::pair<interval, interval> int32Operands(const interval& x, const interval& y)
{
    if (x.lsb() >= 0 && y.lsb() >= 0) return {int32Hull(x), int32Hull(y)};
    return {x, y};
}

// Metadata must not overflow signed int when an underflowed power already has
// an extreme grid exponent. INT_MIN is reserved by interval's constructor.
inline int sumLSB(int a, int b)
{
    return static_cast<int>(std::clamp<int64_t>(int64_t(a) + int64_t(b), INT_MIN + 1LL, INT_MAX));
}

double directedBinary(BinaryOp op, double x, double y, Direction direction);
double directedUnary(UnaryOp op, double x, Direction direction);
double directedPi(Direction direction);

// Internal exact power-of-two scaling, rounded outward including subnormals.
// The exponent is a bounded internal range-reduction exponent, not arbitrary int.
double directedScale(double x, int exponent, Direction direction);

// Numeric image on the valid real domain: outward binary64 reference steps,
// then interval construction narrows outward in float mode. NaN is
// still represented by the library's historical empty convention, not tracked.
// The LSB stamp is a fallback estimate, separate from the bound inclusion proof.
interval doubleUnaryBounds(UnaryOp op, const interval& x);
interval doubleAtan2Bounds(const interval& y, const interval& x);
interval doublePowBounds(const interval& x, const interval& y);
interval doubleFmodBounds(const interval& x, const interval& y);

}  // namespace itv::detail
