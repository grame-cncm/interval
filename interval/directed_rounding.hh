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

// LSB is the legacy integer marker, but an out-of-range or infinite bound must
// never reach a C++ int cast. Such corridors are handled numerically instead.
inline bool hasInt32Bounds(const interval& x)
{
    return x.lsb() >= 0 && x.lo() >= -2147483648.0 && x.hi() <= 2147483647.0;
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

// Numeric image on the valid real domain, with outward binary64 bounds. NaN is
// still represented by the library's historical empty convention, not tracked.
// The LSB stamp is a fallback estimate, separate from the bound inclusion proof.
interval doubleUnaryBounds(UnaryOp op, const interval& x);
interval doubleAtan2Bounds(const interval& y, const interval& x);
interval doublePowBounds(const interval& x, const interval& y);
interval doubleFmodBounds(const interval& x, const interval& y);

}  // namespace itv::detail
