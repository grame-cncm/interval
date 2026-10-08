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
#include "interval_algebra.hh"
#include "interval_def.hh"

namespace itv {
//------------------------------------------------------------------------------------------
// Interval IntNum

// Public API: inject the known binary64 literal as a point. In double mode retain
// floating nature even for integral-looking literals; other modes keep legacy
// injection. A known literal needs no widening for analyzer rounding error.
// NaN keeps the library's historical empty representation.
interval interval_algebra::FloatNum(double x) const
{
    const interval value = singleton(x);
    if (programPrecision() != 2) return value;
    return {value.lo(), value.hi(), std::min(value.lsb(), -1)};
}
}  // namespace itv
