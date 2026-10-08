![C/C++ CI](https://github.com/orlarey/interval/actions/workflows/ubuntu.yml/badge.svg)![C/C++ CI](https://github.com/orlarey/interval/actions/workflows/macos.yml/badge.svg)![C/C++ CI](https://github.com/orlarey/interval/actions/workflows/windows.yml/badge.svg)

# Intervals

The interval library of the Faust compiler. Since July 2026 it is *the* interval
domain of the compiler's type system: every signal's range is produced here, by the
fixpoint engine running over these algebras.

Licensed under Apache-2.0 (see LICENSE). Copyright 2020-2026 Yann Orlarey, Agathe
Herrou, Stéphane Letz.

## The two roles of the interval computation

- **Correctness** — the program runs right: delay lines sized from provable bounds,
  table indices certified within `[0, size)`, operation-domain and conversion
  validity, compile-time error reporting. For reliable static validation, an
  uncertified table index must cause a diagnostic rejection; explicit index
  clamping belongs in the author's program. These are the intended guarantees,
  subject to the limitations described below.
- **Sound quality** — the precision of the computations: in fixed point (FPGA
  targets), a signal's format takes its integer bits (msb) from its range and its
  fractional bits from the lsb analysis, so a tighter interval converts directly into
  signal-to-noise ratio. See BACKGROUND.md and PRECISION.md.

## Correctness criterion (normative)

**Every interval implementation must satisfy the following inclusion and
validity requirements.** They define the mathematical result to approximate,
independently of the algorithm used to compute it.

For a unary operation, let:

- `X` denote the set of numeric values represented by the input interval;
- `b_X` denote its `mayBeInvalid` flag;
- `F` denote the operation under the chosen numeric contract;
- `D` denote its admissible domain: inputs for which `F` is defined and produces
  no NaN;
- `S` denote the permitted representable bounds, for example binary32 or binary64
  values, completed with negative and positive infinity. NaN is not a bound.

The contract fixes the meaning of `F` and `D`, including numeric types,
rounding, wrapping and conversions. Thus `F` need not be an exact real-valued
mathematical function. Any further restriction, such as excluding infinite
results, must also be part of this contract.

Define the projections onto representable bounds by:

$$
\lfloor a\rfloor_S = \max\{s\in S\mid s\le a\},
\qquad
\lceil a\rceil_S = \min\{s\in S\mid s\ge a\}.
$$

The set of valid results and its ideal enclosure are:

$$
V = \{F(x)\mid x\in X\cap D\},
$$

$$
Y =
\begin{cases}
\varnothing & \text{if } V=\varnothing,\\
[\lfloor\inf V\rfloor_S,\;\lceil\sup V\rceil_S] & \text{otherwise}.
\end{cases}
$$

`Y` is the smallest interval with bounds in `S` containing every valid result.
The validity rule retains input invalidity and records inputs outside the domain:

$$
b_Y = b_X \lor (X\setminus D\ne\varnothing).
$$

If an implementation returns an interval $\widehat Y$ and a flag
$\widehat b_Y$, **its correctness criterion is**:

$$
Y\subseteq\widehat Y,
\qquad
b_Y\Rightarrow\widehat b_Y.
$$

An implementation may return a wider interval or conservatively report a
possible invalidity. It must never exclude a required value or report `false`
when the validity rule requires `true`. Failure to establish domain membership
must not be silently treated as proof of validity.

Empty numeric bounds with a true flag describe possible invalid execution with
no valid numeric result. They are distinct from an empty set of possibilities
with a false flag. Numerical enclosure alone cannot erase an invalidity.

For operations with several arguments, use the set of possible argument tuples
(the Cartesian product when only separate intervals are known), intersect it
with the joint domain `D`, and combine input flags with OR. Checking each
argument independently is not sufficient for restrictions such as `0 / 0` or
`0 * infinity`.

For Faust, certifying an index requires both a false invalidity flag and a
numeric enclosure entirely inside its admissible integer range. A dangerous or
uncertified index must be rejected with a diagnostic; silently clamping it does
not satisfy that validation policy.

**This section specifies the required contract, not a claim that the current
implementation fully satisfies it.** Both interval layers carry `mayBeInvalid`
and include it in their ordering. Domain checks cover NaN-producing numeric
operations and undefined int32 conversions, remainders and shift counts under
the declared contract. Target-libm guarantees and recurrence verification remain
obligations described below. The integer nature convention remains `lsb >= 0`;
the validity flag is independent of `lsb`.

## Two layers

1. **Ordinary intervals** (`interval_def.hh`, `interval_algebra.hh`):
   `<lo, hi, lsb, mayBeInvalid>` — numeric bounds, least-significant-bit precision
   and possible invalidity — with one
   implementation file per Faust primitive (`intervalXXX.cpp`). `interval_algebra`
   implements the full `FaustAlgebra<interval>` interface (the primitive set of the
   Faust signal language, vendored in `FaustAlgebra/`).
2. **Affine-in-time intervals** (`affint.hh`, `affine_ops.hh`): bounds of the form
   `[a0 + a1·t, b0 + b1·t]` over a lifetime `[0, T]` — the smallest useful fragment
   of the polyhedral domain. Linear operations work on the coefficients (a delay
   shifts the intercept by `-n·rate`, which is what makes accumulators stationary);
   nonlinear operations collapse to the hull over `[0, T]` and delegate to
   `interval_algebra` as the semantic oracle. This layer is what lets the compiler
   DATE the failure modes of a program (int32 wrap, float absorption) instead of
   merely bounding its values — the *horizon* analysis.

## Semantics of the COMPUTATIONAL interval

The interval describes what the *compiled program* can produce, not the mathematical
range:

- **Integer operations** model signed 32-bit arithmetic with two's-complement
  wrapping. Both the analyzer and generated C/C++ integer arithmetic must preserve
  that contract, for example with `-fwrapv` on GCC/Clang.
- **Integer modulo** (both operands with `lsb >= 0`) follows C semantics: the result
  has the sign of `x`, its magnitude is at most `min(|x|max, m-1)` with `m` the
  largest divisor magnitude, and `x % y = x` whenever `|x|` stays below the smallest
  nonzero divisor magnitude. The float path models `fmod`, whose bound closes at
  `|y|`.

- **Floating-point operations** use `itv::programPrecision()`: 2 double by default,
  1 single, 3 quad, 4 fixed point. Single and double use directed reference
  calculations. Single narrows computed bounds outward to binary32 after each
  operation; explicitly known literals/conversions use nearest rounding.
  Quad and fixed point retain their existing rules.
- **Target-libm compensation** uses `itv::libmCompensation()`, enabled by default.
  The historical margin is two ULPs of the program precision, clipped to known
  function images. It is an assumption, not a universal bound for arbitrary
  libms. In single and double mode, the native reference kernel is used even when
  compensation is disabled; computed singletons no longer bypass the margin.
  Quad/fixed retain their previous singleton exemption and host-libm evaluation.

Compiler reassociation, target approximations and recurrence invariants require
an explicit execution contract. A local bound calculation alone cannot certify
all of those properties. The reference arithmetic is kept separate from these
remaining obligations.

## Conservative binary64 bound calculations

With `itv::programPrecision() == 2`, the native kernel computes outward binary64
bounds **without MPFR/GMP**. This covers the analyzer's own rounding error: if the
desired lower bound belongs to `[L-, L+]` and the upper bound to `[U-, U+]`, the
retained interval is `[L-, U+]`. It is not enough to widen only the final result
after intermediate calculations have already rounded inward.

The elementary kernel uses TwoSum for nearest-rounded addition/subtraction and
exact comparisons of binary significands for multiplication, division, square
root, and power-of-two scaling. A product needs only two 64-bit integer limbs,
not a multiprecision library or hardware FMA. These comparisons remain valid
when a floating residual would underflow to zero. The rounded result is moved
by one `nextafter` only when the requested direction requires it. Non-nearest
addition/subtraction instead use a conservative neighbour fallback.

The `*Bounds` transcendental reference images use interval polynomial/series
calculations with analytic remainder bounds, independently of the host libm:

- `exp`: reduce by enclosing `ln(2)`, sum through degree 32 on `|r| <= 1`,
  bound the tail by three times the absolute degree-33 term, then scale outward.
- `log`: exact `frexp` reduction to `m` in `[1, 2)`, then 32 terms of the
  atanh series in `z = (m-1)/(m+1)`. Bound the remaining positive tail by
  `2*z^65 / (65*(1-z^2))` using outward arithmetic.
- `atan`: reciprocal/quadrant identities reduce to `|z| <= 1/2`; the 32-term
  alternating series has a tail bounded by its first omitted term.
- `sin`/`cos`: reduce with an enclosing pi and a bounded integer quadrant,
  evaluate 16 Taylor terms on `|r| <= 1`, and bound the alternating tail by
  the first omitted term. `tan` divides these images; inverse/hyperbolic
  functions use identities built from the same reference bounds.
- Integer powers use outward exponentiation by squaring; other powers use
  enclosing logarithms and exponentials. `log10` divides by an enclosing `ln(10)`.

These are conservative reference algorithms, **not correctly-rounded libm
replacements**. Bounds can span several ULPs and cancellation can widen them
further. Trigonometric arguments with magnitude greater than `2^20`, or an
uncertain reduction, use the full function image (`[-1, 1]` or an unrestricted
tangent corridor). Extrema and tangent poles are tested using enclosing pi
quotients; uncertain large phases deliberately include critical points.

`fmod`/`remainder` scalar values use exact binary shift/subtract division,
including huge quotients and ties-to-even. Interval remainder bounds round
half-divisors outward. Affine coefficients, line evaluation, chords, shifts,
and inclusion comparisons use the same elementary kernel. Integer operations
retain their int32 wrapping contract.

An injected literal is already a known machine value: `FloatNum(0.1)` remains a
singleton. `FloatNum` preserves the floating nature of integral-looking literals,
so `FloatNum(65536) * FloatNum(65536)` does not take the int32 path. Exact
operations such as `0.5 * 0.25` remain singletons. An inexact operation such as
`1 + 2^-54` is enclosed by `[1, nextafter(1, +inf)]`; a rounded equality of
endpoints must not create a false constant. Exact transcendental identities
such as `sin(0)`, `exp(0)` and `log(1)` also remain points in `*Bounds`.

This is the first step towards a reliable computational analysis. It does not
certify arbitrary target libms or reassociation, track signed zeros separately,
prove recurrence invariants, or certify the LSB precision estimates.
Those obligations need separate verification before bounds justify safety
decisions. The quad/fixed paths retain their previous rules and need
their own audit. The new transcendental reference is intended for inclusion;
its performance and tightness in the complete Faust compiler need measurement.

Single precision now uses the same native reference kernel: enclose each
calculation in binary64, then take the greatest binary32 value below its lower
bound and the least binary32 value above its upper bound. This also applies to
computed singleton intervals, subnormal results, overflow, and integral-looking
float values above `2^24`. Outward composition includes both separate rounding
and FMA contraction, under a contract that forbids reassociation and preserves
gradual underflow. Explicit `FloatNum` literals are rounded once and remain
points; compiler-folded values must be injected through this literal API.
`FloatCast` describes a known monotone conversion. Mixed arithmetic converts its
integer operands before floating evaluation. Float zero and integral-valued
floating functions retain a negative LSB, so later operations cannot wrap as int32.

`*Bounds` provides a conservative mathematical reference, independently of the
target libm. Public transcendental operations still add the existing two-ULP
margin when `libmCompensation()` is enabled. Inclusion of **executed libm values**
therefore requires a target whose result is at most two representable steps
from the correctly rounded value, with the exact special values and function
ranges used by the clamp. This is the precise contract of the historical
`nextafter` margin; a published absolute error in ULPs must be translated into
that contract, especially at binade boundaries. Arbitrary libms remain uncertified.
Disabling the margin requires a correctly-rounded target. This distinction
applies equally to float and double. `mayBeInvalid` separately records possible
NaN or undefined execution; NaN endpoints encode only numeric emptiness.
Reciprocal and negative powers of merged signed zero use wider float/double
hulls to include both infinite signs.

Float affine arithmetic collapses floating add/subtract/multiply/divide and
casts to a constant interval hull over the horizon. Rounded staircases cannot
in general be bounded by chords through their rounded endpoints. This sacrifices
rate precision; recurrence invariants still require separate verification.

[FLOAT_COUNTEREXAMPLES.md](FLOAT_COUNTEREXAMPLES.md) records the six original
failures and their corrected bounds. The negative-delay example now yields
`[-1, 8195]`. All seven witnesses (including huge sine) are ordinary inclusion
tests in `FloatConservativenessTests`; CTest has no expected-gap mode.
`FloatBoundsTests` broadens coverage, and optional MPFR targets verify both the
mathematical enclosure and binary32 rounding, including subnormals and overflow.
Passing finite tests supports these rules; it is not a proof for every Faust
program, unrestricted compiler transformations, recurrence or target libm.

## Integer and floating nature

The sign of `lsb` is the nature marker in both simple and double precision:
`lsb >= 0` means integer, and `lsb < 0` means floating. Integral-valued floating
results, including zero, `round` and `rint`, must retain a negative marker.
Inject typed literals through `IntNum` and `FloatNum`; the scalar
`interval(double)` constructor infers its marker from the value and cannot
identify the source type.

An integer corridor leaving int32, whether through infinite widening or finite
horizon overflow, is normalized to `[INT_MIN, INT_MAX]` before integer arithmetic,
implicit floating conversion or comparison. This normalization preserves the
integer marker; leaving int32 never selects floating semantics. Initial zeros
in affine delays inherit the signal's marker. Soundfile metadata remains integer.
Bitwise operations and shifts produce integers; arithmetic right shift rounds
negative quotients downward, rather than keeping real-valued scaled bounds.
Shift counts outside `[0,31]` receive a full hull, which does not define an
invalid runtime shift.

`FaustIntegrationTests` checks these rules in both precisions, including a
post-widening LCG invariant and independent runtime witnesses for the second
noise channel of `rain`. The permanent reduced DSP and a deterministic output
driver live in `tests/faust/`; they can also be tested with a real Faust compiler:

```sh
python3 tests/faust/run.py --faust /path/to/faust/build/bin/faust \
  --reference /path/to/trusted/faust --faust-root /path/to/faust
```

The runner checks both the reduced program and the full example, in `cpp` and
`ocpp`, single and double, with 4096 frames per case. It rejects non-finite values,
a silent output channel, and runtime differences from the optional reference.
Use `--cxx` and `--cxxflags` to select the C++ compiler or sanitizers.

The separate [rejection corpus](tests/faust/REJECTIONS.md) requires Faust to
reject an unbounded `rwtable` write index even when its read index is safe:

```sh
python3 tests/faust/run-rejections.py --faust /path/to/faust/build/bin/faust
```

This test currently exposes a missing compiler rejection; it does not pass by
silently adding a runtime clamp and is separate from the standalone library tests.

## Validity API and compiler integration

`interval::mayBeInvalid()` and `AffItv::mayBeInvalid` are the same attribute.
`withInvalid(true)` adds an alert; `withInvalid(false)` never clears one.
The ordinary four-argument constructor explicitly accepts an incoming flag.
For example:

```cpp
itv::interval_algebra algebra;
const auto y = algebra.Sqrt(itv::interval(-1, 4));
// y contains [0, 2] and y.mayBeInvalid() is true.
const auto z = algebra.IntCast(y);
// z contains [0, 2], has an integer LSB, and retains the alert.
```

`empty()` / `aempty()` construct numeric bottom without an alert. An injected
NaN constructs an invalid-only value. `isEmpty()` tests the numeric part only;
ordinary `isValid()` additionally requires a false flag. An alerted point is not
a foldable constant. Equality, inclusion, joins and affine widening include the
flag, so a flag-only change cannot masquerade as fixpoint convergence. Numeric
intersection preserves existing alerts, rather than proving away earlier invalid
execution. `libmBounds`, casts and both affine bridges retain the attribute.

`IntCast` considers truncation toward zero: in double mode, `2147483647.75`
converts validly to `INT_MAX`; `2147483648` does not. An invalid-only input has
no valid numeric integer result; a partly valid input keeps its truncated valid
image and raises the flag. Returning a hardware-dependent sentinel or full int32
does not make an undefined C++ cast valid.

The wrapping integer-power path requires nonnegative exponents; a possible
negative exponent is uncertified rather than silently treated as zero. Floating
negative powers use the floating path and its own domain.

The default numeric contract allows IEEE infinities when they are not NaN:
`log(0)` and `1/0` need not raise this flag, whereas `0/0`, `0*infinity`, and
`infinity-infinity` do. An application requiring finite results must also check
finiteness. Integer wrapping is valid under the declared wrapping contract.

Unknown foreign floating storage and uncontracted foreign functions remain
uncertified. Table reads and indexed soundfile accesses cannot establish their
full domain in this value-only interface, which lacks resource extents; they
conservatively raise the flag. The compiler must provide domain-aware rules using
the actual table/channel/part dimensions. Writes check their supplied table size.
Inputs, UI zones and soundfile samples still rely on their declared host domains;
this flag does not enforce those domains at runtime. Widget declarations assume
finite ordered bounds, an initial value within them and a nonnegative finite step.

**Faust must preserve the flag in type construction, copying, serialization,
joins and recurrence verification, and consult it before validating indices or
folding constants.** Copying only `lo`, `hi` and `lsb` loses information. Reconstruct
numeric bottom with `empty(lsb).withInvalid(storedFlag)`, not a NaN-valued literal.
Invalid constant conversions must also be checked before folding erases their
source. Adding the field to this library does not update those compiler consumers
or its diagnostic rejection policy automatically.

`InvalidityTests` checks these rules in float and double, including ordinary and
affine transfers, partial domains, invalid-only values, flag-only convergence,
joint infinite arguments, cast thresholds and independent NaN witnesses. Optional
`InvalidityOracleTests` adds MPFR domain and truncation checks. See
[VALIDITY.md](VALIDITY.md) for the French implementation and migration notes.

## Building and optional MPFR oracle

The library and its public/affine headers support C++17, with **no MPFR/GMP
production dependency**. The complete test build uses C++20; the integration
regressions deliberately compile as C++17 to exercise that consumer contract:

```sh
cmake -S . -B build -DNOTIDY=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

Enable the independent 256-bit MPFR oracle explicitly:

```sh
cmake -S . -B build-oracle -DNOTIDY=ON -DINTERVAL_ENABLE_MPFR_TESTS=ON
cmake --build build-oracle
ctest --test-dir build-oracle --output-on-failure
```

Only the optional oracle executables link MPFR/GMP; the library and the
dependency-free test executables remain independent of them. For this optional test, install
`libmpfr-dev libgmp-dev` on Debian/Ubuntu, `mpfr` through Homebrew or MacPorts on
macOS, or `mpfr:x64-windows` through vcpkg on Windows. Use `CMAKE_PREFIX_PATH`
for a nonstandard installation or the vcpkg CMake toolchain on Windows. An
explicit oracle request fails if its dependencies are unavailable.

For WebAssembly, with Emscripten and Node available:

```sh
emcmake cmake -S . -B build-wasm -DNOTIDY=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-wasm
ctest --test-dir build-wasm --output-on-failure
```

The CMake target disables fast-math and implicit contraction for the analyzer
and its header-based helpers. Direct consumers must preserve this compilation
contract as well. The kernel does not change the CPU rounding mode and has no
mutable numerical state. IEEE gradual underflow is required (FTZ/DAZ disabled
in the analyzer); target FTZ/DAZ remains a separate model. Only the optional
oracle needs a thread-safe MPFR build or caller serialization if used concurrently.

## Special intervals

Division computes the four endpoint quotients directly, avoiding the extra rounding
of multiplication by a reciprocal. If the denominator contains zero, or an endpoint
quotient is indeterminate (`inf/inf`), it conservatively returns `[-inf, +inf]`.
The affine layer preserves rates for division by a constant corridor; a denominator
with a rate is collapsed over the full horizon before division, since reciprocal
curves and interior zero crossings cannot be bounded by endpoint interpolation.

- `interval()` is the historical default interval. It contains every finite value
  representable by a `double` and has an LSB of `-24`.
- `fullFinite()` constructs the same interval explicitly, optionally with a different LSB.
- `empty()` constructs numeric bottom with a false validity flag, using `NaN`
  endpoints as storage sentinels. `empty().withInvalid()` is invalid-only.
- `interval(NAN)` injects an invalid-only value; it is distinct from `empty()`.
- `interval(-HUGE_VAL, HUGE_VAL)` also contains infinities and is therefore distinct from
  `fullFinite()`.

## Organization of the code

All the code lives in the namespace `itv`:

- `interval_def.hh` — the interval data structure and its basic accessors
- `interval_algebra.hh/cpp` — all operations on intervals, as defined by the Faust
  primitives; one `intervalXXX.cpp` per operation
- `intervalValidity.cpp`, `validity.hh` — public validity-aware transfers and
  internal operation-domain predicates, separate from the numeric kernels
- `bitwiseOperations.hh/cpp` — helpers for the bitwise operations
- `affint.hh` — the affine-in-time interval (`AffItv`) and its numeric core
- `affine_ops.hh` — `AffineOps<Base>`, the FaustAlgebra operations over `AffItv`;
  `affine_algebra` is its standalone instantiation
- `check.hh/cpp` — test helpers (exact checks and sampled numerical analyses)
- `FaustAlgebra/` — vendored copy of the FaustAlgebra interface
- `tests/` — the deterministic regression suite (see TESTING.md)
