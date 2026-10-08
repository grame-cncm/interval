![C/C++ CI](https://github.com/orlarey/interval/actions/workflows/ubuntu.yml/badge.svg)![C/C++ CI](https://github.com/orlarey/interval/actions/workflows/macos.yml/badge.svg)![C/C++ CI](https://github.com/orlarey/interval/actions/workflows/windows.yml/badge.svg)

# Intervals

The interval library of the Faust compiler. Since July 2026 it is *the* interval
domain of the compiler's type system: every signal's range is produced here, by the
fixpoint engine running over these algebras.

Licensed under Apache-2.0 (see LICENSE). Copyright 2020-2026 Yann Orlarey, Agathe
Herrou, Stéphane Letz.

## The two roles of the interval computation

- **Correctness** — the program runs right: delay lines sized from provable bounds,
  table accesses provably within `[0, size)` (clamps elided when proven, kept when
  not), no division by zero, no NaN, no infinity, compile-time error reporting.
- **Sound quality** — the precision of the computations: in fixed point (FPGA
  targets), a signal's format takes its integer bits (msb) from its range and its
  fractional bits from the lsb analysis, so a tighter interval converts directly into
  signal-to-noise ratio. See BACKGROUND.md and PRECISION.md.

## Two layers

1. **Ordinary intervals** (`interval_def.hh`, `interval_algebra.hh`): a triplet
   `<lo, hi, lsb>` — bounds plus least-significant-bit precision — with one
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
  wrapping. Generated C/C++ programs must be compiled with `-fwrapv` (GCC/Clang) —
  the option applies to the generated program, not to the Faust compiler itself.
- **Integer modulo** (both operands with `lsb >= 0`) follows C semantics: the result
  has the sign of `x`, its magnitude is at most `min(|x|max, m-1)` with `m` the
  largest divisor magnitude, and `x % y = x` whenever `|x|` stays below the smallest
  nonzero divisor magnitude. The float path models `fmod`, whose bound closes at
  `|y|`.

- **Floating-point operations** use `itv::programPrecision()`: 2 double by default,
  1 single, 3 quad, 4 fixed point. Double mode uses the directed reference
  calculations described below. The other modes retain their existing rules:
  notably, single-precision bounds are converted to the nearest float when an
  interval is constructed, including subnormal rounding and underflow to zero.
  Those rules are not yet a general conservative analysis of float programs.
- **Target-libm compensation** uses `itv::libmCompensation()`, enabled by default.
  The historical margin is two ULPs of the program precision, clipped to known
  function images. It is an assumption, not a universal bound for arbitrary
  libms. In double mode, the reference calculation uses the native kernel even when this
  compensation is disabled; computed singletons no longer bypass the margin.
  The other modes retain their previous singleton exemption and host-libm
  evaluation, which need a separate audit.

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
such as `sin(0)`, `exp(0)` and `log(1)` also remain points.

This is the first step towards a reliable computational analysis. It does not
certify arbitrary target libms or reassociation, add separate NaN/signed-zero
tracking, prove recurrence invariants, or certify the LSB precision estimates.
Those obligations need separate verification before bounds justify safety
decisions. The single/quad/fixed paths retain their previous rules and need
their own audit. The new transcendental reference is intended for inclusion;
its performance and tightness in the complete Faust compiler need measurement.

Six reproducible missing inclusions in the current float mode are documented in
[FLOAT_COUNTEREXAMPLES.md](FLOAT_COUNTEREXAMPLES.md), including an analyzed
nonnegative integer delay whose strict float execution produces `-1`.
`FloatConservativenessTests` requires inclusion by default and currently fails.
CTest labels these as `known-float-gaps` and uses an explicit expected-gap mode;
a passing known-gap test confirms the defect, not float conservativeness.

## Building and optional MPFR oracle

A normal build needs a C++20 compiler and CMake, with **no MPFR/GMP dependency**:

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
- `empty()` constructs the empty interval, represented internally with `NaN` bounds.
- `interval(-HUGE_VAL, HUGE_VAL)` also contains infinities and is therefore distinct from
  `fullFinite()`.

## Organization of the code

All the code lives in the namespace `itv`:

- `interval_def.hh` — the interval data structure and its basic accessors
- `interval_algebra.hh/cpp` — all operations on intervals, as defined by the Faust
  primitives; one `intervalXXX.cpp` per operation
- `bitwiseOperations.hh/cpp` — helpers for the bitwise operations
- `affint.hh` — the affine-in-time interval (`AffItv`) and its numeric core
- `affine_ops.hh` — `AffineOps<Base>`, the FaustAlgebra operations over `AffItv`;
  `affine_algebra` is its standalone instantiation
- `check.hh/cpp` — test helpers (exact checks and sampled numerical analyses)
- `FaustAlgebra/` — vendored copy of the FaustAlgebra interface
- `tests/` — the deterministic regression suite (see TESTING.md)
