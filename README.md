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

- **Floating-point operations** follow two settings of the user:
  - `itv::programPrecision()`, the precision of the program: 2 double by default,
    1 single, 3 quad, 4 fixed point. The Faust compiler declares it from `-single`,
    `-double`, `-quad`. In single precision, every bound of a float-carried interval
    (`lsb < 0`) is rounded to the nearest float when the interval is built
    (`programBound`). Round to nearest is monotone, and for `+`, `-`, `*`, `/` and
    `sqrt` of floats the double rounding is innocuous (53 >= 2×24 + 2: rounding in
    double, then in float, gives the float the program computes, even when the double
    result is not exact): an elementary operation gives the very bound the program
    computes, and a constant stays a point (`0.1 + 0.2` is the float `0.3f`). Not
    rounded: the integer bounds beyond 2^24 (an integer value may carry a float
    precision by default). Subnormal bounds are rounded too, including underflow
    to zero; `pow` uses nonzero domain separators representable at the program's
    precision.
  - `itv::libmCompensation()`, true by default. IEEE 754 requires a correct rounding
    of `+`, `-`, `*`, `/` and `sqrt` only: the functions of the libm (`sin`, `cos`,
    `tan`, the inverse and hyperbolic functions, `exp`, `log`, `log10`, `pow`) may be
    off by an ulp (macOS: `sinf(1.00000596f)` is one ulp above `(float)sin(1.00000596)`),
    and the program calls the libm of its target, not the libm of the machine that
    computes the intervals. Their bounds widen by 2 ulps of the program's precision
    (`libmBounds`), within the image of the function (`sin` stays in `[-1, 1]`, `exp`
    stays `>= 0`), never for a point (a function of a constant is folded at compile
    time), an exact 0 (the C standard makes the libm exact there, annex F) or an
    integer result. Turn it off only when both libms are correctly rounded (CORE-MATH
    for instance): the libm of the target and the libm of the machine that computes the
    intervals.

  Not covered by the bounds: the FMA contraction and the reassociation of the C++
  compiler. The decisions that read the intervals keep a margin for those: the
  compiler keeps the guard of a table access whose index is within 4 ulps of an edge.
  The rules that reason on reals (the hull of a convex combination, where
  `(1 - t)*m + t*m` is `100.0` in float for `m = 99.9999924`) add
  `ulpMargin(lo, hi, k)`, k ulps of the program's precision at the magnitude of the
  bounds.

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
