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

## Special intervals

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
