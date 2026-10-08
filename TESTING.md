These notes on testing the intervals results are destined to future library contributors who might want to add tests to validate their additions to the library.

# Test suites

The project has five dependency-free executables and three optional oracle targets:

- `IntervalTests` contains fast, deterministic regression tests. It checks interval bounds and
  LSB values strictly and is the suite executed by CTest and CI.
- `TestInterval` is an experimental numerical-analysis program. It can run large sampled studies
  of interval operations and is not part of the default CTest suite.
- `DirectedBoundsTests` contains dependency-free rounding regressions, including
  half-ULP errors, false singleton underflow, overflow, exact constants, poles,
  affine interpolation, and preservation of the host rounding mode. It also runs
  under Emscripten/Node (WebAssembly has a fixed rounding mode).
- `MPFROracleTests`, built only with `INTERVAL_ENABLE_MPFR_TESTS=ON`, compares
  native bounds against independent 256-bit directed MPFR calculations. Fixed
  seeds cover all binary exponents and ordinary arguments for elementary and
  transcendental functions, exact modulo/remainder, powers, scaling and interval
  images. MPFR is linked to this test executable, never to the interval library.

- `FloatConservativenessTests` requires inclusion for the six original float
  counterexamples and the large-sine diagnostic. All seven cases must now pass;
  CTest uses the ordinary driver, with no expected-failure/known-gap mode.
- `FloatBoundsTests` checks literals versus computed points, mixed conversions,
  float nature (including zero and integral-valued functions), signed-zero
  infinities, UI domains, and the affine bridge. Fixed bit-pattern samples cover
  all binary32 exponents, point/nonpoint arithmetic and 4000 FMA contractions;
  sixteen unary functions also cover domain and overflow/underflow thresholds.
- `FloatConservativenessOracleTests`, enabled by the MPFR option, verifies the
  seven original witnesses at precision 24, rounding after each target operation.
  Their target results are normal or zero and need no subnormal emulation.
- `FloatBoundsOracleTests` adds independent 256-bit directed real references and
  correctly rounded binary32 references. It sets MPFR precision 24, exponent
  limits -148/128, and uses `mpfr_subnormalize` with each operation's ternary
  result to emulate gradual underflow without double rounding. Limits are restored
  before the wide reference or any interval-library call. See the
  [MPFR subnormal-emulation documentation](https://www.mpfr.org/mpfr-current/mpfr.html#Exception-Related-Functions).

See [FLOAT_COUNTEREXAMPLES.md](FLOAT_COUNTEREXAMPLES.md) for the original failures,
corrected bounds and reproduction commands. Both direct execution and CTest
require inclusion and return zero on success:

```sh
./build/FloatConservativenessTests --require-inclusion
./build/FloatConservativenessTests --case delay
./build/FloatConservativenessTests --case sine-large
ctest --test-dir build -L float-inclusion --output-on-failure
```

The optional witness executable has the same arguments. Invalid arguments return
2; the old `--expect-known-gaps` mode is removed. UBSan must also pass the large
sine case, which formerly exposed an infinity-to-int metadata conversion.

No external numerical library is required for the default tests. For the oracle,
install the MPFR/GMP development files described in README.md and add
`-DINTERVAL_ENABLE_MPFR_TESTS=ON` when configuring a separate build directory.
CTest runs the oracle when it is built. CI checks both dependency-free and oracle
configurations, and also executes the default tests in WebAssembly under Node.
An oracle comparison detects implementation errors; it does not replace the
analytic inclusion arguments and assumptions documented in README.md/source.

Configure, build and run the regression suite with:

```sh
cmake -S . -B build -DNOTIDY=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

Use `checkExact` when both the bounds and the LSB are part of the expected result. The historical
`check` overload for intervals compares bounds and reports an LSB difference as a warning; it is
kept for exploratory numerical analyses.

Sampled analyses use fixed seeds so a failure can be reproduced. They remain complements to,
not replacements for, deterministic boundary and regression tests.

The regression suite (`tests/interval_tests.cpp`) covers both layers: the ordinary
interval algebra (including the integer paths of Add and Mod) and the affine-in-time
layer (`affine_algebra`: rate lifecycle, delays, joins, widening, the int32 cap at the
`toItv` bridge).

# Testing functions

How the various test functions work.

## `check`

# Number of samples

In a number of cases, the number of samples needs to be higher than usual to give valid results.

## Large intervals
If the interval is larger than usual (for example, if we're testing the limits of number formats), the number of samples needed to have a dense covering of the interval is higher.
On these interval, the `analyzeUnaryMethod` and `analyzeBinaryMethod` functions can measure a precision that is coarser than the actual ground-truth precision, because it didn't draw enough samples to land on a minimal gap.
This should not be taken as an indicator that the implemented interval function is incorrect. 
Increasing the number of samples should solve the issue.

## Binary functions
A binary function has two input intervals, and thus the domain to cover with the sampling is quadratically larger than for a unary function with similar input interval.
Thus, the number of samples needed is the square of the number needed to have the same results for a unary function.
