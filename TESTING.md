These notes on testing the intervals results are destined to future library contributors who might want to add tests to validate their additions to the library.

# Test suites

The project has two complementary executables:

- `IntervalTests` contains fast, deterministic regression tests. It checks interval bounds and
  LSB values strictly and is the suite executed by CTest and CI.
- `TestInterval` is an experimental numerical-analysis program. It can run large sampled studies
  of interval operations and is not part of the default CTest suite.

Configure, build and run the regression suite with:

```sh
cmake -S . -B build
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
