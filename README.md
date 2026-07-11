![C/C++ CI](https://github.com/orlarey/interval/actions/workflows/ubuntu.yml/badge.svg)![C/C++ CI](https://github.com/orlarey/interval/actions/workflows/macos.yml/badge.svg)![C/C++ CI](https://github.com/orlarey/interval/actions/workflows/windows.yml/badge.svg)

# Intervals

A simple interval library for the Faust compiler. 

## Purposes of intervals as used by the Faust compiler

- determine the size of a delay line
- check that a Faust program will execute correctly, without producing infinite or NAN values
- compute the precision needed for fix-point and fpga implementations
- detect and report errors at compile time

## Definition of an interval

An interval is a triplet <lo,hi,lsb> where:
- lo: is the lowest value
- hi: is the highest value
- lsb: is the least significant bit
- The msb can be deduced from lo and hi

An interval represent integer values if lo and hi are integers and if lsb >= 0

## Special intervals

- `interval()` is the historical default interval. It contains every finite value
  representable by a `double` and has an LSB of `-24`.
- `fullFinite()` constructs the same interval explicitly, optionally with a different LSB.
- `empty()` constructs the empty interval, represented internally with `NaN` bounds.
- `interval(-HUGE_VAL, HUGE_VAL)` also contains infinities and is therefore distinct from
  `fullFinite()`.

## Integer overflow semantics

Integer interval operations model signed 32-bit arithmetic with two's-complement wrapping.
Generated C/C++ programs must be compiled with `-fwrapv` when using GCC or Clang. The option
must apply to the generated program itself, not only to the Faust compiler.

## Organization of the code

All the code is encapsulated in the namespace 'itv'. It is organized as follows:

- interval_def.hh : defines intervals as data structures with some very basic methods to access the fields
- interval_algebra.hh/cpp: class gathering all operations on intervals as defined by Faust primitives.
- intervalXXX.cpp: implementation of the XXX operation on intervals.
