# Coverage of the Faust primitives

`interval_algebra` implements the complete `FaustAlgebra<interval>` interface — every
primitive of the Faust signal language has an interval interpretation. The list below
is the current inventory (from `interval_algebra.hh`); each operation `XXX` is
implemented in `intervalXXX.cpp`, except the structural ones grouped in
`interval_algebra.cpp`.

## Numbers, casts, external values

`IntNum`, `Int64Num`, `FloatNum`, `Label`, `IntCast`, `FloatCast`, `BitCast`,
`ForeignConst`, `ForeignVar`, `ForeignFunction`

Foreign values have no computable bounds: they return the full finite interval
(`interval()`), which is sound and join-neutral — never the empty interval.

## Arithmetic

`Add`, `Sub`, `Mul`, `Div`, `Neg`, `Inv`, `Abs`, `Min`, `Max`, `Mod` (integer path +
`fmod` path), `Remainder`, `Pow`, `Rint`, `Round`, `Floor`, `Ceil`

Integer `Add` models int32 wrapping (see README). `Mod` has a dedicated integer path
with C semantics.

## Comparisons and logic

`Eq`, `Ne`, `Gt`, `Ge`, `Lt`, `Le`, `And`, `Or`, `Xor`, `Not`, `Lsh`, `ARsh`, `LRsh`

## Transcendental

`Exp`, `Exp10`, `Log`, `Log10`, `Sqrt`, `Sin`, `Cos`, `Tan`, `Asin`, `Acos`, `Atan`,
`Atan2`, `Sinh`, `Cosh`, `Tanh`, `Asinh`, `Acosh`, `Atanh`

## Time and memory

`Mem`, `Delay`, `Prefix`, `FixPointUpdate`

## Tables, waveforms, soundfiles

`Gen`, `WRTbl`, `RDTbl`, `Waveform`, `SoundFile`, `SoundFileLength`,
`SoundFileRate`, `SoundFileBuffer`

## User interface and control

`Button`, `Checkbox`, `HSlider`, `VSlider`, `NumEntry`, `HBargraph`, `VBargraph`,
`Attach`, `Enable`, `Control`

## Selection, bounds, I/O

`Select2`, `AssertBounds`, `Highest`, `Lowest`, `Input`, `Output`

## The affine layer

`AffineOps<Base>` (affine_ops.hh) lifts every operation above to affine-in-time
intervals: linear operations act on the coefficients, nonlinear ones collapse to the
hull over `[0, T]` and delegate to `interval_algebra` as the oracle, mixed cases
(multiplication or division by a rate-0 operand) chord the oracle's endpoint values.
