# The interval of mix

`mix(c, a, b) = (1 - c) * a + c * b` is the interpolation of a smoother. This note
explains how its interval is computed (`numericMix` and `mixEnvelope` in
`interval/intervalCompound.cpp`), and in particular how a recursive smoother
`y = mix(c, x, y@1)` whose weight `c` is a few rounding errors below 1 is bounded.

## What the program computes

The program interpolates from the operand whose weight may be near 1, chosen on
the interval of `c` (`mixFromB`):

| interval of `c` | computation |
| :-- | :-- |
| `c >= 1/2` everywhere | `b + (1 - c) * (a - b)`, where `1 - c` is exact (Sterbenz) |
| otherwise | `a + c * (b - a)` |

An interpolation keeps its fixed point exactly (`a = b` gives `a`), and costs about
one rounding of the result when the weight of the difference is small.

## Three bounds

When `c` lies in `[0, 1]` and `a`, `b` are bounded, the interval of the mix is the
intersection of three bounds:

1. **The hull with a rounding margin.** For a fixed `c` the mix is linear in `c`,
   so the extremes of the exact mix are reached at the two ends of `c`. They are
   widened by the rounding of the program: about `2u w W + u |r|`, with `u` the
   unit roundoff, `w` the largest weight of the difference, `W` the largest
   `|b - a|` and `|r|` the magnitude of the result.
2. **The sign.** A mix of nonnegative operands stays nonnegative, whatever the
   rounding (and symmetrically).
3. **The power-of-two envelope.** `[-M, M]`, with `M` the smallest power of two
   containing `a` and `b`, as stated below. It applies in float and double, away
   from the limits of the format: `2^-102 <= M <= 2^126` in float,
   `2^-969 <= M <= 2^1022` in double.

Outside of these premises (`c` beyond `[0, 1]`, an unbounded operand), the
interval is the one of the definition, computed operation by operation.

## The power-of-two envelope

**Theorem.** In a binary floating-point format with `p` bits of significand,
rounding to nearest with ties to even, let `M` be a power of two, at least the
smallest normal number, with `2M` representable. For all floats `a, b` in
`[-M, M]` and `c` in `[0, 1]`, the computed `a + c * (b - a)` lies in `[-M, M]`,
with the product and the sum rounded separately or fused (FMA). So does
`b + (1 - c) * (a - b)`.

**Proof.** Rounding `fl` is nondecreasing and odd (`fl(-t) = -fl(t)`).

*Upper bound.* Let `d = fl(b - a)` and `D = fl(M - a)`. Since `b <= M`, `d <= D`.
The product does not exceed `D`: if `d >= 0`, `c d <= d` and `d` is a float, so
`fl(c d) <= d <= D`; if `d < 0`, `fl(c d) <= 0 <= D`. By monotonicity the result is
at most `fl(a + D)`, fused or not. Let `t = M - a`, in `[0, 2M]`. If `t` is a
float, `a + D = M`. Otherwise `t < 2M`, its ulp is at most `2uM` (the ulp of
`[M, 2M)`), and `D <= t + uM`, so `a + D <= M (1 + u)`. The float after `M` is
`M (1 + 2u)`: `M (1 + u)` is their midpoint. A value below the midpoint rounds to
`M`, and so does the midpoint itself, `M` being even (significand `1.00...0`).

*Lower bound.* Rounding is odd: the computation on `(-a, -b)` gives exactly the
opposite result, and the upper bound applies to it.

*The other form.* `b + w (a - b)` with `w = fl(1 - c)` in `[0, 1]` is the same
interpolation with `a` and `b` exchanged. ∎

The power of two is necessary: it puts the possible overshoot exactly on a
midpoint, that ties-to-even brings back to `M`. In float, with `a` the float just
above `-12`, `b = 12` and `c = 1`, the difference `23.999999` rounds to `24`, and
the result is `12 + 2^-20`. With `b = M = 16`, no `a` in `[-16, -4]` gives more
than 16.

## Why a recursive smoother needs it

The compiler bounds a recursion by assumption and check: it assumes `y` in `Y`,
computes one step `mix(c, x, Y)`, and accepts `Y` when the step stays in `Y`. The
candidates form a ladder of powers of two (`[0, 1]`, `[-1, 1]`, `[0, 16]`,
`[-16, 16]`, ...). With the first bound alone the check compares the pull of one
step towards `x`, `(1 - c)(M - 1)`, with the margin, which grows like `u M`. When
`1 - c` is only a few rounding errors the margin wins, for every `M`. Measured in
float, the envelope disabled, for `x` in `[-1, 1]`:

| assumption `Y` | `c` in `[0, 1 - 4u]` | `c` in `[0.4, 1 - 2u]` |
| :-- | :-- | :-- |
| `[-1, 1]` | `1 + 6u` | `1 + 6u` |
| `[-16, 16]` | `16 + 32u` | `16 + 64u` |
| `[-256, 256]` | `256 + 512u` | `256 + 1024u` |

The envelope cuts each step back to `[-M, M]`: the first rung `[-1, 1]` maps into
itself and the recursion is bounded. The weight `[0, 1 - 4u]` is the one computed
for `exp(-6.91 / (SR * T60))` with `SR <= 192000` Hz and `T60 <= 100` s, the
smoother of the SAM echo examples, which are accepted in float.

## Tests

`tests/compound_tests.cpp`, in float and double:

- both forms, plain and fused, stay in `[-M, M]` for `M` in `{1, 16, 1024, 2^-60}`,
  over 200000 draws per `M`, the extremes `±M` included and weights taken near 0,
  1/2 and 1;
- a smoother of weight `[0.4, 1 - 2u]` maps `[-1, 1]` into itself;
- a mix of operands in `[-12, 12]` stays in `[-16, 16]`;
- the interval contains the computed value on the cases of the rounding margin.

## The recursive mix `Irm`

`Irm(x, c)` is the signal `y` with `y(n) = mix(c(n), x(n), y(n - 1))` and
`y(-1) = 0`: a smoother of `x` whose pole is `c`. Its interval is built from an
*invariant*: an interval `Y` that holds `0` and such that `Mix(c, x, Y)` lies in
`Y`. By induction every `y(n)` then lies in `Y`, rounding included. The rounding
is the one of either form of the interpolation, from `x` with the pole `c` or
from the state with the weight `1 - c` of `x`: a code generator may choose one or
the other (the weight of `x` is sometimes exact in the program, the pole computed
from it is not). The samples are then mixes of `x` and of a state in `Y`: the
result is `Y` met with `Mix(c, x, Y)`, which leaves out the initial state `0` (a
smoother of a positive `x` is positive).

Two invariants are tried, and their intersection, itself an invariant (the mix of
intervals is monotone), is the result:

1. **The hull of `x` and `0`, with a margin.** The exact mix of a value of `Y` and
   a value of `x` moves a bound of `Y` inwards by `1 - c` times its distance to
   `x`. That pays for the rounding of the mix once the margin reaches about one
   rounding divided by `1 - c`. The margin starts at one rounding of the result
   and doubles until the mix maps `Y` into itself; a bisection then narrows it to
   within an eighth of the last doubling. A bound at `0` that `x` does
   not cross needs no margin (the sign of the mix).
2. **The power-of-two envelope** of `x` and `0`, which needs no margin and holds
   for a pole up to `1`.

| `x` | `c` | `Irm(x, c)` |
| :-- | :-- | :-- |
| `[0, 384000]` | `[0.99, 0.9999]` | `[0, 384000]` plus a margin below 1 % |
| `[-1, 1]` | `[0, 1]` | `[-1, 1]` |
| `[-3, -2]` | `[0.5, 0.9]` | about `[-3, -0.2]`, the lower bound with a small margin |
| `[1, 20]` | `[0.9, 0.99]` | about `[0.01, 20]` |

A pole that may leave `[0, 1]` (a recursion that may diverge) or an unbounded `x`
gives anything, NaN included.
