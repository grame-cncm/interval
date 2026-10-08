// Expected: reject this program, rather than silently clamp its write index.
// With the modulo outside the recurrence, the counter eventually wraps int32.
// Its full int32 range gives a signed remainder in [-199, 199], then an index
// in [-249, 149]. The explicitly bounded read must not hide the write error.
i = (+(1) ~ _) : %(200) : -(50);
process = rwtable(100, 0.0, i, _, max(0, min(i, 99)));
