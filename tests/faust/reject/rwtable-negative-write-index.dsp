// Expected: reject this program, rather than silently clamp its write index.
// The modulo is inside the recurrence: the counter stays in [0, 199].
// Its shifted write index is in [-50, 149], outside the table's [0, 99].
// The read index is explicitly bounded; that must not hide the write error.
i = ((+(1) : %(200)) ~ _) : -(50);
process = rwtable(100, 0.0, i, _, max(0, min(i, 99)));
