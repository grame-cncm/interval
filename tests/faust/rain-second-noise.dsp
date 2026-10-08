// The second LCG stage used to lose integer nature in the affine domain,
// making its always-true predicate false and replacing the output by zero.
import("stdfaust.lib");
process = no.multinoise(2) : !, (_ <: @(1), (abs < 1.5) : *);
