The `interval` library was initially developed in the context of the FAST project,
which aims at porting the Faust language to FPGA. The fixed-point precision material
below (msb/lsb formats, pseudo-injectivity, forward and backward propagation) comes
from that work and remains the reference for the `lsb` field carried by every
interval.

*Update, July 2026:* the library has since become the interval domain of the Faust
compiler's type system itself — every signal range is produced by a fixpoint over
`interval_algebra` — and was extended with affine-in-time intervals (`affint.hh`,
`affine_ops.hh`) that date a program's failure modes (int32 wrap, float absorption)
instead of merely bounding its values. The fixed-point role described here is
unchanged: it is one of the two roles of the computation (see README.md).

# Faust on FPGAs

Faust is a functional programming language for signal processing, typically used to create audio applications such as virtual instruments, or plugins for DAWs. 

FPGAs are reprogrammable logical circuits. 
Because a program on an FPGA is running on *"bare metal"* without being mediated by an operating system, very high performance and low latency are possible.
This is particularly valuable in audio applications where latency can be critical, and where very low latency allows for almost instantaneous processing of the signal.

Programming on FPGA is normally done using a low-level Hardware Description Language (HDL) such as VHDL or Verilog. 
However, programming in these languages is not an easy task, which is why a common approach to programming FPGAs is High Level Synthesis (HLS), whereby the program to run on the FPGA is described in a higher-level language and then automatically translated to a HDL.

This is the motivation behind the FAST project, whose goal is to allow programming on FPGA for audio applications using the Faust language.
Compilation of Faust programs for FPGA is done with the Syfala toolchain[^2].

# Fixed-point precision in Faust programs

When programming on FPGAs, using fixed-point numbers is much more efficient than using floating-point numbers, both in terms of computation time and resources used.

To ensure numerical precision, we need to determine, for each operator and wiring present in the final design, the number of bits used for the numbers that they manipulate. 
This number of bits (as well as the position of the binary point) is referred to as the *fixed point format* of the numbers.

## Fixed-point formats

A fixed-point format is defined by two parameters: the Most Significant Bit (MSB), which we will denote $m$ in this document, and the Least Significant Bit (LSB), denoted $l$. 
These parameters inform us of the positions between which a given number is guaranteed to have non-null bits.

For example, `101.1101` has MSB 2 and LSB -4.

These notions can be extended to signals (i.e. infinite sequences of numbers), where the MSB and the LSB will respectively be the maximum MSB and the minimum LSB of any number appearing in the signal.

Determining the MSB is fairly straightforward: it requires determining the order of magnitude of a signal, i.e., the maximum absolute value it is susceptible to take.
This can be done by determining the interval in which the signal is guaranteed to be constrained, through application of the rules of interval algebra.

Determining the LSB is a more subtle task. 
On the one hand, we would like to minimise the use of resources, which means using the coarsest precision we can afford (i.e., the LSB with the highest value), but on the other hand, this should not be done at the cost of signal quality. 
The property we have chosen as criterion to link the input and output LSBs is the pseudo-injectivity property defined as follows.

**Notation**: $u = 2^l$ denotes the gap between two consecutive fixed-point numbers with precision $l$.

**Notation**: $\lfloor x \rfloor_l = u \cdot \lfloor \frac{x}{u} \rfloor\$ denotes the truncation of the real number $x$ up to $l$ bits.

**Notation**: $[a;b]_l$ denotes the interval $[a;b]$ with LSB $l$: $[a;b]_l = [a;b]\cap u\mathbb{Z}$

**Pseudo-injectivity property**:
We say that the fixed-point representation of a unary, real function $f : [a; b] \rightarrow \mathbb{R}$ with input LSB $l_{in}$ and output LSB $l_{out}$ respects the pseudo-injectivity property if
$\forall x, y \in [a;b], \lfloor x \rfloor_{l_{in}} \neq \lfloor y \rfloor_{l_{in}} \Rightarrow \lfloor f(x) \rfloor_{l_{out}} \neq \lfloor f(y) \rfloor_{l_{out}} \vee f(x) = f(y)$.

The condition $f(x) = f(y)$ serves to capture the case where function $f$ is not injective and there are legitimate distinct pairs $x \neq y$ that yield the same image. 
These pairs should not be considered for testing pseudo-injectivity, as their equality in fixed-point comes not from a truncation artefact but from a real mathematical equality.

This property can also be formulated as its contraposition:
$\forall x, y \in [a;b],  \lfloor f(x) \rfloor_{l_{out}} = \lfloor f(y) \rfloor_{l_{out}} \wedge f(x) \neq f(y)  \Rightarrow \lfloor x \rfloor_{l_{in}} = \lfloor y \rfloor_{l_{in}}$.

In more intuitive terms, a function is pseudo-injective for a pair of input and output fixed-point formats if all distinguishable outputs in the real function correspond to distinguishable outputs in its fixed-point counterpart.

A property that will show useful when propagating determined LSBs is the composability of the pseudo-injectivity property.

**Composability of the pseudo-injectivity property**: 
If $f: [a;b]\_{l\_{in}} \rightarrow [a';b']\_{l\_{mid}}$ and $g: [a';']\_{l\_{mid}} \rightarrow [a'';b'']\_{l\_{out}}$ are pseudo-injective, then $g\circ f:[a;b]\_{l\_{in}} \rightarrow [a'';b'']\_{l\_{out}}$ is pseudo-injective as well.

**Proof**
Let $f$ and $g$ be such pseudo-injective functions.
Let $x, y \in [a; b]\_l$ such that $\lfloor x \rfloor\_{l\_{in}}\neq \lfloor y\rfloor\_{l\_{in}}$.

Let us first assume that $f(x) \neq f(y)$. 
Then, by pseudo-injectivity of $f$, $\lfloor f(x) \rfloor_{l_{mid}} \neq \lfloor f(y) \rfloor_{l_{mid}}$.
If we also assume that $g\circ f (x) \neq g \circ f(y)$, by pseudo-injectivity of $g$ we can deduce that $\lfloor g\circ f (x) \rfloor_{l_{out}} \neq \lfloor g\circ f (y) \rfloor_{l_{out}}$, which verifies the right hand of the implication of the pseudo-injectivity of $g\circ f$.

In the case where $f(x) = f(y)$, $g\circ f(x) = g\circ f(y)$, which verifies the other part of the disjunction.

In the case where $f(x) \neq f(y)$ but $g\circ f(x) = g\circ f(y)$, it is verified it as well.
□

Thus, if each of the individual components of a circuit are pseudo-injective, by composition, the circuit as a whol is pseudo-injective as well.

# Determining formats

## Format propagation

In order to determine all the formats in a given circuit, Beata Burreau[^1] suggests that the associated signal graph is scanned using a traversal algorithm.
This first scan computes, for each node (representing a function or an operator in the Faust program), its output LSB given the input LSB(s). 
For that purpose, the algorithm is designed in such a way that a node is not treated before all its predecessors have been treated.

This gives us a final output LSB for the sinks of the program.
This LSB can exceed the precision called for by the output (typically a LSB of -24): this means that at least some of the computed internal LSBs need not be so precise, and that we can thus trim some of them by propagating the output LSB backwards.

A second scan is thus performed, starting from the outputs and working its way through the signal graph towards the inputs. 
This calls for a procedure capable of inferring the input LSBs of an operator given its output LSB. 
When the inferred LSB from the backwards propagation is coarser than the previously computed value, it is updated.
Note that a finer LSB cannot possibly happen because LSB inference is monotonic.

Issues can arise with this approach when the signal graph presents a cycle due to the presence of a feedback loop: in the case where an operator that has higher output LSB than input LSB is looped on itself, the prescribed LSB might get finer and finer at each computation and never settle on a finite value.
In addition, even in the case where a sound LSB value is determined, there is a real risk of accumulation of rounding errors in the loop.
One proposed solution to deal with this is to keep a floating-point format for operators caught in such loops.

Thus, we need, for each potential operator we might encounter, to determine how to propagate LSBs both from the input to the output, and from the output to the input, while respecting the pseudo-injectivity property. 

## LSB computation

To determine the output LSB of a unary function, given an input LSB and an input interval, we need to determine how close two outputs can get.
In the case of a continuous, unary, monotonic function, we can restrict our study to that of consecutive outputs: in that case, the closeness of such outputs is directly correlated to the slope of the 
function at that point, i.e. to the amplitude of its derivative.
See PRECISION.md for details of the computations.

In the forward direction, we search for the coarsest precision that will preserve the pseudo-injectivity condition: two distinguishable input numbers should produce two distinguishable outputs.

In the backwards direction, we search for the finest precision that will preserve the reverse pseudo-injectivity condition: two indistinguishable outputs should have been produced by two indistinguishable inputs.

## Testing the precision

Once a theoretical precision bound has been determined, we need to test it empirically.
 
In the forward direction (i.e. when propagating the LSB from inputs to outputs), this is done by randomly polling the discrete input interval, keeping an ordered record of the samples' images through the function tested while discarding equal values. 
This record is then scanned to determine how close two images can get and determine the corresponding LSB. 
The precision is obtained from these two closest images by computing their longest common prefix and taking `log2` of the result in order to count its bits.
This empirical LSB is expected to be higher (i.e. coarser) or equal to the theoretical LSB (assuming it has been correctly determined): it is equal if two samples corresponding to the smallest gap are drawn, higher otherwise.

An interesting question to explore would be to try to determine what gap between theoretical and empirical precisions can be expected, given the number of samples drawn.

The backwards direction has not been implemented yet, but an algorithm has been imagined as follows:

To empirically test the backwards direction, we poll the input interval (considered continuous and represented by floating-point numbers) in a way that the sampling order respects the order on the interval, and compute their images.
We group preimages that have indistinguishable images into bins, and for each bin compute their longest common prefix.
This bin-wise prefix corresponds to a canonical representative of the bin.
Its length yields a precision that corresponds to the finest precision we can "afford" to give to the inputs in the bin without jeopardising the pseudo-injectivity property.
For the final result over the interval, we take the maximum over all these bin-wise precisions.

### Detail of the backwards polling method

We want to draw $n$ samples $(x\_i)\_{1≤i≤n}$ such that $x\_i < x\_{i+1}$ for all $i$ over the span of interval $[a;b]$.

Unlike for forward polling, we don't draw the samples directly, but instead start from the lowest part of the interval $x\_0 = a$ and draw the step $δ\_i = x\_{i+1} - x\_i$ to get from this sample to the next. 
This step is drawn with a gaussian distribution with mean $μ = \frac{b-a}{n}$, so that we will cover the whole interval on average in $n$ steps. 
We take the standard deviation to be $σ = μ$ for no particularly good reason. We'll probably have to adjust that empirically. 
If the step $δ_i$ we have drawn is $≤ 0$, we reject it, since we want the samples to be in strictly increasing order.

## Fuzzy precision

**TODO** What is done in `analyze-precision`.

## Expanding to non-unary functions

**TODO** e.g. pow, atan2, division

[^1]: https://raw.githubusercontent.com/beataburreau/FAUST-signal-graphs-goes-Haskell/main/ProjectReport.pdf

[^2]: https://github.com/inria-emeraude/syfala
