This document describes how we compute the output precision of a function over an interval with certain input precision.

# General principle

We consider a function $f$ that we study over the interval $[\underline{x};\overline{x}]$, with $\underline{x} < \overline{x}$. The fixpoint format for the input is given, we denote its lsb $l$ and $u = 2^l$ the gap between two consecutive numbers. The goal is to determine the minimum of $|f(x) - f(y)|$ with $x$, $y$ ranging over $[\underline{x};\overline{x}]$, which will inform us of the minimum number of digits needed to distinguish between all the images of $f$ over $[\underline{x};\overline{x}]$.

This minimum number of digits needed to distinguish between two numbers $x$ and $y$ is linked to the notion of their Least Common Bit (LCB), that is the position of the first bit at which they differ.
If two non-negative numbers $x$ and $y$ have LCB $-k$, then $|x-y|\le 2^{-k}$ and $log_2|x-y| \le LCB(x, y)$: thus, $log_2|x-y|$ is a sound approximation of their LCB.
This is particularly useful in situations where we do not have direct access to the digits of the numbers, for example when studying them in a general abstract setting.

In the case of monotonous functions, this minimum will be attained for two consecutive arguments: $|f(x) - f(x±u)|$. This $x$ is determined as the point where function $f$ has lowest slope over $[\underline{x};\overline{x}]$.

When a function is not monotonous, it can happen that two non-consecutive fixpoint arguments have images closer than any two consecutive fixpoints. The usual functions subjected to this phenomenon are the periodic trigonometric functions $sin$, $cos$ and $tan$. For these functions, we relax the constraint and only seek the lowest difference between two consecutive arguments.

The overall goal is to find the proper $x$ and $±u$ for each usual function $f$. After that, the target lsb is given by $l' = prec(f, x, ±u) = ⌊log₂(|f(x±u) - f(x)|)⌋$.

## Possible pitfalls

If the input precision $l$ is too fine, then $f(x+u)$ and $f(x)$ might be so close that they will be computed as the same floating-point number. In that case, the computation of $f(x+u) - f(x)$ will yield zero and the computed output precision $-\infty$. 

When a situation of this type occur, we fall back on another method for precision computation, based on a Taylor expansion of $f$. We approximate $|f(x+u) - f(x)|$ with $u\cdot|f'(x)|$ and convert it into the output precision $l' = l + \log_2(|f'(x)|)$.

## Implementation conventions

The point at which the precision is computed, called $x$ above, is called `v` in the implementation. 
Whether the next point of the computation $x \pm \epsilon$ is the one before or after $x$ is encoded in `sign`. This mostly depends on whether $x$ is the higher or lower bound of the interval: we chose `sign` so that $x \pm u$ remains in the interval.

# Faust primitives

Here we explain the logic behind the implementation of precision inference in each Faust primitive. The interval $[\underline{x}; \overline{x}]$ is the input interval of the function, possibly restricted to its definition domain.

## Abs

abs is a function that is equal to the identity function on positive values and its opposite on negative values. As such, its output precision should be the same as its input precision.

### Integers

When applied to integers, Abs is implemented using the `std::abs` C++ function. 
This function has undefined behaviour when applied to `INT_MIN`, since `-INT_MIN = INT_MAX + 1` is not representable.
We chose to interpret this fact as the indication that `INT_MIN` is out of bounds for Abs. 
We deal with it the same way we deal with other domain violations: 
we restrict the input interval to legal values and apply Abs to the restricted interval.


**TODO**: Modifications will have to be made to the code generator as well to warn about domain violations.

## Acos 

acos is the reciprocal of the cosine function.

It is defined over the interval $[-1;1]$, so we restrict the input interval by intersecting it with this domain. If this intersection is empty, the output has no valid value and thus the output interval is empty as well.

Its derivative attains a global minimum at point 0 (and diverges towards $\pm \infty$ at the bounds of the interval). If 0 is included in the interval $[\underline{x}; \overline{x}]$, we compute $prec(f, 0, ±u)$. Otherwise, if $0 < \underline{x}$, the minimum derivative is at $\underline{x}$ and we compute $prec(f, \underline{x}, +u)$, and if $0 > \overline{x}$, the minimum derivative is at $\overline{x}$ and we compute $prec(f, \overline{x}, -u)$.

If the precision computed with this method is not sound (ie, the inferred precision is infinite), then we fall back on a less accurate precision determination method, based on the Taylor expansion of $acos$:
$|acos(x + u) - acos(x)| = |u\cdot\frac{-1}{\sqrt{1 - x^2}}|$.

We take the $\log_2$ of this expression to retrieve the relevant number of bits:
$l' = l - \frac{1}{2}\log_2(1 - x^2)$.

## Acosh

acosh is the reciprocal of the cosh (hyperbolic cosine) function.

It is defined over the domain $[1; + \infty[$: we restrict the input interval by intersecting it with this domain. If the intersection is empty, the output interval is empty as well.

It is a concave function : its derivative is decreasing, so the lowest slope is attained at the high end of the interval. We compute $prec(f, \overline{x}, -u)$, with $-u$ is in order to always evaluate $f$ at a point that is in the interval.

If this method computes an infinite precision, we fall back on a Taylor-based method: 
$|acosh(x + u) - acosh(x)| = |u\cdot \frac{1}{x^ -1}|$,
so $l' = l - \frac{1}{2}\log_2(\overline{x}^2 - 1)$.

## Add 

Addition is a binary operator, and the above method, designed for unary, real functions, does not apply to it. We have to devise another precision computation technique. 

Let us consider the addition of arguments $x$ and $y$, with respective precisions $l_x$ and $l_y$. Without loss of generality, let's assume that $l_x < l_y$, ie, $x$ is represented with more precision than $y$. In order to be able to distinguish the two sums $x + y$ and $x' + y$, where the operands $x$ and $x'$ only differ by their last bit, bits up to $l_x$ have to be conserved in the output.
Thus, we establish that $l' = min(l_x, l_y)$ is the coarsest output precision that still respects pseudo-injectivity.

### Integers

Wrapping situations due to integer overflow have to be taken care of. In the case where both operands are integers, if their sum is higher than the biggest representable integer or lower than the smallest representable integer, the default behaviour is for the value to "wrap" around and get to the other end of the integer range. 

The boundaries of the output interval are computed as `double`s in order to test whether the output interval wraps. 

If `INT_MAX` and `INT_MAX+1`, or `INT_MIN-1` and `INT_MIN` are both present in the double-value output interval, then the output interval contains both wrapped value and unwrapped values, and thus spans `[INT_MIN; INT_MAX]`.

Otherwise, the bounds of the output interval are computed by adding the bounds of the input intervals as integers, with potential wrapping.

**Note:** The gap between `INT_MAX` and the double right below is smaller than one: all values in an integer interval are representable as doubles.

## And 

And is a bitwise boolean operator.

Example: `And(1011, 1100) = 1000`. Only the positions where both operands have bits set to `1` get set to `1` in the result.  

`And` behaves like a mask: one of the arguments will "mask" the other in that only the bits where the first argument has a `1` will be preserved in the other, while the remaining bits will be indiscriminately set to `0`.

The correct precision to use then depends on whether we want to prioritise pseudo-injectivity or the concrete use case of bitwise boolean operators.

### Bitwise boolean operator

If we consider that all bits are significant, including the trailing zeroes, then the output precision precision should be the finest of both input precisions, floored to $0$ as the output is an integer.

This solution is the one that is used in the code, but the alternative is implemented in comments.

### Pseudo-injectivity

If we prioritise pseudo-injectivity (i.e., consider that we can discard bits as long as they're not useful to distinguish outputs), then the reasoning goes as follow:

Considering that all bits beyond the LSB of an argument have value `0`, all output bits in these positions will have value `0` as well. 
From the point of view of pseudo-injectivity, this means that these bits cannot be used to distinguish outputs, and can be ditched as such. 
As a consequence, the output LSB will coarser or equal than the coarsest input LSB: $l' \ge max(l_x, l_y)$.

In the case where neither of the intervals are singletons, ie when they both containt at least two values, this bound is attained.
Indeed, these intervals are guaranteed to contain values with the last bit set to both `0` and `1` values. 
This means that the "mask of 0s" will not expand prior to the LSB.

*Example:* $[000; 010]_0 \wedge [100; 110]_1 = \{000; 001; 010\} \wedge \{100; 110\} = \{000; 010\} = [000; 010]_1$

If one (or both) of the interval is a singleton, the mask can extend prior to the LSB. 
We detect the position of bit with value `1` and the least weight, and treat it as the "actual" LSB of the singleton-interval. 
The output LSB is the maximum of both input LSBs. 

*Example:* $[000; 110]_0 \wedge [100]_1 = \{000; 001; 010; 011; 100; 101; 110\} \wedge {100} = \{000; 100\} = [000; 100]_2$

Even though the input interval are specified with LSBs $0$ and $1$, the actual output interval can be represented with LSB $2$ without loss of information.

**TODO** Explain the bounds of the interval.

## Asin 

asin is the reciprocal of the sine function. It is very similar to the acos function in its analysis.

It is defined over the interval $[-1; 1]$: we restrict its input interval to this domain. Applying the function to an empty interval yields an empty interval.

Its derivative has a minimum in $0$: this is where we compute the precision if $0$ is included in the interval. Otherwise, the boundary of the interval closest to $0$ is used: $\underline{x}$ if $0 < \underline{x}$, $\overline{x}$ if $0 > \overline{x}$.

In case of an infinite precision, the Taylor formula is:
$|asin(x\pm u) - asin(x)| = |u \cdot \frac{1}{\sqrt{1- x^2}}|$, 
hence $l' = l - \frac{1}{2}\log_2(1- x^2)$.

## Asinh

asinh is the reciprocal of the sinh (hyperbolic sine) function. 

It is defined over the whole $\mathbb{R}$ set.

It is an increasing function, its derivative tends to 0 in both $+\infty$ and $-\infty$. Thus, the precision is computed at the point of the interval of highest magnitude: $\overline{x}$ if $|\overline{x}| > |\underline{x}|$, $\underline{x}$ otherwise.

The Taylor fallback is:
$|asinh(x+u) - asinh(x)| = |u \cdot \frac{1}{\sqrt{1 + x^2}}|$,
so $l' = l - \frac{1}{2} \log_2(1 + x^2)$.

## Atan

atan is the reciprocal of the tangent function. 

It is defined over the whole $\mathbb{R}$ set.

It is an increasing function, whose derivative tends to 0 in $+\infty$ and $-\infty$. Thus, its precision is computed at the point of highest magnitude: $\overline{x}$ if $|\overline{x}| > |\underline{x}|$, $\underline{x}$ otherwise.

The Taylor fallback is:
$atan(x+u) - atan(x) = |u\cdot \frac{1}{1+x^2}|$, 
so $l' = l - \log_2(1+ x^2)$.

## atan2

$\textrm{atan2}(y,x)$ represents the angle the vector of coordinates $(x,y)$ makes with the $(Ox^-)$ axis in the Cartesian plane. 

It is computed using the following formula:

$$\textrm{atan2}(y, x) =
\begin{align}
\arctan\left(\frac{y}{x}\right) &\text{ if } x > 0, \\
\arctan\left(\frac{y}{x}\right) + \pi &\text{ if } x < 0 \text{ and } y \ge 0, \\
\arctan\left(\frac{y}{x}\right) - \pi &\text{ if } x < 0 \text{ and } y < 0, \\
+\frac{\pi}{2} &\text{if } x = 0 \text{ and } y > 0, \\
-\frac{\pi}{2} &\text{if } x = 0 \text{ and } y < 0, \\
\text{undefined} &\text{ if } x = 0 \text{ and } y = 0.
\end{align}
$$

### Interval

To determine its output interval, we have to compute the maximum and minimum value it takes over the input intervals $[\underline{x};\overline{x}]\times [\underline{y};\overline{y}]$. These values depend on the quadrant(s) with which this domain intersects.

This function presents a discontinuity on the negative part of the x-axis, where the angle leaps from $\pi$ to $-\pi$.

If the domain over which we study $\textrm{atan2}$ intersects this discontinuity, we separate it between positive and negative values of $y$, so that $\textrm{atan2}$ is continuous over each of these sub-domains.

Then, when studying a domain over which $\textrm{atan2}$ is continuous, the maximum and the minimum are determined in function of the quadrant(s) intersected by the domain, considered in clockwise order for the maximum and counter-clockwise for the minimum.

* If the domain intersects $\{(x, y) | x \le 0, y\ge0\}$, then the maximum of $atan2$ over the domain is attained at $(\underline{x}, \underline{y})$.
* The next quadrant to consider is $\{(x, y) | x \ge 0, y\ge0\}$, where the maximum would lie at $(\underline{x},\overline{y})$.
* If there is an intersection with $\{(x, y) | x \ge 0, y< 0\}$ but not the previous quadrants, look for the maximum at $(\overline{x}, \overline{y})$.
* Finally, if the domain is entirely contained in the last quadrant $\{(x, y) | x > 0, y<0\}$, the maximum will be at $(\overline{x}, \underline{y})$.

The location of the minimum is determined with a similar method.

### Precision

Computing the exact minimum precision verifying pseudo-injectivity is difficult for the $\textrm{atan2}$ function, as it has two arguments and is non-linear. We use the formula above to compute a sound precision, by applying the pseudo-injectivity closed formula on function $\arctan$ and an input precision corresponding to the output precision of $\frac{y}{x}$.

We were not able to verify that this solution is optimal.

This is not the most satisfying solution, but since $\textrm{atan2}$ is not a function that tends to be used a lot in Faust programs, the impact should be limited.

## Atanh

atanh is the reciprocal of the tanh (hyperbolic tangent) function.

It is defined over the $]-1;1[$ open interval. It tends towards $\pm \infty$ at the bounds of the interval. To account for the exclusion of $-1$ and $1$ from the definition interval, the `std::nexttoward` function is used to generate bounds as close to the excluded values as the `double` format permits.

An empty input interval gives an empty output interval.

Its derivative attains a minimum at $0$: we compute the precision at the point $x$ of smallest magnitude in the interval ($0$ if it is included).

The Taylor fallback is computed as follows:
$|atanh(x+u) - atanh(x)| = |u \cdot \frac{1}{1 - x^2}|$,
so $l' = l - \log_2(1 - x^2)$.

## Button

A button is a Faust UI element that takes values $0$ when it is up and $1$ when it is down.
It thus outputs a boolean signal, that only needs $1$ bit to be represented, and by consequence a LSB of $l'=0$.

## Ceil

Ceil is the ceiling function $x \mapsto \lceil x \rceil$, returning the integer right above its argument. 

It returns an empty interval on an empty argument.

Although its result is mathematically an integer, the implementation of the primitive uses `std::ceil`, whose return value is a floating-point number.
We thus chose to set the output precision to $-1$, which might be under-optimal, but will force the interval library and the code generator to treat the result as floating-point.

## Checkbox

Checkbox functions similarly to Button intervals-wise. Its boolean output warrants a precision of $l'=0$.

## Cos

Cos is the cosine trigonometric function.

It is a periodic function of period $2\pi$. 
Its derivative periodically attains $0$ at points of the form $k\cdot\pi$, $k\in\mathbb{Z}$. 
These values are irrationnal numbers, which are not representable using floating-point or fixed-point numbers.
The smallest gap between two consecutive images is attained at one of the multiple of $\pi$ present in the interval.
While it is possible to bound how close to a multiple of $\pi$ a fixed-point number will be, 
it is too costly to compute for this application, and we will consider that the difference computed at the smallest multiple of $\pi$ (i.e. 0) gives sufficient precision.

## Cosh

Cosh is the hyperbolic cosine function. 

It is defined by $cosh(x) = \frac{e^x + e^{-x}}{2}$ over the whole $\mathbb{R}$ set.

The derivative of $cosh$ is $sinh$, which attains its lowest absolute value $0$ in $0$. If $0 \in [\underline{x}; \overline{x}]$, then the precision is computed at 0, otherwise it is computed at the point of lowest magnitude: $\overline{x}$ if $|\overline{x}| < |\underline{x}|$, $\underline{x}$ otherwise.

The Taylor fallback is:
* if $x = 0$, we use the second-order Taylor expansion of $cosh$: $cosh(u) = \frac{u^2}{2} + o(u^2)$, so $l' = 2*l -1$.
* if $x\neq 0$, $cosh(x + u) - cosh(x) = |u \cdot sinh(x)|$ so $l' = l + \log_2(|sinh(x)|)$.

## Delay

Delay is a Faust primitive used to temporally shift a signal. It takes two arguments: the signal that is delayed, and the signal indicating the amount of delay.  

The output signal is composed of the same samples as the input signal, so it should be represented with the same precision. 

## Div

Div is the division operator. 
It is expressed as the composition of multiplication and the inverse. 
Both the bounds and the precision of the output interval are computed from the result of the Mul function.

## Eq

Eq is the boolean equality test between two signals. 
It can return either $0$ or $1$, so no more than 1 bit is needed to represent the result: output precision is $0$.

## Exp 

Exp is the exponential function. 

It is a convex function that is defined on the whole $\mathbb{R}$ set. 
It is equal to its own derivative and the derivative has a minimum in $\underline{x}$.

Due to the very high dynamic range of the exponential function, it actually has the potential to push usual floating-point formats to their limits (ie generate overflows and underflows).
To manage that phenomenon, precision computation is implemented differently for the exponential than for other real functions.

**TODO** Give a concrete example as to how a straightforward implementation could fail.

The general formula for the computation of the output precision is $\log_2(e^{x+u} - e^x)$.
Due to the properties of the exponential function, this expression can be rewritten as $\log_2(e^x\cdot(e^u -1)) = x\cdot\log_2(e) + \log_2(e^u - 1)$.
We denote $\delta = e^u - 1$, $p_1 = x \cdot \log_2(e)$ and $p_2 = \log_2(\delta)$ to simplify proofs.

If $u$ is too small, then $e^u$ is computed to be $1$ and $\delta = 0$.
The threshold for such a $u$ is $2\cdot 10^{-16}$ (the difference between `1` and the next representable `double`), corresponding to an input LSB of $l = -52$. 

If $u$ has large enough magnitude, $\delta$ can be computed directly as $e^u-1$ and $p_2$ as $\log_2(\delta)$.

If $\delta$ is computed to be $0$, we replace it by the Taylor approximation of $e^u - 1 \approx u$, and then $p_2 = l$.

## FloatCast

FloatCast is the operator casting values into floats.

If we cast an `int` value in `float` without modifying the precision, it will remain above $0$ and the output value will be interpreted as `int` because of that, defeating the purpose of a `float` cast.
Having observed that, we reassign precisions greater than $0$ to $-1$ in order to force the `float` typing.

If the input interval is empty, then the constructed ouput interval is empty as well.

## FloatNum

FloatNum is the operator constructing a singleton interval from a floating-point number.
The creation of the interval, including its precision, is deferred to the `singleton` function.

## Floor

Floor is the floor function $x \mapsto \lfloor x \rfloor$, returning the integer right above its argument.

It returns an empty interval on an empty argument.

Although its result is mathematically an integer, the implementation of the primitive uses `std::floor`, whose return value is a floating-point number.
We thus chose to set the output precision to $-1$, which might be under-optimal, but will force the interval library and the code generator to treat the result as floating-point.

## Ge

Ge is the "greater than or equal" ($\ge$) boolean test.
It can return either $0$ or $1$, so no more than 1 bit is needed to represent the result: output precision is $0$.

## Gt 

Gt is the "greater than" ($>$) boolean test.
It can return either $0$ or $1$, so no more than 1 bit is needed to represent the result: output precision is $0$.

## HSlider

HSlider is the horizontal slider primitive, used to input numbers contained in a given range, with a given step between two values of the slider.

The output interval represents all the values that can be represented by the slider.
The input intervals are those representing the lower bound of the slider's range, its higher bound, its initial value and its step.

If one of the numerical parameters of the slider is an empty interval, then the slider cannot generate values and its ouput is empty as well. 

The bounds of the output interval are computed from the bounds of the parameters' intervals:
if the minimum value of the slider lies in $[\underline{x_{lo}}; \overline{x_{lo}}]$ and its maximum value in $[\underline{x_{hi}}; \overline{x_{hi}}]$, its output values lie in $[\underline{x_{lo}}; \overline{x_{hi}}]$.

Discrete values that can be taken by a slider with range $[x_{lo}; x_{hi}]$ and step $s$ are of the form $x_{lo} + k \cdot s \le x_{hi}$ with $k\in \mathbb{Z}$.
The values of the slider are aligned with the lower bound and regularly placed according to the step.
In order to be able to correctly represent these values, the output precision should be the finer between the precision of the lower bound interval, the precision of the step interval, and the precision implied by the step (given by $\lfloor \log_2(s) \rfloor$).

## IntCast

IntCast is the operator casting values into integers.

The result is an integer and thus has precision 0.

## IntNum

IntNum is the operator constructing a singleton interval from an integer.
As an integer interval, its precision is set to $0$.

## Inv

Inv is the inverse function $x \mapsto \frac{1}{x}$.
It is defined over the $\mathbb{R} - \{0\}$.

The point of lowest slope is the point of highest magnitude.

The Taylor fallback is 
$\frac{1}{x + u} - \frac{1}{x} = \frac{-u}{x^2}$,
so $l' = l - 2 \cdot \log_2(|x|)$.

**Remarks:** 
* We explored whether we could use the 31 MSB cap on fixed-point formats to cap the LSB of the inverse to -31.
The experiments were not concluding, and we found that it is unwise to use properties of fixed-point numbers to infer properties that will be used by all manners of number representations.
* Similarly, we wondered whether we could assume that no actual division by 0 will occur to give a lower boundary on the magnitude of numbers that will be inverted, and by consequence a higher boundary on the result.
The experiments were not concluding either.

## Label

Label is used to create labels for UI elements.
It is a string and doesn't represent a numerical value, thus its assigned interval is the empty interval.

## Le

Le is the "lower than or equal" ($\le$) boolean test.
It can return either $0$ or $1$, so no more than 1 bit is needed to represent the result: output precision is $0$.
It is implemented by inverting the arguments of `Ge`.

## Log

Log is the base $e$ logarithm, or natural logarithm, denoted $\ln$ in mathematics.
It is defined over the $\mathbb{R}_+^*$, and has a limit of $-\infty$ in $0$.

Its lowest slope is attained at the higer bound of the interval.

The Taylor fallback is expressed as:
$\log(x+u) - \log(x) = \frac{u}{x}$,
from which $l' = l - \log_2(|x|)$.

## Log10 

Log10 is the base $10$ logarithm. 

Like for the natural logarithm, its lowest slope is attained at the higher bound of the interval.

The Taylor fallback is expressed as:
$\log_{10}(x+u) - \log_{10}(x) = \frac{u}{x \cdot \ln(10)}$, 
thus $l' = l - \log_2(|x|) - \log_2(\ln(10))$.

## Lsh

Lsh is the binary left-shift operator, also denoted `<<` in various languages.

*Example:* `0b00001010 << 2 = 0b00101000`

In terms of precision, this will shift the LSB by the number of positions $k$ indicated by the bitshift: 
$l' = l + k$.

The change in MSB will be reflected by a change of mangnitude in the bounds of the intervals: 
output interval will be $[\underline{x}\cdot 2^k; \overline{x}\cdot 2^k]$.

## Lt  

Lt is the "lower than" ($<$) boolean test.
It can return either $0$ or $1$, so no more than 1 bit is needed to represent the result: output precision is $0$.
It is implemented by inverting the arguments to `Gt`.

## Max

Max is the function returning the maximum of two values.

The interval representing the result must be able to represent values from both arguments' intervals, and thus should be as precise as the most precise of both:
$l' = min(l_x, l_y)$.

## Mem

Mem is the one-sample delay operator.

The output signal contains the same value as the input signal, and thus should be represented with the same precision.

## Min

Min is the function returning the minimum of two values.

The interval representing the result must be able to represent values from both arguments' intervals, and thus should be as precise as the most precise of both:
$l' = min(l_x, l_y)$.

## Missing

The `intervalMissing.cpp` file implement placeholders for primitives that have not been dealt with yet.

## Mod

Mod is the modulo operator. 

It is implemented using C++ function `std::fmod`, which has the property that $\mathrm{fmod}(x, y) = x - y \cdot \lfloor \frac{x}{y} \rfloor$. 
The value $\lfloor \frac{x}{y}\rfloor$ will be called the quotient of $x$ by $y$.

Let's denote the input intervals $X = [\underline{x}; \overline{x}]$ for the dividend and $Y =[\underline{y}; \overline{y}]$ for the divisor.
The case where both intervals are positive is a particular case, implemented with function `positiveFMod`, and the general case can be reduced to the positive case by splitting the arguments intervals between positive and negative values, which is done in function `Mod`.

In the positive case, the result interval is automatically included in $[0; \overline{y}[$.
This is a direct consequence of the definition of $\lfloor~\rfloor$:
$\lfloor \frac{x}{y} \rfloor \le \frac{x}{y} < \lfloor \frac{x}{y} \rfloor + 1$,
so $y\cdot \lfloor \frac{x}{y} \rfloor \le x < y \cdot \lfloor \frac{x}{y} \rfloor + y$
and $0 \le x - y\cdot \lfloor \frac{x}{y} \rfloor  < y \le \overline{y}$.

The $\mathrm{fmod}$ function is a piecewise linear function, with discontinuities along the lines of equation $x = k\cdot y$, $k \in \mathbb{Z}$. 
These discontinuities correspond to the places where $x$ is a multiple of $y$ (and the remainder is thus null).
- On the $x \ge ky$ side of the line, the limit of $\mathrm{fmod}(x, y)$ is $0$ as we approach the line and the quotient $\lfloor\frac{x}{y}\rfloor = k$.
- On the $x < ky$ side, the limit of $\mathrm{fmod}(x, y)$ is $y$ and $\lfloor \frac{x}{y} \rfloor = k-1$.

We have to distinguish two cases, depending on whether the interval product $X \times Y$ spans a discontinuity or not.
Let us pose $n = \lfloor \frac{\underline{x}}{\overline{y}}\rfloor$, which is the smallest value that can be taken by the integral quotient between $x$ and $y$.


If it doesn't span a discontinuity, the quotient $\lfloor \frac{x}{y}\rfloor$ is constant with value $n$.
The $\mathrm{fmod}$ function is thus linear over the domain: $\mathrm{fmod}(x, y) = x - n\cdot y$, and has maximum value $\overline{x} - n \cdot\underline{y}$ and minimum $\underline{x} - n \cdot \overline{y}$.
This is the case represented by the last `return` statement in the implementation.

If the domain does span discontinuities, the minimum value is $0$ and is attained on the discontinuities.
Regarding the maximum, the function $\mathrm{fmod}$ presents local suprema along these discontinuities of equation $x = k\cdot y, ~k \ge n+1$.
We exclude the discontinuity $x = n \cdot y$ because $n = \lfloor \frac{\underline{x}}{\overline{y}} \rfloor$ (by definition of $n$) means that $\underline{x} \ge n \overline{y}$ (by definition of $\lfloor ~\rfloor$), from which we can deduce that $\forall x \in X, y \in Y, x > n \cdot y$, except maybe in the case $x = \underline{x}, y = \overline{y}$, where the discontinuity line intersects the rectangle domain in the single point $(\underline{x}, \overline{y})$.

We are looking for the largest value of $\mathrm{fmod}(x, y)$. 
We have seen that the limit of this value is $y$ along lines of equation $x = k\cdot y, ~k\ge n+1$.
We are thus looking for the largest $y$ that verify $\exists x \in X, k \ge n+1, y = \frac{1}{k} x$.

The value of the global upper bound we are looking for is thus in $Y \cap \bigcup_{k \ge n+1} \frac{1}{k} X = \bigcup_{k\ge n+1} Y \cap \frac{1}{k}X$, and is the maximum of this set.
Since all $Y \cap \frac{1}{k}X$ are (possibly empty) intervals, $max(\bigcup_{k\ge n+1} Y \cap \frac{1}{k}X) = max_{k\ge n+1} (max (Y \cap \frac{1}{k} X))$.

For $k \ge n+1$, $max(Y\cap \frac{1}{k} X) = min(\overline{y}, \frac{1}{k}\overline{x}) \le min(\overline{y}, \frac{1}{n+1}\overline{x}) = max(Y\cap \frac{1}{n+1} X)$, so $max_{k\ge n+1} (max (Y \cap \frac{1}{k} X)) = min(\overline{y}, \frac{1}{n+1}\overline{x})$.

Thus, the upper bound of the output interval is:
- $\overline{y}$ if $\overline{y} \ge \frac{1}{n+1} \overline{x}$
- $\frac{1}{n+1} \overline{x}$ if $\overline{y} < \frac{1}{n+1} \overline{x}$ (and $\underline{y} \le \frac{1}{n+1}\overline{x}$ because the domain $X \times Y$ has an intersection with the $x = (n+1)\cdot y$ discontinuity).

### Precision

Let us assume that $x$ is represented with precision $l_x$, from which the ulp $u_x$, and $y$ has precision $l_y$ and ulp $u_y$.
If we assume that $u_x < u_y$, we can write $x = p \cdot u_x$ and $y = q \cdot u_x$, where $p, q \in \mathbb{Z}$.
Then, $\mathrm{fmod}(x, y) = x - \lfloor\frac{x}{y}\rfloor \cdot y = (p - \lfloor \frac{x}{y} \rfloor \cdot q)\cdot u_x$, thus $\mathrm{fmod}(x, y)$ can be represented with precision $l_x$.
A similar reasoning takes place when $u_y < u_x$.

Output precision is the finest precision between those of the two arguments: this will be the one to determine which discrete grid the result aligns with.

### Remark on the sampling tests

Some inputs can measure an output interval that excludes $0$ as a lower bound.
This is normal and due to the fact that the discontinuities $x = k\cdot y$ do not always intersect the discrete fixed-point grid from which `analyzeBinaryMethod` samples.
Validity of the result on these inputs should instead be checked using `check` (which implies computing the output interval by hand).

## Mul 

$0$ is considered an absorbant element when it comes to interval bounds.
Indeed, even though $0\times\infty$ is an indeterminate form in the realms of mathematics, a bound of $\infty$ does not mean that values in the interval can be infinite, but rather that they are unbounded.
They nonetheless remain finite, and hence their multiplication by $0$ yields $0$.

### Overflow

When both input intervals represent integers, there is a possibility for overflow to occur.
These happen when multiplying two values in the intervals $x$ and $y$ yields a result $|x\times y| \ge$ `INT_MAX, |INT_MIN|`.
In practice, the test detecting overflows is $\max(|\underline{x}|, |\overline{x}|) \times \max(|\underline{y}|, |\overline{y}|) \ge$ `INT_MAX`.

If an overflow does occur, then the output value can pass from `INT_MIN` to `INT_MAX` (or conversely), thus the output interval will be between those two values (cast to `double`).
This interval is not optimal, but it seems difficult to compute the actual interval.

If an overflow does not occur, the boundaries are computed by casting the boundaries to `int`, computing all the possible multiplications combinations, and setting the minimum and the maximum as boundaries for the output interval.

### Precision

The output precision is the sum of the input precisions.

**TODO** Proof that it is optimal

## Ne

Boolean operator testing whether two arguments are not equal.

If the intervals do not intersect (the lower bound of one is higher than the higher bound of the other), the result is necessarily $1$ (`True`).
If the intervals are singletons and are equal, then the result is necessarily $0$ (`False`).
Otherwise, both values are possible.

Being a boolean function, its output has precision $0$.

## Neg

Computes the opposite of a value: $x \mapsto -x$.

The precision remains the same.

## Not 

Bitwise negation on integers.

The output interval is computed by brute-force iterating through the integer elements of the input interval.

The output is made of integers, so its precision is at least $0$.
If the input precision was above $0$, it remains identical.

## NumEntry 

UI primitive used to input a number.

Like for `HSlider`, the precision is the finest of those induced by the lower bound and the step.

## Or

Or is a bitwise boolean "or" function.

Example: `Or((0101)b, (0110)b) = 0111b`.

The reasoning is similar to that for `And`, except that the trailing *zeroes* are here trailing *ones*.

The final implementation prioritises keeping all bits, as they are deemed to all be significant; however, an implementation prioritising pseudo-injectivity is present in comments.

**TODO:** Document `bitwiseSignedOr`.

## Pow 

Pow is the binary "power" function: $(x, y) \mapsto x^y$.

It implementation is split into two functions: `fPow`, that deals with elevating a real base value to a power that is a real number as well, and `iPow`, that deals with elevating a real value to an integer power.

This distinction is necessary because while it is possible to elevate a negative base value to an integer power, it is a domain violation for a real-valued power.

### fPow

`fPow` is expressed as the composition of `Exp`, `Mul` and `Log`, following the formula $x^y = e^{y\cdot \ln(x)}$.

The interval of the first argument $x$ is constrained to positive values.

### iPow

`iPow` is expressed with the help of the auxiliary function `ipow`, which computes the result of an real-valued interval elevated to a single integer power (instead of an interval of integers).

A case disjunction is operated in `ipow`, depending on whether the exponent `k` is even () or odd ().

The precision is computed using the same principle as `exactPrecisionUnary`, but without actually using it, since the function we are computing the precision of is parameterised by `k`.
The finest precision is attained at the lower bound $\underline{x}$ of the interval: we are computing $\log_2(|(\underline{x}+u)^k - \underline{x}^k|)$.
- If $\underline{x} = 0$, this formula is reduced to $\log_2(u^k) = k\cdot l$.
- Otherwise, $(\underline{x} + u)k - \underline{x}^k \approx u\cdot k \cdot x^{k-1}$, thus $\ell_{out} = \log_2(k) + \ell + (k-1)\cdot \log_2(\underline{x})$, which we split into $p_1 = (k-1)\cdot \log_2(\underline{x})$ and $p_2 = $.

## Rem

Rem is the function giving the remainder of the floating-point division. 
It is defined as $\mathrm{rem}(x, y) = x - \lfloor \frac{x}{y} \rceil \cdot y$, where $\lfloor \frac{x}{y} \rceil$ is the closest integer to $\frac{x}{y}$ (chosed to be an even number if $\frac{x}{y}$ is exactly halfway between two integers).

The reasoning is very similar as the one for Mod.

The discontinuities are lines of equation

## Rint

Rint is the function rounding its argument to an integer value, using the current rounding mode.
It is implemented using the `std::rint` C++ function.

The boundaries of the output interval are the result of rounding the input boundaries, and since the output is an integer, precision is set to $0$ (unless the input precision was already $\ge 0$, in which case it remains unchanged).

**Remark:** The rounding mode might be different when compiling the Faust compiler and when compiling the produced code. 
This might cause incoherence.

## Round 

Round is the function rounding its argument to the *nearest* integer value, regardless of the current rounding mode.
It is implemented using the `std::round` C++ function.

The boundaries of the output interval are the result of rounding the input boundaries to nearest, and since the output is an integer, precision is set to $0$ (unless the input precision was already $\ge 0$, in which case it remains unchanged).

## Rsh 

Rsh is the binary right-shift operator, also denoted `>>` in various languages.

*Example:* `0b00001010 >> 2 = 0b00000010`

In terms of precision, this will shift the LSB by the number of positions $k$ indicated by the bitshift: 
$l' = l - k$.

The change in MSB will be reflected by a change of mangnitude in the bounds of the intervals: 
output interval will be $[\underline{x}/2^k; \overline{x}/2^k]$.

## Sin

Sin is the cosine trigonometric function.

It is a periodic function of period $2\pi$. 
Its derivative periodically attains $0$ at points of the form $k\cdot 2\pi + \frac{\pi}{2}$, $k\in\mathbb{Z}$. 
These values are irrationnal numbers, which are not representable using floating-point or fixed-point numbers.
The smallest gap between two consecutive images is attained at one of these values present in the interval.
While it is possible to bound how close to such a value a fixed-point number will be, 
it is too costly to compute for this application, and we will consider that the difference computed at the smallest value of the form $k\cdot 2\pi + \frac{\pi}{2}$ (i.e. $\frac{pi}{2}$) gives sufficient precision.

Taylor fallback:

**Remark:** This might cause some discrepancy between the precision computed by the test function and that computed by the interval version of `Sin`. 
A LSB difference $\gtrsim -10$ should not be taken to mean that the function implementation is false.

## Sinh

Sinh is the hyperbolic sine function.
It is defined as $sinh(x) = \frac{e^x - e^{-x}}{2}$ over the whole $\mathbb{R}$ set.

The derivative of $sinh$ is $cosh$, which attains its lower value at $0$.
Like for cosh, the precision is computed at the point $x$ in the interval of lowest magnitude.

The Taylor fallback is:
- if $x = 0$, we use the second-order Taylor approximation of $sinh(u) \approx $
- if $x\neq 0$, $\ell' = \ell + \log_2(cosh(x))$ 

## Sqrt 

Sqrt is the square root function. 
It is defined over the $\mathbb{R}_{+}$ set.

Its derivative is decreasing and thus attains its lower value at the highest bound of the interval.

The Taylor fallback is:
- if $\overline{x} = 0$ (the definition interval is a singleton), then the second-order Taylor development is used: $\sqrt{u} \approx$
- otherwise, $\sqrt{x+u} - \sqrt{x}$

## Sub

Sub is the subtraction function.

It is equivalent to the composition of the addition and the opposite function, and as such, its output precision is the finest of the input precisions.

## Tan

Tan is the trigonometric tangent function. 

It is periodic of period $\pi$.
It periodically present vertical asymptotes at points $k\cdot \pi + \frac{\pi}{2}$, $k \in \mathbb{Z}$, the asymptote tending to $+\infty$ before these points and to $- \infty$ after.
Its derivative periodically attains $0$ at points of the form $k\cdot \pi$, $k \in \mathbb{Z}$.

Similarly to what has been done for sin and cos, we do not compute the tightest precision respecting pseudo-injectivity, but an approximation computing the precision at multiples of $\pi$.

Taylor fallback:

**Remark:** This might cause some discrepancy between the precision computed by the test function and that computed by the interval version of `Sin`. 
A LSB difference $\gtrsim -10$ should not be taken to mean that the function implementation is false.

## Tanh 

Tanh is the hyperbolic tangent function.
It is defined as $tanh(x) = \frac{sinh(x)}{cosh(x)}$.

Its derivative tends to $0$ in $\pm \infty$, thus the LSB is computed at the interval boundary with the highest absolute value.

Taylor fallback:


## VSlider

VSlider is the horizontal slider primitive, used to input numbers contained in a given range, with a given step between two values of the slider.

The output interval represents all the value that can be represented by the slider.
The input intervals are those representing the lower bound of the slider's range, its higher bound, its initial value and its step.

The bounds of the output interval are computed from the bounds of the parameters' intervals:
if the minimum value of the slider lies in $[\underline{x_{lo}}; \overline{x_{lo}}]$ and its maximum value in $[\underline{x_{hi}}; \overline{x_{hi}}]$, its output values lie in $[\underline{x_{lo}}; \overline{x_{hi}}]$.

Discrete values that can be taken by a slider with range $[x_{lo}; x_{hi}]$ and step $s$ are of the form $x_{lo} + k \cdot s \le x_{hi}$ with $k\in \mathbb{Z}$.
The values of the slider are aligned with the lower bound and regularly placed according to the step.
In order to be able to correctly represent these values, the output precision should be the finer between the precision of the lower bound interval, the precision of the step interval, and the precision implied by the step (given by $\lfloor \log_2(s) \rfloor$)

## Xor

Xor is the bitwise "exclusive or" function.

Similarily to And and Or, the output precision preserving all the bits is the finest of both input precisions.

If we chose to prioritise pseudo-injectivity, 

**TODO:** Document `bitwiseSignedXOr`.



# Avoiding absorption

In the case of $pow$, we chain the computations of LSBs associated to $ln$, multiplication and $exp$.
These can quickly get very fine precisions, and result in computations of the form $f(x + u) - f(x)$ where $u$ is so small compared to $x$ that in machine representation $x + u = x$, resulting in an infinite computed precision.

We thus need to rewrite precision computations so that these forms do not appear.
We need to rewrite them in the case of $log$ and $exp$.



## Integer power

When both input intervals represent integers, Faust uses 32-bit integer
arithmetic with two's-complement wrapping semantics. Generated C/C++ code must
be compiled with `-fwrapv` when using GCC or Clang so signed addition,
subtraction and multiplication have the same semantics as the interval model.
This option has to be applied when compiling the generated DSP code, not only
when compiling the Faust compiler itself.

If any value in the input intervals can overflow while computing the power,
the interval implementation currently returns `[INT32_MIN; INT32_MAX]`. This
is a conservative envelope: wrapped powers can form a disconnected set that a
single interval cannot in general represent precisely.

Let's study $x \mapsto x^n$ for $n \in \mathbb{N}$.

The absolute lowest value of the slope is attained at zero, so if the input interval contains zero, the precision is $\lfloor log_2((2^l)^n - 0)\rfloor = l\cdot n$.

If zero is not in the interval, the precision is attained at $\tilde x$, the boundary of the interval at the lowest absolute value. 

Let's pose $\delta = (\tilde x+2^l)^n - \tilde x^n = \tilde x^n\cdot\left((1 + \frac{2^l}{\tilde x})^n - 1\right)$, with the precision being $\lfloor log\_2(\delta)\rfloor$. 
If $2^l$ is very small, $(1+ \frac{2^l}{\tilde x})^n \approx 1 + n\cdot \frac{2^l}{\tilde x}$, and thus the precision is $\lfloor log\_2(\tilde x^n \cdot n \cdot \frac{2^l}{\tilde x})\rfloor = \lfloor (n-1)\cdot log\_2(\tilde x) + log\_2(n) + l\rfloor$.

The logarithmic terms must be added before applying `floor`. Rounding each
term separately can produce an output LSB that is too coarse, notably when
$|\tilde x| < 1$.

For an interval of integer exponents $[n_0;n_1]$, results obtained with two
different exponents must also be distinguishable. Since a fixed-point input is
of the form $m\,2^l$, the powers lie on a common grid. We use $2^{l n_1}$ when
$l < 0$ and $2^{l n_0}$ when $l \ge 0$ as conservative precision bounds across
exponents.

## 12/05 update

In a shocking turn of events, it appears that computing Pow is not the only context in which functions are composed. 
In fact, this happens about everywhere in any slightly context program.
Thus, any function can potentially get an argument with a precision arising from a previous argument computation, which can get low enough that $x+u$ and $x$ have the same representation. 
In consequence, we apply the same principle as above to every primitive real function.

### Arccos

If $u << x$, $$arccos(x+u) - arccos(x) = u \cdot arccos'(x) + o(u) = \frac{u}{\sqrt{1 - x^2}} + o(u)$$

If the input interval contains $0$, at which the slope of the function is the lowest, 
then $l\_{out} = \lfloor log\_2(u)\rfloor = l$.

If it does not, then it is computed at point $x$ plosest to $0$ and the precision is 
$l\_{out} = \lfloor log\_2(\frac{u}{\sqrt{1 - x^2}}) \rfloor = \lfloor l\_{in}- \frac{1}{2} \cdot log\_2(1 - x^2) \rfloor$.

Note that one could prove that if $|x| < \frac{\sqrt{3}}{2}$, then $l\_{out} = l\_{in}$ 
(but since this is not the case for all x we won't use this property in the implementation).

### SinPi

If $u << x$, 
$$sinPi(x + u) - sinPi(x) = u \cdot \pi \cdot sinPi'(x) + \frac{u^2}{2}\cdot \pi^2 \cdot sinPi''(x) + o(u^2) = u \cdot \pi \cdot cosPi(x) - \frac{u^2}{2}\cdot \pi^2 \cdot sinPi(x) + o(u^2)$$.

If the input interval contains a half-integer (for instance $\frac{1}{2}$), 
then precision is computed at that point and becomes $l\_{out} = \lfloor log\_2(\frac{u^2}{2}\cdot \pi^2) \rfloor = 2\cdot(l_in + \lfloor log\_2(\pi)\rfloor) - 1 = 2\cdot l - 1$.

If it doesn't, then precision is computed at point $x$ closest to a half-integer and becomes $l\_{out} = \lfloor log\_2(u\cdot\pi\cdot cosPi(x)) \rfloor = l + \lfloor log\_2(\pi) + log\_2(cosPi(x))\rfloor$

# Backwards precision computation

i.e., given an output precision $l_{out}$, determine an input precision $l_{in}$ such that $l_{in}$ would induce output precision $l_{out}$ in the direct direction.

We suppose that we know at which point $x_0$ to compute the precision in the direct direction if given an input precision.
So, we know that the relationship between $l_{in}$ and $l_{out}$ is as follows:
$l_{out} = ⌊log_2(f(x_0+2^{l_{in}}) - f(x_0))⌋$
(we do away with the absolute value and the $±$ without loss of generality)
This can be translated as a pair of inequalities:
$l_{out} ≤ log_2(f(x+2^{l_{in}}) - f(x)) < l_{out} + 1$
By inverting these inequalities and inverting $f$ (assuming that it has the right properties for that) we get:
$log_2(f^{-1}(2^{l_{out}} + f(x_0)) - x_0) ≤ l_{in} < log_2(f^{-1}(2^{l_{out}+1} + f(x_0)) - x_0)$

The first inequality is the most important: in the direct direction, it ensures that the choice of precision is sound, ie that it won't cause any collision.
In the backwards direction, it means that the input precision shouldn't be chosen lower than this limit, lest we end up with inputs that the chosen output precision won't be able to discriminate.

The other inequality is less crucial: it is only here to ensure that the precision chosen is tight, ie that we don't choose an output precision that is too precise wrt what we need (or, in the reverse sense, that we are not too sloppy with the input precision).

So, in conclusion, the formula of $l_{in}$ in function of $l_{out}$ is:
$l_{in} = ⌈log_2(f^{-1}(2^{l_{out}} + f(x_0)) - x_0)⌉$.
By design, it always fulfills the important inequality, and it may or may not fulfill the other, but we shouldn't be too concerned about this.
