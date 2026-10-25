// Each check computes an interval and a value the compiled program produces : the value must
// be in the interval (the soundness criterion Y within Y^ of the README). One check per
// witness of a defect found by a review of the library.

#include <climits>
#include <cmath>
#include <cstdint>

#include "interval/affine_ops.hh"
#include "interval/affint.hh"
#include "interval/check.hh"
#include "interval/double_double.hh"
#include "interval/interval_algebra.hh"
#include "interval/interval_def.hh"

using namespace itv;

namespace {

bool contains(const interval& r, double v)
{
    return r.lo() <= v && v <= r.hi();
}

std::int32_t wrapAdd(std::int32_t a, std::int32_t b)
{
    return std::int32_t(std::uint32_t(a) + std::uint32_t(b));
}
std::int32_t wrapMul(std::int32_t a, std::int32_t b)
{
    return std::int32_t(std::uint32_t(a) * std::uint32_t(b));
}

}  // namespace

int main()
{
    programPrecision() = 2;
    const interval_algebra A;
    const double           T = 1000;
    const affine_algebra   F(T);
    const AffItv           counter{0, 1, 0, 1, 0};  // an int counter : exactly t

    // a divisor beyond int32 has wrapped : 1.0 / 1 is reached
    check("witness : 1.0 / [1, +inf] integer holds 1", true,
          contains(A.Div(A.FloatNum(1.0), interval(1, HUGE_VAL, 0)), 1.0));
    // exp10 of an int computes a floating power
    check("witness : affine exp10(10) holds 1e10", true,
          contains(toItv(F.Exp10(F.IntNum(10)), T), std::pow(10.0, 10.0)));
    // highest of a float keeps the nature of the float
    {
        const interval x = A.Mul(interval(0, 1, -24), A.FloatNum(1e10));
        check("witness : highest(float [0, 1e10]) + 0.5 holds 1e10 + 0.5", true,
              contains(A.Add(A.Highest(x), A.FloatNum(0.5)), 1e10 + 0.5));
    }
    // an int product that wraps over [0, T], then delayed
    {
        const AffItv q = F.Delay(F.Mul(counter, F.IntNum(3000000)), F.IntNum(284));
        check("witness : affine (t * 3000000)@284 holds its value at t = 797", true,
              contains(toItv(q, T), double(wrapMul(797 - 284, 3000000))));
    }
    // the int conversion of an int that has wrapped is the identity
    {
        const AffItv c = F.Add(F.IntNum(-2000000000), F.IntNum(-1000000000));
        const AffItv k = F.IntCast(F.Add(F.Mem(c), AffItv{0, 0, 0, 1, 0}));
        check("witness : affine int(c@1 + t) holds its wrapped value", true,
              contains(toItv(k, T), double(wrapAdd(wrapAdd(-2000000000, -1000000000), 1))));
    }
    // 1 / ceil(-0.3) = 1 / -0 = -inf
    check("witness : 1 / ceil([-0.5, 3]) holds -inf", true,
          contains(A.Inv(A.Ceil(interval(-0.5, 3, -24))), 1.0 / std::ceil(-0.3)));
    // the phase of a counter that has wrapped is negative
    {
        const std::int32_t xw = std::int32_t(std::uint32_t(3000000000u));
        check("witness : phase([0, +inf] integer, 5) holds a wrapped phase", true,
              contains(A.Phase(interval(0, HUGE_VAL, 0), interval(5, 5, 0)), double(xw % 5) / 5.0));
    }
    // a select2 joined with a fixed corridor keeps the motion of the counter past T
    {
        const AffItv s = F.Select2(F.Checkbox(F.Label("s")), counter, fromItv(interval(0, 5000, 0)));
        check("witness : select2(s, t, [0, 5000]) + 1 may wrap past T", true, F.Add(s, F.IntNum(1)).mayOverflow);
    }
    // the zero of pow(-inf, y < 0) is a float
    {
        const interval p = A.Pow(interval(-HUGE_VAL, -2, -24), interval(-0.7, -0.2, -24));
        check("witness : pow(-inf, -0.5) + INT_MAX + 1 holds 2^31", true,
              contains(A.Add(A.Add(p, A.IntNum(2147483647)), A.IntNum(1)),
                       std::pow(-HUGE_VAL, -0.5) + 2147483647.0 + 1.0));
    }
    check("witness : 2.5 is not a power of 2", false, interval(2.5).ispowerof2());
    // frac of an unbounded product contracted with its subtraction : fma(a, b, -floor(a * b))
    {
        volatile double a = std::ldexp(1.0, 30) - 1, b = std::ldexp(1.0, 30) + 1;
        check("witness : frac([0, +inf]) holds a contracted remainder", true,
              contains(A.Frac(interval(0, HUGE_VAL, -24)), std::fma(a, b, -std::floor(a * b))));
    }
    // in double, the program rounds (t + 0.1) - t at every sample
    {
        const AffItv x{0, 1, 0, 1, -24};
        const AffItv v = F.Sub(F.Add(x, F.FloatNum(0.1)), x);
        volatile double t1 = 1, u1 = t1 + 0.1, v1 = u1 - t1;
        check("witness : affine (t + 0.1) - t holds its rounded value at t = 1", true, v.lo(1) <= v1 && v1 <= v.hi(1));
        check("witness : affine ((t + 0.1) - t) == 0.1 is not folded to true", true,
              contains(toItv(F.Eq(v, F.FloatNum(0.1)), T), v1 == 0.1 ? 1.0 : 0.0));
    }
    // a double accumulator y = y' + 0.1 : its widened corridor is stationary and holds the
    // rounded sums of the program
    {
        const AffItv inc = F.FloatNum(0.1);
        AffItv       y   = aempty();
        bool         done = false;
        for (int it = 1; !done && it <= 100; it++) {
            AffItv f = F.Add(F.Delay(y, F.IntNum(1)), inc);
            if (it > 8) f = awiden(y, f, T);
            done = !y.isEmpty() && aleq(f, y, T);
            y    = f;
        }
        volatile double acc = 0;
        bool            in  = true;
        for (int t = 1; t <= T; t++) {
            acc = acc + 0.1;
            in  = in && y.lo(t) <= acc && acc <= y.hi(t);
        }
        check("witness : affine double accumulator is stationary", true, done && std::isfinite(y.b0) && y.b1 > 0);
        check("witness : affine double accumulator holds the rounded sums", true, in);
    }
    // in quad precision, exp of an unbounded value reaches +inf and stays a real
    {
        const int saved    = programPrecision();
        programPrecision() = 3;
        const interval e   = A.Exp(interval(-HUGE_VAL, HUGE_VAL, -24));
        programPrecision() = saved;
        check("witness : quad exp of an unbounded value holds +inf", true, contains(e, HUGE_VAL));
        check("witness : quad exp of an unbounded value stays a real", true, e.lsb() < 0);
    }
    // round, as the generated code : halves away from zero
    check("witness : affine round(2.5) holds 3", true, contains(toItv(F.Round(fromItv(interval(2.5, 2.5, -1))), T), 3.0));
    // prefix(1000, t) - t is -1 from t = 1 on
    {
        const AffItv d = F.Sub(F.Prefix(F.IntNum(1000), counter), counter);
        check("witness : affine prefix(1000, t) - t holds -1 at t = 1", true, d.lo(1) <= -1 && -1 <= d.hi(1));
        const AffItv w = F.WRTbl(F.IntNum(16), F.IntNum(1000), fromItv(interval(0, 1, 0)), counter);
        const AffItv e = F.Sub(F.RDTbl(w, F.IntNum(0)), counter);
        check("witness : a table cell written at 0 read at 500, minus t, holds -500", true, contains(toItv(e, T), -500.0));
    }
    // --- second review ---
    // the join is idempotent : a double accumulator driven through FixPointUpdate is stationary
    {
        AffItv y    = aempty();
        bool   done = false;
        for (int it = 1; !done && it <= 200; it++) {
            AffItv f = F.FixPointUpdate(y, F.Add(F.Mem(y), F.FloatNum(0.1)));
            if (it > 8) f = awiden(y, f, T);
            done = !y.isEmpty() && aleq(f, y, T);
            y    = f;
        }
        const AffItv j = ajoin(y, y, T);
        check("witness : a joined double accumulator is stationary", true, done);
        check("witness : ajoin(y, y) is y", true, j.a0 == y.a0 && j.a1 == y.a1 && j.b0 == y.b0 && j.b1 == y.b1);
        // a join that returns the larger corridor keeps the finer precision : a float zero
        // joined to an int-stamped seed is a float
        const AffItv seed = fromItv(interval(-HUGE_VAL, 16777216.0, 0));
        check("witness : a join within a corridor keeps the finer precision", true,
              ajoin(seed, F.FloatNum(0.0), T).lsb < 0);
    }
    // a widening holds fresh : the int recursion y = y' + y' + 1000 wraps, and its fixpoint holds
    // the wrapped value at t = 21
    {
        const double         T64 = 64;
        const affine_algebra F64(T64);
        AffItv               y    = aempty();
        bool                 done = false;
        for (int it = 1; !done && it <= 100; it++) {
            AffItv f = F64.Add(F64.Add(F64.Mem(y), F64.Mem(y)), F64.IntNum(1000));
            if (it > 8) f = awiden(y, f, T64);
            done = !y.isEmpty() && aleq(f, y, T64);
            y    = f;
        }
        std::int32_t v = 1000;
        for (int t = 1; t <= 21; t++) v = wrapAdd(wrapAdd(v, v), 1000);
        check("witness : the fixpoint of y = 2 y' + 1000 (int) holds its wrapped value", true,
              done && contains(toItv(y, T64), double(v)));
    }
    // a widening never narrows the side that does not widen : y = -(y' + [89041, 91976]) settles
    {
        const double         T20 = 20;
        const affine_algebra F20(T20);
        const AffItv         c    = fromItv(interval(89041, 91976, 0));
        AffItv               y    = aempty();
        bool                 done = false;
        for (int it = 1; !done && it <= 200; it++) {
            AffItv f = F20.Neg(F20.Add(F20.Mem(y), c));
            if (it > 8) f = awiden(y, f, T20);
            done = !y.isEmpty() && aleq(f, y, T20);
            y    = f;
        }
        check("witness : y = -(y' + [89041, 91976]) (int) settles", true, done);
    }
    // the meet set by a moving bound moves : [0, 100] meet min(t, 100) holds 11 at t = 11
    {
        const double         T10 = 10;
        const affine_algebra F10(T10);
        const AffItv         a = fromItv(interval(0, 100, -24));
        const AffItv         b = F10.Min(AffItv{0, 0, 0, 1, -24}, F10.FloatNum(100));
        const AffItv         r = meetDrift(fromItv(intersection(toItv(a, T10), toItv(b, T10))), a, b, T10);
        check("witness : a meet set by a rising bound holds its value past T", true, contains(forever(r, T10), 11.0));
    }
    // a leaky integrator y = 0.999 y' + 1 keeps a finite corridor
    {
        AffItv y    = aempty();
        bool   done = false;
        for (int it = 1; !done && it <= 100; it++) {
            AffItv f = F.Add(F.Mul(F.Delay(y, F.IntNum(1)), F.FloatNum(0.999)), F.FloatNum(1.0));
            if (it > 8) f = awiden(y, f, T);
            done = !y.isEmpty() && aleq(f, y, T);
            y    = f;
        }
        volatile double acc = 0;
        bool            in  = true;
        for (int t = 1; t <= T; t++) {
            acc = acc * 0.999 + 1.0;
            in  = in && y.lo(t) <= acc && acc <= y.hi(t);
        }
        check("witness : a leaky integrator keeps a finite corridor", true, done && std::isfinite(y.b0));
        check("witness : a leaky integrator corridor holds the program", true, in);
    }
    // a filter of null float coefficients is a float zero : fir(0.5, {0.0}) - INT_MIN is 2^31
    {
        const affine_algebra F10(10);
        const AffItv         r = F10.Sub(F10.Fir(F10.FloatNum(0.5), {F10.FloatNum(0.0)}), F10.IntNum(INT32_MIN));
        check("witness : fir(0.5, {0.0}) - INT_MIN holds 2^31", true, contains(toItv(r, 10), 2147483648.0));
        check("witness : interval fir of a null float coefficient is a float", true,
              A.Fir(A.FloatNum(0.5), {A.FloatNum(0.0)}).lsb() < 0);
    }
    // a delayed bound collapsed to a constant still moves past the horizon
    {
        const AffItv cross{-1610612736.0, 1, -1610612736.0, 1, -24};
        check("witness : mem of a rising lower bound moves past T", true, contains(forever(F.Mem(cross), T), -2e9));
    }
    // a table cell holds its generator at its index, past the horizon too
    {
        const affine_algebra F10(10);
        const AffItv         w = F10.WRTbl(F10.IntNum(100), F10.Gen(AffItv{0, 0, 0, 1, -24}), F10.IntNum(0), F10.FloatNum(0));
        check("witness : rdtable(100, g = t, 99) holds 99", true, contains(toItv(F10.RDTbl(w, F10.IntNum(99)), 10), 99.0));
    }
    // double-double at its edges : the rounding to float reads lo, an overflow is an infinity,
    // an intermediate overflow is computed again, a zero keeps its sign, a subnormal product
    // is the correctly rounded plain product
    {
        const double DMAX = 0x1.fffffffffffffp+1023;
        check("witness : ddToFloat(16777217 + 1e-10) is 16777218", true,
              ddToFloat(ddAdd({16777217, 0}, {1e-10, 0})) == 16777218.0f);
        check("witness : ddToFloat(16777217 - 1e-10) is 16777216", true,
              ddToFloat(ddAdd({16777217, 0}, {-1e-10, 0})) == 16777216.0f);
        check("witness : ddMul(DBL_MAX, 1 + 2^-53) overflows", true, std::isinf(ddMul({DMAX, 0}, {1, 0x1p-53}).hi));
        check("witness : ddAdd(DBL_MAX, 2^970 - 2^960) is DBL_MAX", true, ddAdd({DMAX, 0}, {0x1p970, -0x1p960}).hi == DMAX);
        check("witness : ddAdd(-0, -0) is -0", true, std::signbit(ddAdd({-0.0, 0}, {-0.0, 0}).hi));
        check("witness : ddMul(-1, 0) is -0", true, std::signbit(ddMul({-1, 0}, {0, 0}).hi));
        check("witness : ddDiv(-0, 1) is -0", true, std::signbit(ddDiv({-0.0, 0}, {1, 0}).hi));
        check("witness : a subnormal ddMul is the rounded product", true,
              ddMul({0x1.352c06d3f408bp-512, 0}, {0x1.18938fc7c6d0ap-508, 0}).hi == 0x1.52da5d649bc45p-1020);
    }
    // the libm within 2 ulps of the exact result, the ulp of its binade : just above a power
    // of two, 4 floats below it ; the image widened by 2 ulps (values at 1.9 to 1.98 ulps,
    // checked at 300 bits)
    {
        const int saved = programPrecision();
        programPrecision() = 2;
        check("witness : double log at 1.955 ulps below 16", true,
              contains(A.Log(A.FloatNum(0x1.0f2ebd0a80028p+23)), 0x1.ffffffffffffdp+3));
        check("witness : double asinh at 1.960 ulps below 32", true,
              contains(A.Asinh(A.FloatNum(0x1.1f43fcc4b663dp+45)), 0x1.ffffffffffffdp+4));
        // cosh(tiny) is just above 1 : 1 - 3 * 2^-53 is 1.5 ulps below it (sin, whose exact
        // result stays below 1, stays within 1 : 1 + 2^-52 is more than 2 ulps of the binade
        // below away from it)
        check("witness : double cosh at 1.5 ulps below 1", true,
              contains(A.Cosh(A.FloatNum(0x1p-30)), 1 - 3 * 0x1p-53));
        check("witness : double sin stays within 1", true, A.Sin(A.FloatNum(0x1.921fb54442d19p+0)).hi() == 1);
        programPrecision() = 1;
        check("witness : float exp at 1.984 ulps below 1/2", true,
              contains(A.Exp(A.FloatNum(-0x1.62e42ep-1)), 0x1.fffffap-2));
        check("witness : float log at 1.981 ulps below 1", true,
              contains(A.Log(A.FloatNum(0x1.5bf0aap+1)), 0x1.fffffap-1));
        programPrecision() = saved;
    }
    // a subnormal result may be flushed to zero (FTZ, DAZ) : x * 1e-20 for x in [1e-18, 1] holds
    // 0 in float, and so does its comparison to 0 ; in the corridors too
    {
        const int saved    = programPrecision();
        programPrecision() = 1;
        const interval x   = interval(1e-18, 1, -24);
        const interval y   = A.Mul(x, A.FloatNum(1e-20));
        check("witness : float x * 1e-20 may be a flushed zero", true, contains(y, 0.0));
        check("witness : float x * 1e-20 > 0 may be false", true, contains(A.Gt(y, A.FloatNum(0)), 0.0));
        const affine_algebra F10(10);
        check("witness : a float corridor x * 1e-20 may be a flushed zero", true,
              contains(toItv(F10.Mul(fromItv(x), F10.FloatNum(1e-20)), 10), 0.0));
        check("witness : a float subnormal constant may be read as 0 (DAZ)", true,
              contains(toItv(F10.FloatNum(1e-40), 10), 0.0));
        // the flush is idempotent : a recursion adding a subnormal settles
        AffItv r    = aempty();
        bool   done = false;
        for (int it = 1; !done && it <= 100; it++) {
            AffItv f = F10.Add(F10.Mul(F10.Mem(r), F10.FloatNum(0.5)), F10.FloatNum(1e-40));
            if (it > 8) f = awiden(r, f, 10);
            done = !r.isEmpty() && aleq(f, r, 10);
            r    = f;
        }
        check("witness : a recursion adding a float subnormal settles", true, done && contains(toItv(r, 10), 0.0));
        programPrecision() = saved;
    }
    // a probe seed certified over [0, T] only drifts : the image of y = select2(sel, 10 - y',
    // c), c = t, from the seed [0, 10] holds -10, the value of the program at t = 21
    {
        const affine_algebra F5(5);
        const AffItv         sel   = F5.Checkbox(F5.Label("sel"));
        const AffItv         seed  = probeSeed(interval(0, 10, -24));
        const AffItv         image = F5.Select2(sel, F5.Sub(F5.FloatNum(10), F5.Mem(seed)), AffItv{0, 0, 0, 1, -24});
        check("witness : the image of a probe seed is within it over [0, T]", true, aleqWithin(image, seed, 5));
        check("witness : the image of a probe seed holds the program past T", true, contains(forever(image, 5), -10.0));
    }
    // a bound of infinite intercept has one writing : the sum of a ramp, y = y' + t, settles
    {
        const double         T48 = 48000;
        const affine_algebra F48(T48);
        AffItv               y    = aempty();
        bool                 done = false;
        for (int it = 1; !done && it <= 100; it++) {
            AffItv f = F48.Add(F48.Mem(y), AffItv{0, 1, 0, 1, -24});
            if (it > 8) f = awiden(y, f, T48);
            done = !y.isEmpty() && aleq(f, y, T48);
            y    = f;
        }
        check("witness : the sum of a ramp settles", true, done);
    }
    // a constant -0.0 is -0 : 1 / -0 = -inf
    check("witness : the reciprocal of the constant -0.0 is -inf", true, A.Inv(A.FloatNum(-0.0)).is(-HUGE_VAL));
    // an interval built outside the algebra enters the corridors flushed
    {
        const int saved    = programPrecision();
        programPrecision() = 1;
        check("witness : fromItv flushes a subnormal bound", true, contains(toItv(fromItv(interval(1e-40, 1, -24)), 10), 0.0));
        programPrecision() = saved;
    }
    // a corridor times the constant 1 is that corridor ; y = y' * 1.0 + 0.1 settles as y' + 0.1
    {
        const AffItv c{-1, 0.25, 2, 0.5, -24};
        const AffItv p = F.Mul(c, F.FloatNum(1.0));
        check("witness : a corridor times 1.0 is itself", true, p.a0 == c.a0 && p.a1 == c.a1 && p.b0 == c.b0 && p.b1 == c.b1);
        const double         T20 = 1048576;
        const affine_algebra F20(T20);
        AffItv               y    = aempty();
        bool                 done = false;
        for (int it = 1; !done && it <= 100; it++) {
            AffItv f = F20.Add(F20.Mul(F20.Mem(y), F20.FloatNum(1.0)), F20.FloatNum(0.1));
            if (it > 8) f = awiden(y, f, T20);
            done = !y.isEmpty() && aleq(f, y, T20);
            y    = f;
        }
        check("witness : y = y' * 1.0 + 0.1 settles on a finite corridor", true, done && std::isfinite(y.b0));
    }
    // an accumulator of positive values stays positive : y = y' + t (double) keeps y >= 0
    {
        const double         T48 = 48000;
        const affine_algebra F48(T48);
        AffItv               y    = aempty();
        bool                 done = false;
        for (int it = 1; !done && it <= 100; it++) {
            AffItv f = F48.Add(F48.Mem(y), AffItv{0, 1, 0, 1, -24});
            if (it > 8) f = awiden(y, f, T48);
            done = !y.isEmpty() && aleq(f, y, T48);
            y    = f;
        }
        check("witness : the sum of a ramp stays >= 0", true, done && y.lo(0) >= 0 && y.lo(T48) >= 0);
    }
    // a product by a power of two is exact where its values stay normal : no rounding margin
    {
        const AffItv c{1, 0.25, 2, 0.5, -24};
        const AffItv h = F.Mul(c, F.FloatNum(0.5));
        const AffItv q = F.Div(c, F.FloatNum(4.0));
        check("witness : a corridor times 0.5 is exactly halved", true, h.a0 == 0.5 && h.a1 == 0.125 && h.b0 == 1 && h.b1 == 0.25);
        check("witness : a corridor over 4 is exactly quartered", true, q.a0 == 0.25 && q.a1 == 0.0625 && q.b0 == 0.5 && q.b1 == 0.125);
    }
    // an operation on a moving interval gives a moving result where the operation on its
    // values for all time passes it : c = t is [0, 10] over the horizon T = 10, its upper
    // bound rising past it ; c + 1 is 21 at t = 20, sin(c) stays within [-1, 1]
    {
        const interval c = hullWithMotion(AffItv{0, 0, 0, 1, -24}, 10);
        const interval k = A.FloatNum(0.5);
        const interval s = A.Add(c, A.FloatNum(1));
        check("witness : a moving interval plus 1 moves up", true, s.hiMoves() && !s.loMoves());
        check("witness : the negation of a rising bound falls", true, A.Neg(c).loMoves() && !A.Neg(c).hiMoves());
        check("witness : sin of a moving interval holds [-1, 1] for all time", false, A.Sin(c).moves());
        check("witness : tanh of a rising interval moves up only", true, A.Tanh(c).hiMoves() && !A.Tanh(c).loMoves());
        check("witness : a mix with a moving operand moves", true, A.Mix(k, c, k).moves());
        check("witness : an IRM of a moving input moves", true, A.Irm(c, k).moves());
        check("witness : a sum with a moving term moves", true, A.Sum({k, c}).moves());
        check("witness : a FIR of a moving input moves", true, A.Fir(c, {k, k}).moves());
        check("witness : an IIR of a moving input moves", true, A.Iir(c, {k}).moves());
        check("witness : a bpf whose input reaches a region past T moves", true, A.Bpf(c, {A.FloatNum(20)}, {k, A.FloatNum(2)}).hiMoves());
        check("witness : an operation on intervals that do not move does not move", false, A.Add(k, A.Mix(k, k, k)).moves());
        // x = 10 - t over T = 5 is [5, 10], falling past T : sqrt(x) is NaN at t = 11
        const interval x = hullWithMotion(AffItv{10, -1, 10, -1, -24}, 5);
        check("witness : sqrt of a falling interval may be invalid past T", true, A.Sqrt(x).mayBeInvalid());
        check("witness : log of a falling interval may be invalid past T", true, A.Log(x).mayBeInvalid());
    }
    // a line that reaches 0 by an exact cancellation bounds either zero : y = -(t - T) is -0 at
    // t = T, and 1 / y is -inf there
    {
        const affine_algebra F1k(1024);
        const AffItv         y = F1k.Neg(F1k.Sub(AffItv{0, 1, 0, 1, -24}, F1k.FloatNum(1024)));
        check("witness : -(t - T) holds -0 at t = T", true, leTotal(y.lo(1024), -0.0));
        check("witness : 1 / -(t - T) holds -inf", true, contains(toItv(F1k.Inv(y), 1024), -HUGE_VAL));
    }
    // a flushed corridor holds the flushed values at every sample time : the lower line
    // 2^-1023 - 2^-1024 t (T = 2) is the subnormal 2^-1024 at t = 1, read as 0 ; the flush is
    // idempotent, and a corridor of ordinary values is left as it is
    {
        const double   m  = 0x1p-1022;
        const AffItv   X  = flushedToZero(AffItv{0x1p-1023, -0x1p-1024, 1, 0, -24}, 2);
        const AffItv   X2 = flushedToZero(X, 2);
        const AffItv   c{0.5, -0.25, 1, 0.5, -24};
        const AffItv   fc = flushedToZero(c, 2);
        check("witness : a flushed lower line holds 0 where it was subnormal", true, X.lo(1) <= 0);
        check("witness : the flush of a corridor is idempotent", true, X2.a0 == X.a0 && X2.a1 == X.a1 && X2.b0 == X.b0 && X2.b1 == X.b1);
        check("witness : the flush leaves an ordinary corridor as it is", true, fc.a0 == c.a0 && fc.a1 == c.a1 && fc.b0 == c.b0 && fc.b1 == c.b1);
        // the upper line from -m to +0 (T = 2), times 1/4 : the chord ends at +0, not -0
        const affine_algebra F2(2);
        const AffItv         R = F2.Mul(AffItv{-1, 0, -m, m / 2, -24}, F2.FloatNum(0.25));
        check("witness : a product that reaches +0 holds +0", true, leTotal(0.0, R.hi(2)));
    }
    return reportCheckResults();
}
