// The int32 operations of the program wrap in two's complement : the interval of a result is
// the set of the values the program computes, after the wrap, and its flag mayOverflow()
// says whether the operation may wrap. Each operation is checked without overflow, at the
// exact edge of the int32 range, crossing the range by either end (the whole range), and
// with a certain overflow (the range shifted by 2^32). A floating operation never wraps.
// The int conversion saturates the image of the valid inputs.

#include <climits>
#include <string>

#include "interval/affine_ops.hh"
#include "interval/affint.hh"
#include "interval/check.hh"
#include "interval/interval_algebra.hh"
#include "interval/interval_def.hh"

using namespace itv;

namespace {

const double M = INT_MAX;
const double m = INT_MIN;

interval I(double lo, double hi)
{
    return interval(lo, hi, 0);
}

// an int32 interval computed by an operation that may wrap
interval W(double lo, double hi)
{
    return I(lo, hi).withOverflow();
}

const interval kAll  = I(m, M);
const interval kWrap = W(m, M);

void addTests(const interval_algebra& A)
{
    checkExact("add : no overflow", A.Add(I(0, 10), I(0, 10)), I(0, 20));
    checkExact("add : reaches INT_MAX exactly", A.Add(I(M - 10, M - 5), I(0, 5)), I(M - 10, M));
    checkExact("add : reaches INT_MIN exactly", A.Add(I(m + 5, m + 10), I(-5, 0)), I(m, m + 10));
    checkExact("add : crosses INT_MAX", A.Add(I(M - 10, M), I(0, 1)), kWrap);
    checkExact("add : crosses INT_MIN", A.Add(I(m, m + 5), I(-1, 0)), kWrap);
    checkExact("add : certain overflow above", A.Add(I(M - 1, M), I(10, 20)), W(m + 8, m + 19));
    checkExact("add : certain overflow below", A.Add(I(m, m + 1), I(-20, -10)),
               W(M - 19, M - 8));
    checkExact("add : the whole range plus zero", A.Add(kAll, I(0, 0)), kAll);
    checkExact("add : the whole range plus one", A.Add(kAll, I(1, 1)), kWrap);
    checkExact("add : the flag of an operand is not the flag of the sum",
               A.Add(kWrap, I(0, 0)), kAll);
}

void subTests(const interval_algebra& A)
{
    checkExact("sub : no overflow", A.Sub(I(0, 10), I(0, 10)), I(-10, 10));
    checkExact("sub : reaches INT_MIN exactly", A.Sub(I(m + 5, m + 10), I(0, 5)), I(m, m + 10));
    checkExact("sub : crosses INT_MIN", A.Sub(I(m, 0), I(0, 1)), kWrap);
    checkExact("sub : crosses INT_MAX", A.Sub(I(0, M), I(-1, 0)), kWrap);
    checkExact("sub : 0 - INT_MIN wraps to INT_MIN", A.Sub(I(0, 0), I(m, m)), W(m, m));
    checkExact("sub : INT_MAX - (-1) wraps to INT_MIN", A.Sub(I(M, M), I(-1, -1)), W(m, m));
}

void negTests(const interval_algebra& A)
{
    checkExact("neg : no overflow", A.Neg(I(-5, 3)), I(-3, 5));
    checkExact("neg : reaches INT_MAX exactly", A.Neg(I(m + 1, 0)), I(0, M));
    checkExact("neg : -INT_MIN within the operand", A.Neg(I(m, 0)), kWrap);
    checkExact("neg : -INT_MIN is INT_MIN", A.Neg(I(m, m)), W(m, m));
}

void mulTests(const interval_algebra& A)
{
    checkExact("mul : no overflow", A.Mul(I(-3, 4), I(-5, 6)), I(-20, 24));
    checkExact("mul : below INT_MAX", A.Mul(I(0, 46340), I(0, 46340)), I(0, 2147395600.0));
    checkExact("mul : INT_MAX times one", A.Mul(I(M, M), I(1, 1)), I(M, M));
    checkExact("mul : reaches INT_MIN exactly", A.Mul(I(-65536, -65536), I(32768, 32768)),
               I(m, m));
    checkExact("mul : crosses INT_MAX", A.Mul(I(0, 46341), I(0, 46341)), kWrap);
    const interval p = A.Mul(I(-1, -1), I(m, m));
    check("mul : -1 times INT_MIN wraps to INT_MIN", true, p.has(m) && p.mayOverflow());
}

// the int32 power is computed by products : it wraps where they would
void powTests(const interval_algebra& A)
{
    check("pow : 2^30 does not wrap", false, A.Pow(I(2, 2), I(30, 30)).mayOverflow());
    check("pow : 2^31 may wrap", true, A.Pow(I(2, 2), I(31, 31)).mayOverflow());
    check("pow : [0, 100000]^2 may wrap", true, A.Pow(I(0, 100000), I(2, 2)).mayOverflow());
    check("pow : [0, 1000]^2 does not wrap", false, A.Pow(I(0, 1000), I(2, 2)).mayOverflow());
    check("pow : a float power never wraps", false, A.Pow(interval(2, 2, -24), interval(40, 40, 0)).mayOverflow());
    // an integer recursion is computed by sums and products that wrap
    check("iir : an integer recursion may wrap", true, A.Iir(I(0, 1), {I(1, 1)}).mayOverflow());
    check("iir : a float recursion does not wrap", false,
          A.Iir(interval(0, 1, -24), {interval(0.5, 0.5, -24)}).mayOverflow());
}

void sumTests(const interval_algebra& A)
{
    // the partial sums from left to right : INT_MAX + 1 wraps, then - 1 comes back
    checkExact("sum : a partial sum wraps, the sum does not",
               A.Sum({I(M, M), I(1, 1), I(-1, -1)}), W(M, M));
    checkExact("sum : no partial sum wraps", A.Sum({I(0, 10), I(0, 10), I(-5, 0)}), I(-5, 20));
    checkExact("prod : a partial product wraps", A.Prod({I(65536, 65536), I(65536, 65536),
                                                         I(0, 0)}),
               W(0, 0));
    checkExact("prod : no partial product wraps", A.Prod({I(-3, 4), I(2, 2), I(10, 10)}),
               I(-60, 80));
}

void floatTests(const interval_algebra& A)
{
    // a floating operation never wraps, whatever its magnitude
    const interval s = A.Add(interval(2147483000.5, 2147483647.5), interval(100, 100));
    check("float add : beyond INT_MAX, no wrap", true,
          s.lo() > 2147483000.0 && s.hi() > M && !s.mayOverflow());
    const interval h = A.Add(interval(0, 1e300), interval(0, 1e300));
    check("float add : huge, no wrap", true, h.lo() == 0 && h.hi() >= 2e300 && !h.mayOverflow());
    const interval f = A.Mul(interval(0, 1e200), interval(0, 1e200));
    check("float mul : huge, no wrap", true, !f.mayOverflow());
}

void intCastTests(const interval_algebra& A)
{
    checkExact("intcast : truncation", A.IntCast(interval(-3.8, 4.9)), I(-3, 4));
    const interval up = A.IntCast(interval(0, 3e9));
    check("intcast : saturates above, no wrap", true,
          up.lo() == 0 && up.hi() == M && !up.mayOverflow());
    check("intcast : beyond INT_MAX may be invalid", true, up.mayBeInvalid());
    const interval down = A.IntCast(interval(-3e9, 0));
    check("intcast : saturates below", true, down.lo() == m && down.hi() == 0);
    check("intcast : below INT_MIN may be invalid", true, down.mayBeInvalid());
    checkExact("intcast : an int is unchanged", A.IntCast(kAll), kAll);
}

// The affine domain, over the horizon of 2^31 samples, as the Faust typing reads it
void affineTests()
{
    const double         T = 2147483648.0;
    const affine_algebra F(T);
    auto                 C = [](double lo, double hi) { return fromItv(I(lo, hi)); };

    check("affine add : no overflow", false, F.Add(C(0, 10), C(0, 10)).mayOverflow);
    check("affine add : crosses INT_MAX", true, F.Add(C(M - 10, M), C(0, 1)).mayOverflow);
    check("affine sub : crosses INT_MIN", true, F.Sub(C(m, 0), C(0, 1)).mayOverflow);
    check("affine neg : -INT_MIN", true, F.Neg(C(m, 0)).mayOverflow);
    check("affine neg : no overflow", false, F.Neg(C(m + 1, 0)).mayOverflow);
    check("affine mul : crosses INT_MAX", true, F.Mul(C(0, 46341), C(0, 46341)).mayOverflow);
    check("affine mul : no overflow", false, F.Mul(C(0, 46340), C(0, 46340)).mayOverflow);

    // a unit counter [0, t] passes INT_MAX within the horizon : its next value may wrap
    const AffItv counter{0, 0, 0, 1, 0};
    check("affine add : a unit counter plus one", true, F.Add(counter, C(1, 1)).mayOverflow);
    check("affine mul : a unit counter times two", true, F.Mul(counter, C(2, 2)).mayOverflow);
    // a wrapped operand holds any int32 : adding one may wrap, adding zero cannot
    const AffItv wrapped{0, 0, 0, 2, 0};
    check("affine add : a wrapped operand plus one", true, F.Add(wrapped, C(1, 1)).mayOverflow);
    check("affine add : a wrapped operand plus zero", false, F.Add(wrapped, C(0, 0)).mayOverflow);

    check("affine sum : a partial sum wraps", true,
          F.Sum({C(M, M), C(1, 1), C(-1, -1)}).mayOverflow);
    check("affine sum : no partial sum wraps", false,
          F.Sum({C(0, 10), C(0, 10), C(-5, 0)}).mayOverflow);

    // the flag crosses between the domains and survives the joins
    const AffItv flagged = C(0, 10).withOverflow();
    check("toItv keeps the flag", true, toItv(flagged, T).mayOverflow());
    check("fromItv keeps the flag", true, fromItv(I(0, 10).withOverflow()).mayOverflow);
    check("a join keeps the flag", true, ajoin(flagged, C(20, 30), T).mayOverflow);
    check("a flagged corridor is not below an unflagged one", false,
          aleq(flagged, C(0, 10), T));
    check("an unflagged corridor is below a flagged one", true, aleq(C(0, 10), flagged, T));
}

// A corridor that moves bounds its values up to the horizon only : an int32 operation may
// wrap when it may on the values of its operands for all time, and a nonlinear operation of
// moving operands drifts past the horizon unless it bounds them.
void growthTests()
{
    const double         T = 2147483648.0;
    const affine_algebra F(T);
    auto                 C = [](double lo, double hi) { return fromItv(I(lo, hi)); };

    // an int growing by a quarter per sample, [0, 2^29] at the horizon
    const AffItv slow{0, 0, 0, 0.25, 0};
    check("growth : a slow int counter plus one", true, F.Add(slow, C(1, 1)).mayOverflow);
    check("growth : a slow int counter minus one", false, F.Sub(slow, C(1, 1)).mayOverflow);
    check("growth : a slow int counter negated", false, F.Neg(slow).mayOverflow);
    check("growth : a slow int counter times two", true, F.Mul(slow, C(2, 2)).mayOverflow);

    // the int conversion of a float accumulator keeps its growth
    const AffItv acc{0, 0, 0, 0.5, -24};
    const AffItv n = F.IntCast(acc);
    check("growth : int of a float accumulator moves", true, n.moving());
    check("growth : int of a float accumulator plus one", true, F.Add(n, C(1, 1)).mayOverflow);

    // a modulo, a comparison, a sine bound a moving operand for all time
    const AffItv r = F.Mod(slow, C(100, 100));
    check("growth : a counter modulo 100 does not move", false, r.moving());
    check("growth : a counter modulo 100, plus one", false, F.Add(r, C(1, 1)).mayOverflow);
    check("growth : a comparison does not move", false, F.Lt(slow, C(10, 10)).moving());
    check("growth : the sine of an accumulator does not move", false, F.Sin(acc).moving());

    // a product of two moving operands is nonlinear : it drifts, and may wrap
    const AffItv tiny{0, 0, 0, 0x1p-20, 0};  // [0, 2048] at the horizon
    const AffItv sq = F.Mul(tiny, tiny);
    check("growth : a product of two counters drifts", true, sq.hiDrifts);
    check("growth : a product of two counters may wrap", true, sq.mayOverflow);
    check("growth : a drifting int plus one", true, F.Add(sq, C(1, 1)).mayOverflow);
    const AffItv gentle{0, 0, 0, 0x1p-30, -24};  // [0, 2] at the horizon
    check("growth : the exponential of an accumulator drifts", true, F.Exp(gentle).hiDrifts);
    check("growth : the exponential of a saturated one does not", false, F.Exp(acc).hiDrifts);

    // a drift survives a delay and a join, and is part of the order
    check("growth : a delayed drifting value drifts", true, F.Delay(sq, C(1, 1)).hiDrifts);
    check("growth : a drifting corridor is not below a fixed one", false,
          aleq(sq, C(-1e10, 1e10), T));
    // a fixed corridor stays valid for all time
    check("growth : a fixed corridor does not move", false, F.Add(C(0, 10), C(0, 10)).moving());
    // an int32 value that may wrap may hold any int32 value once it has wrapped : it moves on
    // both sides, and the operations after it, through a delay, may wrap in turn
    const AffItv w = F.Add(slow, C(1, 1));
    check("wrap : a counter plus one that may wrap drifts on both sides", true, w.loDrifts && w.hiDrifts);
    check("wrap : then minus one may wrap", true, F.Sub(w, C(1, 1)).mayOverflow);
    check("wrap : then delayed, minus two may wrap", true, F.Sub(F.Delay(w, C(1, 1)), C(2, 2)).mayOverflow);
    check("wrap : then negated may wrap", true, F.Neg(w).mayOverflow);
    const AffItv k = F.Add(C(M, M), C(1, 1));
    check("wrap : a constant that wraps moves", true, k.mayOverflow && k.moving());
    check("wrap : a fixed corridor that cannot wrap does not drift", false, F.Add(C(0, 10), C(0, 10)).loDrifts);

    // a chord that does not rise while a bound it covers keeps rising : the bound drifts.
    // The join of a rising and a falling bound : flat at 10 over [0, T], unbounded past T.
    const AffItv up{0, 0, 0, 10 / T, 0}, down{0, 0, 10, -10 / T, 0};
    const AffItv j = ajoin(up, down, T);
    check("chord : the join of a rising and a falling bound drifts", true, j.hiDrifts);
    check("chord : its flat chord bounds [0, T]", true, j.hi(0) == 10 && j.hi(T) == 10);
    // a float crossing zero, times a constant of either sign : the product shrinks over
    // [0, T], then grows again past it
    const AffItv cross{-1610612736.0, 1, -1610612736.0, 1, -24};
    const AffItv g = F.Mul(cross, fromItv(interval(-1, 1, -24)));
    check("chord : a product crossing zero by a signed constant drifts", true, g.loDrifts && g.hiDrifts);
    // a rising chord needs no drift : its line already leaves every bound past T
    const AffItv twice = F.Mul(acc, fromItv(interval(2, 2, -24)));
    check("chord : a rising float product moves", true, twice.moving());
    check("chord : a rising float product does not drift", false, twice.loDrifts || twice.hiDrifts);

    // a probe seed claims a bound over [0, T] : its image may drift past it
    AffItv seed = C(0, 16), image = C(0, 8);
    image.hiDrifts = true;
    check("probe : a drifting image is not below a fixed seed", false, aleq(image, seed, T));
    check("probe : it is within the seed over [0, T]", true, aleqWithin(image, seed, T));
    check("probe : a larger image is not within the seed", false, aleqWithin(C(0, 32), seed, T));

    // a corridor rebuilt from hulls : between moving values it drifts where they move ; as
    // the meet of two bounds, where the bound that sets it moves
    const AffItv fixed = C(0, 10), rising{0, 0, 0, 1, -24};
    const AffItv h = hullDrift(C(0, 20), {fixed, rising});
    check("rebuild : between a fixed and a rising value, the upper bound drifts", true, h.hiDrifts && !h.loDrifts);
    const AffItv mt = meetDrift(C(0, 10), fixed, rising, T);
    check("rebuild : the meet set by a fixed bound does not drift", false, mt.hiDrifts);
    check("rebuild : the meet of two rising bounds drifts", true, meetDrift(C(0, 10), rising, rising, T).hiDrifts);
    // set by the rising bound (min(t, 100) is [0, 10] over [0, 10]), the meet rises with it
    {
        const double         T10 = 10;
        const affine_algebra F10(T10);
        const AffItv         wide = C(0, 100), climbing = F10.Min(rising, F10.FloatNum(100));
        check("rebuild : the meet set by a rising bound drifts", true,
              meetDrift(fromItv(intersection(toItv(wide, T10), toItv(climbing, T10))), wide, climbing, T10).hiDrifts);
    }

    // a coefficient null over [0, T] that moves past it still counts
    AffItv zeroNow = C(0, 0);
    zeroNow.hiDrifts = true;
    const AffItv fir = F.Fir(C(1, 1), {C(1, 1), zeroNow});
    check("fir : a null coefficient that drifts makes the sum drift", true, fir.hiDrifts);
    check("fir : a null coefficient that does not drift adds nothing", false,
          F.Fir(C(1, 1), {C(1, 1), C(0, 0)}).hiDrifts);

    // the motion of the bounds of a hull : each bound on its own
    const interval accm = hullWithMotion(acc, T);
    check("motion : a float accumulator, its lower bound holds", false, accm.loMoves());
    check("motion : a float accumulator, its upper bound moves", true, accm.hiMoves());
    const interval wm = hullWithMotion(w, T);
    check("motion : a counter that may wrap moves on both sides", true, wm.loMoves() && wm.hiMoves());
    check("motion : a fixed corridor does not move", false, hullWithMotion(C(0, 10), T).moves());
    check("motion : a wrapped counter, the whole int32 range, does not move", false,
          hullWithMotion(F.Add(AffItv{0, 0, 0, 1, 0}, C(1, 1)), T).moves());
    check("motion : the reunion moves where either moves", true,
          reunion(I(0, 1), I(0, 1).withMotion(false, true)).hiMoves());
    check("motion : the intersection set by a fixed bound does not move", false,
          intersection(I(0, 1), I(0, 2).withMotion(false, true)).hiMoves());
    check("motion : the intersection set by a moving bound moves", true,
          intersection(I(0, 1).withMotion(false, true), I(0, 2)).hiMoves());
}

}  // namespace

int main()
{
    interval_algebra A;
    addTests(A);
    subTests(A);
    negTests(A);
    mulTests(A);
    powTests(A);
    sumTests(A);
    floatTests(A);
    intCastTests(A);
    affineTests();
    growthTests();
    return reportCheckResults();
}
