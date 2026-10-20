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
    // a wrapped operand holds any int32 : the operation may wrap, even adding zero
    const AffItv wrapped{0, 0, 0, 2, 0};
    check("affine add : a wrapped operand", true, F.Add(wrapped, C(0, 0)).mayOverflow);

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

}  // namespace

int main()
{
    interval_algebra A;
    addTests(A);
    subTests(A);
    negTests(A);
    mulTests(A);
    sumTests(A);
    floatTests(A);
    intCastTests(A);
    affineTests();
    return reportCheckResults();
}
