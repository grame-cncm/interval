// The int32 operations of the program wrap in two's complement : the interval of a result is
// the set of the values the program computes, after the wrap. Each operation is checked
// without overflow, at the exact edge of the int32 range, crossing the range by either end
// (the whole range), and with a certain overflow (the range shifted by 2^32). A floating
// operation never wraps. The int conversion saturates the image of the valid inputs.

#include <climits>
#include <string>

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

const interval kAll = I(m, M);

void addTests(const interval_algebra& A)
{
    checkExact("add : no overflow", A.Add(I(0, 10), I(0, 10)), I(0, 20));
    checkExact("add : reaches INT_MAX exactly", A.Add(I(M - 10, M - 5), I(0, 5)), I(M - 10, M));
    checkExact("add : reaches INT_MIN exactly", A.Add(I(m + 5, m + 10), I(-5, 0)), I(m, m + 10));
    checkExact("add : crosses INT_MAX", A.Add(I(M - 10, M), I(0, 1)), kAll);
    checkExact("add : crosses INT_MIN", A.Add(I(m, m + 5), I(-1, 0)), kAll);
    checkExact("add : certain overflow above", A.Add(I(M - 1, M), I(10, 20)), I(m + 8, m + 19));
    checkExact("add : certain overflow below", A.Add(I(m, m + 1), I(-20, -10)),
               I(M - 19, M - 8));
    checkExact("add : the whole range plus zero", A.Add(kAll, I(0, 0)), kAll);
    checkExact("add : the whole range plus one", A.Add(kAll, I(1, 1)), kAll);
}

void subTests(const interval_algebra& A)
{
    checkExact("sub : no overflow", A.Sub(I(0, 10), I(0, 10)), I(-10, 10));
    checkExact("sub : reaches INT_MIN exactly", A.Sub(I(m + 5, m + 10), I(0, 5)), I(m, m + 10));
    checkExact("sub : crosses INT_MIN", A.Sub(I(m, 0), I(0, 1)), kAll);
    checkExact("sub : crosses INT_MAX", A.Sub(I(0, M), I(-1, 0)), kAll);
    checkExact("sub : 0 - INT_MIN wraps to INT_MIN", A.Sub(I(0, 0), I(m, m)), I(m, m));
    checkExact("sub : INT_MAX - (-1) wraps to INT_MIN", A.Sub(I(M, M), I(-1, -1)), I(m, m));
}

void negTests(const interval_algebra& A)
{
    checkExact("neg : no overflow", A.Neg(I(-5, 3)), I(-3, 5));
    checkExact("neg : reaches INT_MAX exactly", A.Neg(I(m + 1, 0)), I(0, M));
    checkExact("neg : -INT_MIN within the operand", A.Neg(I(m, 0)), kAll);
    checkExact("neg : -INT_MIN is INT_MIN", A.Neg(I(m, m)), I(m, m));
}

void mulTests(const interval_algebra& A)
{
    checkExact("mul : no overflow", A.Mul(I(-3, 4), I(-5, 6)), I(-20, 24));
    checkExact("mul : below INT_MAX", A.Mul(I(0, 46340), I(0, 46340)), I(0, 2147395600.0));
    checkExact("mul : INT_MAX times one", A.Mul(I(M, M), I(1, 1)), I(M, M));
    checkExact("mul : reaches INT_MIN exactly", A.Mul(I(-65536, -65536), I(32768, 32768)),
               I(m, m));
    checkExact("mul : crosses INT_MAX", A.Mul(I(0, 46341), I(0, 46341)), kAll);
    check("mul : -1 times INT_MIN wraps to INT_MIN", true,
          A.Mul(I(-1, -1), I(m, m)).has(m));
}

void sumTests(const interval_algebra& A)
{
    // the partial sums from left to right : INT_MAX + 1 wraps, then - 1 comes back
    checkExact("sum : a partial sum wraps, the sum does not",
               A.Sum({I(M, M), I(1, 1), I(-1, -1)}), I(M, M));
    checkExact("sum : no partial sum wraps", A.Sum({I(0, 10), I(0, 10), I(-5, 0)}), I(-5, 20));
}

void floatTests(const interval_algebra& A)
{
    // a floating operation never wraps, whatever its magnitude
    const interval s = A.Add(interval(2147483000.5, 2147483647.5), interval(100, 100));
    check("float add : beyond INT_MAX, no wrap", true, s.lo() > 2147483000.0 && s.hi() > M);
    const interval h = A.Add(interval(0, 1e300), interval(0, 1e300));
    check("float add : huge, no wrap", true, h.lo() == 0 && h.hi() >= 2e300);
}

void intCastTests(const interval_algebra& A)
{
    checkExact("intcast : truncation", A.IntCast(interval(-3.8, 4.9)), I(-3, 4));
    const interval up = A.IntCast(interval(0, 3e9));
    check("intcast : saturates above", true, up.lo() == 0 && up.hi() == M);
    check("intcast : beyond INT_MAX may be invalid", true, up.mayBeInvalid());
    const interval down = A.IntCast(interval(-3e9, 0));
    check("intcast : saturates below", true, down.lo() == m && down.hi() == 0);
    check("intcast : below INT_MIN may be invalid", true, down.mayBeInvalid());
    checkExact("intcast : an int is unchanged", A.IntCast(kAll), kAll);
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
    return reportCheckResults();
}
