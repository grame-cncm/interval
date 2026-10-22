// Each check computes an interval and a value the compiled program produces : the value must
// be in the interval (the soundness criterion Y within Y^ of the README). One check per
// witness of a defect found by a review of the library.

#include <climits>
#include <cmath>
#include <cstdint>

#include "interval/affine_ops.hh"
#include "interval/affint.hh"
#include "interval/check.hh"
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
    return reportCheckResults();
}
