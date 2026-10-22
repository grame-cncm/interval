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
    return reportCheckResults();
}
