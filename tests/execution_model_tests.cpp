// The execution model of the README ("In short"), one rule at a time : every value that an
// execution may produce lies in the interval of the computation.
//   1. basic operations round to nearest in the program's precision ;
//   2. libm functions are within 2 ulps of the exact result ;
//   3. the C++ compiler may contract a * b + c into a fused multiply-add.
// Random operands with a fixed seed, in single and in double precision.

#include <cmath>
#include <cstdint>
#include <random>

#include "interval/check.hh"
#include "interval/interval_algebra.hh"
#include "interval/interval_def.hh"

using namespace itv;

namespace {

bool contains(const interval& r, double v)
{
    return r.lo() <= v && v <= r.hi();
}

// a value of the program : a double, or a float in single precision
double programValue(double x, int precision)
{
    return precision == 1 ? double(float(x)) : x;
}

// the value one step away in the program's precision
double next(double x, int precision, double dir)
{
    return precision == 1 ? double(std::nextafter(float(x), float(dir))) : std::nextafter(x, dir);
}

const char* name(int precision)
{
    return precision == 1 ? "single" : "double";
}

}  // namespace

int main()
{
    std::mt19937_64                        rng(20261023);
    std::uniform_real_distribution<double> mantissa(-1.0, 1.0);
    std::uniform_int_distribution<int>     exponent(-20, 20);
    auto draw = [&](int precision) { return programValue(std::ldexp(mantissa(rng), exponent(rng)), precision); };

    for (int precision : {1, 2}) {
        programPrecision() = precision;
        const interval_algebra A;

        // rule 1 : + - * / sqrt, each rounded to nearest in the program's precision
        bool rounded = true;
        for (int i = 0; i < 4000 && rounded; i++) {
            const double a = draw(precision), b = draw(precision);
            const interval ia = A.FloatNum(a), ib = A.FloatNum(b);
            volatile float  fa = float(a), fb = float(b);
            volatile double da = a, db = b;
            const double add = precision == 1 ? double(fa + fb) : da + db;
            const double sub = precision == 1 ? double(fa - fb) : da - db;
            const double mul = precision == 1 ? double(fa * fb) : da * db;
            const double div = precision == 1 ? double(fa / fb) : da / db;
            const double sq  = precision == 1 ? double(std::sqrt(std::fabs(fa))) : std::sqrt(std::fabs(da));
            rounded = contains(A.Add(ia, ib), add) && contains(A.Sub(ia, ib), sub) &&
                      contains(A.Mul(ia, ib), mul) && (b == 0 || contains(A.Div(ia, ib), div)) &&
                      contains(A.Sqrt(A.Abs(ia)), sq);
        }
        check(std::string("model : ") + name(precision) + " basic operations rounded to nearest", true, rounded);

        // rule 1, at compile time : a literal is a value of the program (0.1 is 0.1f in single)
        const double tenth = programValue(0.1, precision);
        check(std::string("model : ") + name(precision) + " literal 0.1 is the program's 0.1", true,
              A.FloatNum(0.1).lo() == tenth && A.FloatNum(0.1).hi() == tenth);

        // rule 2 : sin, exp and log within 2 ulps of the exact result. In single precision the
        // exact result is taken in double, and every float within 2 ulps of it is checked ; in
        // double precision the libm result and its neighbours, within 2 ulps of the exact one
        bool libm = true;
        for (int i = 0; i < 2000 && libm; i++) {
            const double x = draw(precision);
            const double e = std::fabs(x) > 80 ? 0.0 : x;  // exp of a moderate argument
            struct Case {
                interval r;
                double   exact;
            } cases[] = {{A.Sin(A.FloatNum(x)), std::sin(x)},
                         {A.Exp(A.FloatNum(e)), std::exp(e)},
                         {A.Log(A.FloatNum(std::fabs(x))), std::log(std::fabs(x))}};
            for (const Case& c : cases) {
                if (!std::isfinite(c.exact)) continue;
                if (precision == 1) {
                    const float  r   = float(c.exact);
                    const double ulp = double(std::nextafter(std::fabs(float(c.exact)), INFINITY)) - std::fabs(double(float(c.exact)));
                    float        v   = std::nextafter(std::nextafter(r, -INFINITY), -INFINITY);
                    for (int k = 0; k < 5; k++, v = std::nextafter(v, INFINITY)) {
                        if (std::fabs(double(v) - c.exact) <= 2 * ulp) libm = libm && contains(c.r, double(v));
                    }
                } else {
                    libm = libm && contains(c.r, c.exact) && contains(c.r, next(c.exact, 2, -HUGE_VAL)) &&
                           contains(c.r, next(c.exact, 2, HUGE_VAL));
                }
            }
        }
        check(std::string("model : ") + name(precision) + " libm results within 2 ulps", true, libm);

        // rule 3 : a * b + c, rounded twice or fused
        bool fused = true;
        for (int i = 0; i < 4000 && fused; i++) {
            const double a = draw(precision), b = draw(precision), c = draw(precision);
            const interval r = A.Add(A.Mul(A.FloatNum(a), A.FloatNum(b)), A.FloatNum(c));
            volatile float  fa = float(a), fb = float(b), fc = float(c);
            volatile double da = a, db = b, dc = c;
            const double twice = precision == 1 ? double(float(fa * fb) + fc) : (da * db) + dc;
            const double once  = precision == 1 ? double(std::fma(float(fa), float(fb), float(fc)))
                                                : std::fma(double(da), double(db), double(dc));
            fused = contains(r, twice) && contains(r, once);
        }
        check(std::string("model : ") + name(precision) + " a * b + c, rounded twice or fused", true, fused);
    }
    programPrecision() = 2;
    return reportCheckResults();
}
