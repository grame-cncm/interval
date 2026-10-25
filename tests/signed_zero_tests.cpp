// Signed zeros : an interval's bounds are read in the total order of IEEE 754, where -0 lies
// before +0 : a lower bound -0 holds -0 (and +0), a lower bound +0 holds +0 only ; an upper
// bound +0 holds +0 (and -0), an upper bound -0 holds -0 only. Random intervals whose bounds
// may be signed zeros, subnormals or infinities ; each operation is evaluated in double on
// values of its operands, its subnormal results and operands possibly flushed to a zero of
// their sign (FTZ, DAZ), and every value must lie in the interval computed, in that order.

#include <cmath>
#include <cstdio>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "interval/affine_ops.hh"
#include "interval/check.hh"
#include "interval/interval_algebra.hh"
#include "interval/interval_def.hh"

using namespace itv;

namespace {

// a <= b in the total order of the reals with -0 before +0 (NaN excluded)
bool orderedTotal(double a, double b)
{
    if (a == 0 && b == 0) return std::signbit(a) || !std::signbit(b);
    return a <= b;
}

bool holds(const interval& r, double v)
{
    if (std::isnan(v)) return r.mayBeInvalid();
    if (r.isEmpty()) return false;
    return orderedTotal(r.lo(), v) && orderedTotal(v, r.hi());
}

// the zero of the sign of a subnormal x : what a flush to zero makes of it
double flushed(double x)
{
    return (x != 0 && std::abs(x) < 0x1p-1022) ? std::copysign(0.0, x) : x;
}

const double kTiny = 0x1p-1074;
const std::vector<double> kBounds = {-HUGE_VAL, -3, -1, -0.5, -kTiny, -0.0, 0.0, kTiny, 0.5, 1, 3, HUGE_VAL};

// values of [lo, hi] : its bounds, the zeros it holds, a point inside
std::vector<double> samples(double lo, double hi)
{
    std::vector<double> v;
    for (double c : kBounds) {
        if (orderedTotal(lo, c) && orderedTotal(c, hi) && std::isfinite(c)) v.push_back(c);
    }
    if (std::isfinite(lo)) v.push_back(lo);
    if (std::isfinite(hi)) v.push_back(hi);
    if (std::isfinite(lo) && std::isfinite(hi) && lo < hi) v.push_back(lo + (hi - lo) / 3);
    return v;
}

struct Unary {
    const char* name;
    double (*f)(double);
    interval (interval_algebra::*op)(const interval&) const;
    AffItv (affine_algebra::*aop)(const AffItv&) const;
};
struct Binary {
    const char* name;
    double (*f)(double, double);
    interval (interval_algebra::*op)(const interval&, const interval&) const;
    AffItv (affine_algebra::*aop)(const AffItv&, const AffItv&) const;
};

double neg(double x) { return -x; }
double inv(double x) { return 1 / x; }
double add(double x, double y) { return x + y; }
double sub(double x, double y) { return x - y; }
double mul(double x, double y) { return x * y; }
double dvd(double x, double y) { return x / y; }
double mn(double x, double y) { return std::min(x, y); }
double mx(double x, double y) { return std::max(x, y); }

}  // namespace

int main()
{
    programPrecision() = 2;
    const interval_algebra A;
    const affine_algebra   F(10);  // constant corridors over [0, 10]
    const std::vector<Unary> unary = {
        {"neg", neg, &interval_algebra::Neg, &affine_algebra::Neg},   {"abs", [](double x) { return std::abs(x); }, &interval_algebra::Abs, &affine_algebra::Abs},
        {"inv", inv, &interval_algebra::Inv, &affine_algebra::Inv},   {"sqrt", [](double x) { return std::sqrt(x); }, &interval_algebra::Sqrt, &affine_algebra::Sqrt},
        {"ceil", [](double x) { return std::ceil(x); }, &interval_algebra::Ceil, &affine_algebra::Ceil},
        {"floor", [](double x) { return std::floor(x); }, &interval_algebra::Floor, &affine_algebra::Floor},
        {"rint", [](double x) { return std::rint(x); }, &interval_algebra::Rint, &affine_algebra::Rint},
        {"round", [](double x) { return std::round(x); }, &interval_algebra::Round, &affine_algebra::Round},
        {"exp", [](double x) { return std::exp(x); }, &interval_algebra::Exp, &affine_algebra::Exp},
        {"log", [](double x) { return std::log(x); }, &interval_algebra::Log, &affine_algebra::Log},
        {"log10", [](double x) { return std::log10(x); }, &interval_algebra::Log10, &affine_algebra::Log10},
        {"sin", [](double x) { return std::sin(x); }, &interval_algebra::Sin, &affine_algebra::Sin},
        {"cos", [](double x) { return std::cos(x); }, &interval_algebra::Cos, &affine_algebra::Cos},
        {"tan", [](double x) { return std::tan(x); }, &interval_algebra::Tan, &affine_algebra::Tan},
        {"asin", [](double x) { return std::asin(x); }, &interval_algebra::Asin, &affine_algebra::Asin},
        {"acos", [](double x) { return std::acos(x); }, &interval_algebra::Acos, &affine_algebra::Acos},
        {"atan", [](double x) { return std::atan(x); }, &interval_algebra::Atan, &affine_algebra::Atan},
        {"sinh", [](double x) { return std::sinh(x); }, &interval_algebra::Sinh, &affine_algebra::Sinh},
        {"cosh", [](double x) { return std::cosh(x); }, &interval_algebra::Cosh, &affine_algebra::Cosh},
        {"tanh", [](double x) { return std::tanh(x); }, &interval_algebra::Tanh, &affine_algebra::Tanh},
        {"asinh", [](double x) { return std::asinh(x); }, &interval_algebra::Asinh, &affine_algebra::Asinh},
        {"acosh", [](double x) { return std::acosh(x); }, &interval_algebra::Acosh, &affine_algebra::Acosh},
        {"atanh", [](double x) { return std::atanh(x); }, &interval_algebra::Atanh, &affine_algebra::Atanh},
    };
    const std::vector<Binary> binary = {
        {"add", add, &interval_algebra::Add, &affine_algebra::Add}, {"sub", sub, &interval_algebra::Sub, &affine_algebra::Sub},
        {"mul", mul, &interval_algebra::Mul, &affine_algebra::Mul}, {"div", dvd, &interval_algebra::Div, &affine_algebra::Div},
        {"min", mn, &interval_algebra::Min, &affine_algebra::Min},  {"max", mx, &interval_algebra::Max, &affine_algebra::Max},
        {"pow", [](double x, double y) { return std::pow(x, y); }, &interval_algebra::Pow, &affine_algebra::Pow},
        {"fmod", [](double x, double y) { return std::fmod(x, y); }, &interval_algebra::Fmod, &affine_algebra::Fmod},
        {"mod", [](double x, double y) { return std::fmod(x, y); }, &interval_algebra::Mod, &affine_algebra::Mod},
        {"remainder", [](double x, double y) { return std::remainder(x, y); }, &interval_algebra::Remainder, &affine_algebra::Remainder},
        {"atan2", [](double x, double y) { return std::atan2(x, y); }, &interval_algebra::Atan2, &affine_algebra::Atan2},
    };

    std::mt19937 rng(20261024);
    auto pick = [&]() {
        double a = kBounds[rng() % kBounds.size()], b = kBounds[rng() % kBounds.size()];
        if (!orderedTotal(a, b)) std::swap(a, b);
        return std::make_pair(a, b);
    };
    long tried = 0;
    std::vector<std::string>    failures;
    std::map<std::string, long> perOp;
    auto fail = [&](const std::string& s) {
        const std::string op = s.substr(0, s.find(' '));
        if (perOp[op]++ < 3) failures.push_back(s);
    };
    // an operand is a value of the library, flushed (its subnormals may be zero) ; the program
    // runs in one of three modes : gradual underflow, FTZ (results flushed), FTZ and DAZ
    // (operands and results flushed)
    for (int it = 0; it < 20000; it++) {
        const auto [xl, xh] = pick();
        const interval x    = flushedToZero(interval(xl, xh, -24));
        for (const Unary& u : unary) {
            for (int alg = 0; alg < 2; alg++) {
            const interval r = alg == 0 ? (A.*(u.op))(x) : toItv((F.*(u.aop))(fromItv(x)), 10);
            for (double v : samples(xl, xh)) {
                for (double y : {u.f(v), flushed(u.f(v)), flushed(u.f(flushed(v)))}) {
                    tried++;
                    if (!holds(r, y)) {
                        char buf[256];
                        snprintf(buf, sizeof(buf), "%s%s [%a, %a] at %a : %a outside [%a, %a]%s", u.name, alg ? "~" : "", x.lo(),
                                 x.hi(), v, y, r.lo(), r.hi(), r.mayBeInvalid() ? " +inv" : "");
                        fail(buf);
                    }
                }
            }
            }
        }
        const auto [yl, yh] = pick();
        const interval z    = flushedToZero(interval(yl, yh, -24));
        for (const Binary& b : binary) {
            for (int alg = 0; alg < 2; alg++) {
            const interval r = alg == 0 ? (A.*(b.op))(x, z) : toItv((F.*(b.aop))(fromItv(x), fromItv(z)), 10);
            for (double v : samples(xl, xh)) {
                for (double w : samples(yl, yh)) {
                    for (double y : {b.f(v, w), flushed(b.f(v, w)), flushed(b.f(flushed(v), flushed(w)))}) {
                        tried++;
                        if (!holds(r, y)) {
                            char buf[256];
                            snprintf(buf, sizeof(buf), "%s%s [%a, %a] [%a, %a] at %a, %a : %a outside [%a, %a]%s", b.name,
                                     alg ? "~" : "", x.lo(), x.hi(), z.lo(), z.hi(), v, w, y, r.lo(), r.hi(),
                                     r.mayBeInvalid() ? " +inv" : "");
                            fail(buf);
                        }
                    }
                }
            }
            }
        }
    }
    for (const std::string& s : failures) std::printf("%s\n", s.c_str());
    for (const auto& [op, n] : perOp) std::printf("%-6s %ld values outside\n", op.c_str(), n);
    std::printf("%ld values checked\n", tried);
    check("signed zeros : every value of an operation lies in its interval", true, failures.empty());
    return reportCheckResults();
}
