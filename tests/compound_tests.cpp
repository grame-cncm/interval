// Copyright 2026 Yann Orlarey, Stéphane Letz. Apache-2.0 (see LICENSE).
//
// The compound operations Frac, Cfrac, Mix, Tpt and the recursive mix Irm : precision of the bounds, inclusion of
// the values a program computes (in float and in double, with and without a fused
// multiply-add of the product that computes the operand), and validity.

#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <random>
#include <string>

#include "interval/affine_ops.hh"
#include "interval/interval_algebra.hh"

using namespace itv;

namespace {

int gFailed = 0;

void expect(bool ok, const std::string& what)
{
    std::cout << (ok ? "OK   : " : "FAIL : ") << what << '\n';
    gFailed += !ok;
}

bool within(double v, const interval& I)
{
    return std::isnan(v) ? I.mayBeInvalid() : (I.lo() <= v && v <= I.hi());
}

/// a program value of x : rounded to the precision of the program
template <typename T>
T value(double x)
{
    return static_cast<T>(x);
}

/// every value the program computes for op over random operands of the interval X
/// lies in result ; "contracted" makes x a product a * b fused with the subtraction
template <typename T>
int inclusion(const interval& X, const interval& result, T (*round)(T), std::mt19937& rng)
{
    int bad = 0;
    std::uniform_real_distribution<double> u(X.lo(), X.hi());
    for (int k = 0; k < 20000; k++) {
        const T x = value<T>(k == 0 ? X.lo() : k == 1 ? X.hi() : u(rng));
        // plain : x computed once, then x - round(x)
        if (!within(double(x - round(x)), result)) bad++;
        // contracted : x = a * b, the subtraction fused with the product
        const T a = value<T>(1.0 + std::uniform_real_distribution<double>(0, 1)(rng));
        const T b = x / a;
        const T p = a * b;
        if (p < value<T>(X.lo()) || p > value<T>(X.hi())) continue;
        if (!within(double(std::fma(a, b, -round(p))), result)) bad++;
    }
    return bad;
}

template <typename T>
T floorT(T v)
{
    return std::floor(v);
}
template <typename T>
T rintT(T v)
{
    return std::rint(v);
}

template <typename T>
void remainders(int precision)
{
    programPrecision() = precision;
    const std::string      p = precision == 1 ? " (float)" : " (double)";
    const interval_algebra A;
    std::mt19937           rng(4321);

    // precision
    const interval f1 = A.Frac(interval(0.2, 0.7, -24));
    expect(f1.lo() > 0.199 && f1.hi() < 0.701 && f1.isValid(), "frac of [0.2, 0.7] is [0.2, 0.7]" + p);
    const interval f2 = A.Frac(interval(0, 10, -24));
    expect(f2.lo() >= -1e-5 && f2.hi() < 1 + 1e-5 && f2.isValid(), "frac of [0, 10] is [0, 1[" + p);
    const interval f3 = A.Frac(interval(3.25, 3.75, -24));
    expect(f3.lo() > 0.249 && f3.hi() < 0.751, "frac of [3.25, 3.75] is [0.25, 0.75]" + p);
    const interval f4 = A.Frac(interval(3, 7, 0));
    expect(f4.lo() == 0 && f4.hi() == 0, "frac of an integer is 0" + p);
    const interval r1 = A.Cfrac(interval(-0.3, 0.3, -24));
    expect(r1.lo() > -0.301 && r1.hi() < 0.301, "cfrac of [-0.3, 0.3] is itself" + p);
    const interval r2 = A.Cfrac(interval(-40, 40, -24));
    expect(r2.lo() >= -0.5 - 1e-4 && r2.hi() <= 0.5 + 1e-4, "cfrac of [-40, 40] is [-1/2, 1/2]" + p);

    // validity
    expect(!A.Frac(interval(-HUGE_VAL, HUGE_VAL, -24)).isValid(), "frac of an infinite x may be NaN" + p);
    expect(!A.Frac(interval(0, 1, -24).withInvalid(true)).isValid(), "frac keeps an alert" + p);
    expect(!A.Cfrac(interval(0, HUGE_VAL, -24)).isValid(), "cfrac of an infinite x may be NaN" + p);

    // inclusion, plain and contracted, over operands of every kind
    const interval cases[] = {interval(0, 1, -24),       interval(0, 10, -24),     interval(-10, 10, -24),
                              interval(3.2, 3.9, -24),   interval(-3.9, -3.2, -24), interval(0.999, 1.001, -24),
                              interval(1e6, 1.0001e6, -24), interval(0, 536870912, -24), interval(-0.5, 0.5, -24)};
    int bad = 0;
    for (const interval& X : cases) {
        bad += inclusion<T>(X, A.Frac(X), floorT<T>, rng);
        bad += inclusion<T>(X, A.Cfrac(X), rintT<T>, rng);
    }
    expect(bad == 0, "inclusion of frac and cfrac, plain and contracted" + p +
                         (bad ? " -- " + std::to_string(bad) + " values outside" : ""));
}

template <typename T>
void mixes(int precision)
{
    programPrecision() = precision;
    const std::string      p = precision == 1 ? " (float)" : " (double)";
    const interval_algebra A;
    std::mt19937           rng(1234);

    const interval m1 = A.Mix(interval(0.1, 0.999, -24), interval(-12, 12, -24), interval(-12, 12, -24));
    expect(m1.lo() > -12.01 && m1.hi() < 12.01 && m1.isValid(), "a mix of [-12, 12] stays in [-12, 12]" + p);
    const interval m2 = A.Mix(interval(0, 1, -24), interval(0, 5, -24), interval(2, 3, -24));
    expect(m2.lo() >= 0 && m2.hi() < 5.01, "a mix of nonnegative operands stays nonnegative" + p);
    const interval m3 = A.Mix(interval(-1, 2, -24), interval(0, 1, -24), interval(0, 1, -24));
    expect(m3.lo() < -0.5 && m3.hi() > 1.5, "beyond c in [0, 1], the mix is its definition" + p);
    expect(!A.Mix(interval(0, 1, -24), interval(0, 1, -24).withInvalid(true), interval(0, 1, -24)).isValid(),
           "a mix keeps an alert" + p);

    // the program's mix : an interpolation from the operand whose weight may be near 1
    auto program = [](const interval& ic, T c, T a, T b, bool fused) {
        if (mixFromB(ic)) {
            const T w = T(1) - c;
            return fused ? T(std::fma(w, a - b, b)) : T(b + w * (a - b));
        }
        return fused ? T(std::fma(c, b - a, a)) : T(a + c * (b - a));
    };

    // a case found by search : a + c (b - a) is nearly 4 rounding errors of the largest
    // operand away from the exact mix, in float ; the program interpolates from b
    if (precision == 1) {
        const float    a = 0.5153635740280151f, b = -0.5094743371009827f, c = 0.9999904036521912f;
        const interval ic(c, c, -24);
        const interval R = A.Mix(ic, interval(a, a, -24), interval(b, b, -24));
        expect(within(double(program(ic, c, a, b, false)), R), "the interpolation of the witness is inside the mix" + p);
    }

    // a slow smoother : a weight 1 - c of a few rounding errors keeps a margin of about one
    // rounding of the result (a margin of 5 would exceed the stability margin of the recursion)
    {
        const double   u  = precision == 1 ? std::ldexp(1.0, -24) : std::ldexp(1.0, -53);
        const double   c  = 1 - 8 * u;
        const interval R  = A.Mix(interval(c, c, -24), interval(-1, 1, -24), interval(-1, 1, -24));
        // (plus the 2^-50 relative slack of the ends, computed in double)
        const double   m  = 2.5 * u + std::ldexp(1.0, -50);
        expect(R.hi() <= 1 + m && R.lo() >= -1 - m, "near c = 1, the margin of the mix is one rounding" + p);
    }

    // the power-of-two envelope : both forms, plain and fused, over operands of [-M, M]
    // (the extremes included) and weights of [0, 1] (near 0, 1/2 and 1 included), stays in
    // [-M, M] ; a slow smoother then maps [-M, M] into itself
    {
        int out = 0;
        for (double M : {1.0, 16.0, 1024.0, std::ldexp(1.0, -60)}) {
            std::uniform_real_distribution<double> uv(-M, M), uw(0, 1);
            const double                           u = precision == 1 ? std::ldexp(1.0, -24) : std::ldexp(1.0, -53);
            const double specials[] = {0, u, 0.5 - u, 0.5, 0.5 + u, 1 - 4 * u, 1 - 2 * u, 1 - u, 1};
            for (int s = 0; s < 200000; s++) {
                const T c = value<T>(s % 3 ? uw(rng) : specials[s % 9]);
                const T a = value<T>(s % 5 == 0 ? M : s % 5 == 1 ? -M : uv(rng));
                const T b = value<T>(s % 7 == 0 ? M : s % 7 == 1 ? -M : uv(rng));
                for (const interval& ic : {interval(0, 1, -24), interval(0.5, 1, -24)}) {
                    for (bool fused : {false, true}) {
                        if (std::abs(double(program(ic, c, a, b, fused))) > M) out++;
                    }
                }
            }
        }
        expect(out == 0, "both forms of mix stay in a power-of-two envelope" + p +
                             (out ? " -- " + std::to_string(out) + " values outside" : ""));
        const double   u = precision == 1 ? std::ldexp(1.0, -24) : std::ldexp(1.0, -53);
        const interval R = A.Mix(interval(0.4, 1 - 2 * u, -24), interval(-1, 1, -24), interval(-1, 1, -24));
        expect(R.lo() >= -1 && R.hi() <= 1, "a smoother maps [-1, 1] into itself" + p);
        const interval Q = A.Mix(interval(0, 1, -24), interval(-12, 12, -24), interval(-12, 12, -24));
        expect(Q.lo() >= -16 && Q.hi() <= 16, "a mix of [-12, 12] stays in [-16, 16]" + p);
    }

    // inclusion : the program's interpolation, plain and with the product fused
    struct Case {
        interval c, a, b;
    };
    const Case cases[] = {{interval(0, 1, -24), interval(-12, 12, -24), interval(-12, 12, -24)},
                          {interval(0.999, 0.9999999, -24), interval(0, 1, -24), interval(0, 1, -24)},
                          {interval(0, 0.5, -24), interval(100, 1000, -24), interval(-3, 3, -24)},
                          {interval(0.2, 0.3, -24), interval(-1e6, 1e6, -24), interval(5, 6, -24)},
                          {interval(0.5, 1, -24), interval(-12, 12, -24), interval(-12, 12, -24)},
                          {interval(0, 0.9999999, -24), interval(-1, 1, -24), interval(-1, 1, -24)},
                          {interval(0.9999, 0.99999994, -24), interval(-1, 1, -24), interval(-1, 1, -24)}};
    int bad = 0;
    for (const Case& k : cases) {
        const interval R = A.Mix(k.c, k.a, k.b);
        std::uniform_real_distribution<double> uc(k.c.lo(), k.c.hi()), ua(k.a.lo(), k.a.hi()),
            ub(k.b.lo(), k.b.hi());
        for (int s = 0; s < 20000; s++) {
            const T c = value<T>(s % 2 ? k.c.hi() : uc(rng)), a = value<T>(ua(rng)), b = value<T>(ub(rng));
            if (!within(double(program(k.c, c, a, b, false)), R)) bad++;
            if (!within(double(program(k.c, c, a, b, true)), R)) bad++;
        }
    }
    expect(bad == 0, "inclusion of mix, plain and fused" + p +
                         (bad ? " -- " + std::to_string(bad) + " values outside" : ""));
}

template <typename T>
void tpts(int precision)
{
    programPrecision() = precision;
    const std::string      p = precision == 1 ? " (float)" : " (double)";
    const interval_algebra A;
    std::mt19937           rng(777);

    const interval t1 = A.Tpt(interval(0, 1000, -24));
    expect(t1.lo() >= -1e-12 && t1.lo() <= 0 && t1.hi() < 1 && t1.isValid(), "tpt of [0, 1000] lies in [0, 1[" + p);
    const interval t2 = A.Tpt(interval(0.5, 0.5, -24));
    expect(t2.lo() > 0.3333 && t2.hi() < 0.3334, "tpt of 1/2 is 1/3" + p);
    const interval t3 = A.Tpt(interval(-0.5, 3, -24));
    expect(t3.lo() > -1.0001 && t3.hi() < 0.7501, "tpt of [-1/2, 3] is [-1, 3/4]" + p);
    const interval t4 = A.Tpt(interval(-2, 2, -24));
    expect(!t4.isBounded() || !t4.isValid(), "tpt of [-2, 2] reaches the pole : the definition" + p);
    expect(!A.Tpt(interval(0, HUGE_VAL, -24)).isValid(), "tpt of an infinite x may be NaN" + p);
    expect(!A.Tpt(interval(0, 1, -24).withInvalid(true)).isValid(), "tpt keeps an alert" + p);

    // inclusion : the program's x / (x + 1)
    const interval cases[] = {interval(0, 1, -24), interval(0, 1e6, -24), interval(-0.5, 2, -24),
                              interval(1e-30, 1e-20, -24), interval(0, 1.6e16, -24),
                              interval(-0.999999, -0.99, -24)};
    int bad = 0;
    for (const interval& X : cases) {
        const interval R = A.Tpt(X);
        std::uniform_real_distribution<double> u(X.lo(), X.hi());
        for (int k = 0; k < 20000; k++) {
            const T x = value<T>(k == 0 ? X.lo() : k == 1 ? X.hi() : u(rng));
            if (double(x) < X.lo() || double(x) > X.hi()) continue;
            if (!within(double(x / (x + T(1))), R)) bad++;
        }
    }
    expect(bad == 0, "inclusion of tpt" + p + (bad ? " -- " + std::to_string(bad) + " values outside" : ""));

    // the ends themselves, over many intervals : the rounding of the program decides
    // whether an end value lands above or below the exact image, half of the time
    int ends = 0;
    std::uniform_real_distribution<double> e(-0.9, 100);
    for (int k = 0; k < 2000; k++) {
        const T a = value<T>(e(rng)), b = value<T>(e(rng));
        const interval X(std::min(double(a), double(b)), std::max(double(a), double(b)), -24);
        const interval R = A.Tpt(X);
        for (T x : {T(X.lo()), T(X.hi())}) {
            if (!within(double(x / (x + T(1))), R)) {
                if (ends++ < 3) {
                    std::printf("  x=%.9g X=[%.9g, %.9g] value=%.17g R=[%.17g, %.17g] def=%s\n", double(x),
                                X.lo(), X.hi(), double(x / (x + T(1))), R.lo(), R.hi(),
                                A.Div(X, A.Add(X, A.FloatNum(1.0))).isBounded() ? "bounded" : "unbounded");
                }
            }
        }
    }
    expect(ends == 0, "inclusion of tpt at the ends" + p + (ends ? " -- " + std::to_string(ends) + " values outside" : ""));
}

template <typename T>
void phases(int precision)
{
    programPrecision() = precision;
    const std::string      p = precision == 1 ? " (float)" : " (double)";
    const interval_algebra A;
    std::mt19937           rng(4242);

    // the counter of ba.period(P) / P : x % P below P, P of a variable rate
    const interval c = A.Phase(interval(1, 19200, 0), interval(5, 19200, 0));
    expect(c.lo() >= 0 && c.hi() <= 1 && c.isValid(), "phase of a counter lies in [0, 1]" + p);
    const interval d = A.Div(A.Mod(interval(1, 19200, 0), interval(5, 19200, 0)), interval(5, 19200, 0));
    expect(d.hi() > 1000, "its definition, read as two values, does not" + p);
    const interval n = A.Phase(interval(-10, -1, 0), interval(2, 100, 0));
    expect(n.lo() >= -1 && n.hi() <= 0, "phase of a negative x lies in [-1, 0]" + p);
    const interval m = A.Phase(interval(-10, 10, 0), interval(-100, -2, 0));
    expect(m.lo() >= -1 && m.hi() <= 1, "phase of any x lies in [-1, 1]" + p);
    const interval q = A.Phase(interval(0, 10, 0), interval(-100, -2, 0));
    expect(q.lo() >= -1 && q.hi() <= 0, "phase by a negative p lies in [-1, 0]" + p);
    expect(!A.Phase(interval(0, 10, 0), interval(0, 100, 0)).isValid(), "a divisor that may be 0 alerts" + p);
    expect(!A.Phase(interval(0, 10, 0), interval(2, 4, 0).withInvalid(true)).isValid(), "phase keeps an alert" + p);

    // inclusion : the program's float(x % p) / float(p) on integers, fmod(x, p) / p on reals
    int bad = 0;
    std::uniform_int_distribution<int> ux(0, 1 << 30), up(1, 1 << 30);
    for (int k = 0; k < 20000; k++) {
        const int      xi = ux(rng), pi = k < 100 ? xi + 1 : up(rng);
        const interval R  = A.Phase(interval(0, 1 << 30, 0), interval(1, 1 << 30, 0));
        if (!within(double(T(xi % pi) / T(pi)), R)) bad++;
    }
    std::uniform_real_distribution<double> rx(0, 1000), rp(0.5, 100);
    const interval Rr = A.Phase(interval(0, 1000, -24), interval(0.5, 100, -24));
    for (int k = 0; k < 20000; k++) {
        const T x = value<T>(rx(rng)), q2 = value<T>(rp(rng));
        if (!within(double(T(std::fmod(x, q2)) / q2), Rr)) bad++;
    }
    expect(bad == 0, "inclusion of phase" + p + (bad ? " -- " + std::to_string(bad) + " values outside" : ""));
}

// Irm(x, c) : y(n) = mix(c, x(n), y(n - 1)), y(-1) = 0, the program computing the mix as
// an interpolation (mixFromB on the interval of c), plain or with the product fused
template <typename T>
void irms(int precision)
{
    programPrecision() = precision;
    const std::string      p = precision == 1 ? " (float)" : " (double)";
    const interval_algebra A;
    std::mt19937           rng(4321);
    const double           u = precision == 1 ? std::ldexp(1.0, -24) : std::ldexp(1.0, -53);

    const interval s1 = A.Irm(interval(0, 384000, -24), interval(0.99, 0.9999, -24));
    expect(s1.isValid() && s1.lo() == 0 && s1.hi() >= 384000 && s1.hi() < 384000 * 1.01,
           "a slow smoother of [0, 384000] stays within 1 % of its target" + p);
    const interval s2 = A.Irm(interval(-1, 1, -24), interval(0, 1, -24));
    expect(s2.isValid() && s2.lo() >= -1 && s2.hi() <= 1, "a smoother of [-1, 1] stays in [-1, 1], its pole up to 1" + p);
    const interval s3 = A.Irm(interval(-3, -2, -24), interval(0.5, 0.9, -24));
    expect(s3.isValid() && s3.hi() < 0 && s3.lo() <= -3 && s3.lo() > -3.1, "a smoother of [-3, -2] stays in [-3, 0[" + p);
    const interval s6 = A.Irm(interval(1, 20, -24), interval(0.9, 0.99, -24));
    expect(s6.isValid() && s6.lo() > 0 && s6.hi() < 20.01, "a smoother of [1, 20] is positive : its samples, not its initial state" + p);
    const interval s4 = A.Irm(interval(0, 1, -24), interval(0, 1.5, -24));
    expect(!s4.isValid() && !std::isfinite(s4.hi()), "a pole beyond 1 : anything" + p);
    const interval s5 = A.Irm(interval(0, HUGE_VAL, -24), interval(0.5, 0.9, -24));
    expect(!s5.isValid(), "an unbounded input : anything" + p);
    expect(!A.Irm(interval(0, 1, -24).withInvalid(true), interval(0.5, 0.9, -24)).isValid(), "Irm keeps an alert" + p);

    // the two forms of the interpolation, whichever a code generator chooses : from the
    // state with the weight 1 - c of x (computed, or exact when the program holds it), from
    // x with the pole c
    auto program = [](int form, T c, T a, T b, bool fused) {
        if (form == 0) {
            const T w = T(1) - c;
            return fused ? T(std::fma(w, a - b, b)) : T(b + w * (a - b));
        }
        return fused ? T(std::fma(c, b - a, a)) : T(a + c * (b - a));
    };
    struct Case {
        interval x, c;
    };
    const Case cases[] = {{interval(0, 384000, -24), interval(0.99, 0.9999, -24)},
                          {interval(-1, 1, -24), interval(0, 1, -24)},
                          {interval(-1, 1, -24), interval(1 - 64 * u, 1 - u, -24)},
                          {interval(5, 6, -24), interval(0, 0.4, -24)},
                          {interval(-1e6, 3, -24), interval(0.3, 0.7, -24)},
                          {interval(48000, 48000, -24), interval(0.9999, 0.99999, -24)}};
    int bad = 0;
    for (const Case& k : cases) {
        const interval                         R = A.Irm(k.x, k.c);
        std::uniform_real_distribution<double> ux(k.x.lo(), k.x.hi()), uc(k.c.lo(), k.c.hi());
        for (int form : {0, 1}) {
            for (bool fused : {false, true}) {
                T y = 0;
                T x = value<T>(k.x.hi());
                for (int n = 0; n < 100000; n++) {
                    // held targets, extremes first, then random ones
                    if (n % 5000 == 0) x = value<T>(n == 0 ? k.x.hi() : n == 5000 ? k.x.lo() : ux(rng));
                    const T c = value<T>(n % 3 ? uc(rng) : k.c.hi());
                    y         = program(form, c, x, y, fused);
                    if (!within(double(y), R)) bad++;
                }
            }
        }
    }
    expect(bad == 0, "inclusion of the recursive mix, plain and fused" + p +
                         (bad ? " -- " + std::to_string(bad) + " values outside" : ""));
}

// Bpf(x, t, r) : the regions that x reaches, the last one for a NaN x
void bpfs()
{
    programPrecision() = 1;
    const interval_algebra      A;
    const std::vector<interval> t{interval(10, 10, -24), interval(20, 20, -24)};
    const std::vector<interval> r{interval(0, 0, -24), interval(1, 5, -24), interval(-3, -2, -24)};
    const interval              b1 = A.Bpf(interval(12, 15, -24), t, r);
    expect(b1.isValid() && b1.lo() == 1 && b1.hi() == 5, "bpf : x within one region reads that region only");
    const interval b2 = A.Bpf(interval(5, 20, -24), t, r);
    expect(b2.lo() == -3 && b2.hi() == 5, "bpf : x reaching a threshold reads the region above it");
    const interval b3 = A.Bpf(interval(5, 19.5, -24), t, r);
    expect(b3.lo() == 0 && b3.hi() == 5, "bpf : x below the last threshold does not read the last region");
    const interval b4 = A.Bpf(interval(12, 15, -24).withInvalid(true), t, r);
    expect(b4.lo() == -3 && b4.hi() == 5, "bpf : a NaN x reads the last region");
    std::vector<interval> r2 = r;
    r2[2]                    = r2[2].withInvalid(true);
    expect(A.Bpf(interval(12, 15, -24), t, r2).isValid(), "bpf : an unreached invalid region leaves the result valid");
    expect(!A.Bpf(interval(12, 25, -24), t, r2).isValid(), "bpf : a reached invalid region makes it invalid");
}

void affine()
{
    programPrecision()        = 2;
    const affine_algebra A(1024);
    // a growing corridor [0, t] : frac collapses over the horizon to [0, 1[
    AffItv x = A.FloatNum(0);
    x.a0 = 0, x.a1 = 0, x.b0 = 0, x.b1 = 1, x.lsb = -24;
    const interval f = toItv(A.Frac(x), 1024);
    expect(f.lo() >= -1e-9 && f.hi() < 1 + 1e-9, "affine frac of a growing corridor is [0, 1[");
    const interval m = toItv(A.Mix(A.FloatNum(0.5), A.FloatNum(2), A.FloatNum(4)), 1024);
    expect(m.lo() > 2.99 && m.hi() < 3.01, "affine mix of constants");
}

}  // namespace

int main()
{
    phases<float>(1);
    phases<double>(2);
    remainders<float>(1);
    remainders<double>(2);
    mixes<float>(1);
    mixes<double>(2);
    tpts<float>(1);
    tpts<double>(2);
    irms<float>(1);
    irms<double>(2);
    bpfs();
    affine();
    std::cout << (gFailed ? "FAILED" : "PASSED") << " (" << gFailed << " failure(s))\n";
    return gFailed ? 1 : 0;
}
