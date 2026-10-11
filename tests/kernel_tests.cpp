// Copyright 2026 Yann Orlarey, Stéphane Letz. Apache-2.0 (see LICENSE).
//
// The kernels Sum, Fir and Iir : the certified gain of a recursive filter against
// references computed independently, the refusal of unstable filters, and the inclusion
// of the values a program computes (float and double, random and adversarial inputs).

#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "interval/interval_algebra.hh"

using namespace itv;

namespace {

int gFailed = 0;

void expect(bool ok, const std::string& what)
{
    std::cout << (ok ? "OK   : " : "FAIL : ") << what << '\n';
    gFailed += !ok;
}

/// the L1 norm of the impulse response, summed until the terms vanish (Kahan, double)
double referenceGain(const std::vector<double>& c, long n = 20000000)
{
    std::vector<double> h(c.size() + 1, 0.0);  // h[0] newest
    double              sum = 0, comp = 0;
    for (long k = 0; k < n; k++) {
        double v = (k == 0) ? 1.0 : 0.0;
        for (size_t j = 0; j < c.size(); j++) v += c[j] * h[j];
        for (size_t j = c.size(); j > 0; j--) h[j] = h[j - 1];
        h[0]      = v;
        double y  = std::abs(v) - comp;
        double t  = sum + y;
        comp      = (t - sum) - y;
        sum       = t;
        if (k > 1000 && std::abs(v) < 1e-300) break;
    }
    return sum;
}

void gains()
{
    struct Case {
        std::string         name;
        std::vector<double> c;
    };
    const double pi = std::acos(-1.0);
    const Case   cases[] = {
        {"one pole 0.5", {0.5}},
        {"one pole -0.9", {-0.9}},
        {"one pole 0.999", {0.999}},
        {"pink noise filter", {2.494956002, -2.017265875, 0.522189400}},
        {"resonator r=0.999", {2 * 0.999 * std::cos(0.05 * pi), -0.999 * 0.999}},
    };
    for (const Case& k : cases) {
        FilterGain g;
        const bool ok  = filterGain(k.c, g);
        const double r = referenceGain(k.c);
        expect(ok && g.L >= r * (1 - 1e-9) && g.L <= r * 1.01,
               "gain of " + k.name + " : " + std::to_string(ok ? g.L : -1) + " for " + std::to_string(r));
    }
    // a slow resonator may not be certified (the squarings of a slowly rotating matrix
    // inflate their error bound) : the rule then claims nothing, never too little
    {
        const std::vector<double> c{2 * 0.9999 * std::cos(0.001 * pi), -0.9999 * 0.9999};
        FilterGain                g;
        const bool                ok = filterGain(c, g);
        expect(!ok || g.L >= referenceGain(c) * (1 - 1e-9), "a slow resonator : no bound, or a bound above its gain");
    }
    FilterGain g;
    expect(!filterGain({1.0001}, g), "an unstable one pole has no gain");
    expect(!filterGain({1.0}, g), "an integrator has no gain");
    expect(!filterGain({-1.0}, g), "a marginal pole has no gain");
    expect(!filterGain({2.0, -1.0}, g), "a double pole at 1 has no gain");
}

/// the program's y(n) = x(n) + sum c[j] y(n-1-j) in T, over a given input
template <typename T>
std::vector<double> run(const std::vector<double>& c, const std::vector<double>& x)
{
    std::vector<T>      y(c.size(), T(0));
    std::vector<double> out;
    for (double xv : x) {
        T v = T(xv);
        for (size_t j = 0; j < c.size(); j++) v = v + T(c[j]) * y[j];
        for (size_t j = c.size(); j-- > 1;) y[j] = y[j - 1];
        if (!c.empty()) y[0] = v;
        out.push_back(double(v));
    }
    return out;
}

template <typename T>
void inclusion(int precision)
{
    programPrecision() = precision;
    const std::string      p = precision == 1 ? " (float)" : " (double)";
    const interval_algebra A;
    std::mt19937           rng(99);
    const double           pi = std::acos(-1.0);
    const std::vector<std::vector<double>> filters = {
        {0.99}, {2.494956002, -2.017265875, 0.522189400}, {2 * 0.999 * std::cos(0.05 * pi), -0.999 * 0.999}, {-0.9}};
    int bad = 0;
    for (const auto& cd : filters) {
        std::vector<interval> cs;
        std::vector<double>   cp;  // the coefficients the program uses
        for (double v : cd) {
            cs.push_back(interval(v, v, -24));
            cp.push_back(precision == 1 ? double(float(v)) : v);
        }
        for (const interval X : {interval(-1, 1, -24), interval(0, 1, -24), interval(0.25, 0.75, -24)}) {
            const interval R = A.Iir(X, cs);
            if (!R.isValid() || !R.isBounded()) {
                bad++;
                continue;
            }
            // random input
            std::uniform_real_distribution<double> u(X.lo(), X.hi());
            std::vector<double>                    x(100000);
            for (double& v : x) v = double(T(u(rng)));
            for (double y : run<T>(cp, x)) bad += !(R.lo() <= y && y <= R.hi());
            // adversarial input : the signs of h, reversed, to reach the extremes at the end
            std::vector<double> h = run<double>(cp, std::vector<double>(1, 1.0));
            {
                std::vector<double> imp(200000, 0.0);
                imp[0] = 1;
                h      = run<double>(cp, imp);
            }
            for (int s : {1, -1}) {
                std::vector<double> xa(h.size());
                for (size_t k = 0; k < h.size(); k++) {
                    const double hk = h[h.size() - 1 - k] * s;
                    xa[k]           = hk >= 0 ? X.hi() : X.lo();
                }
                for (double y : run<T>(cp, xa)) bad += !(R.lo() <= y && y <= R.hi());
            }
        }
    }
    expect(bad == 0, "inclusion of iir, random and adversarial inputs" + p +
                         (bad ? " -- " + std::to_string(bad) + " values outside" : ""));

    // the rounding of the program can go beyond the exact bound : a resonator whose float
    // recursion rises 171 ulps above its exact peak (found by a search)
    if (precision == 1) {
        const float    c1 = 0x1.efea18p+0f, c2 = -0x1.e08p-1f, xf = 0x1.b69c92p-1f;
        const interval R  = A.Iir(interval(double(xf), double(xf), -24),
                                  {interval(double(c1), double(c1), -24), interval(double(c2), double(c2), -24)});
        float          y1 = 0, y2 = 0;
        bool           in = true;
        for (long n = 0; n < 300000; n++) {
            const float y = xf + c1 * y1 + c2 * y2;
            y2            = y1;
            y1            = y;
            in            = in && R.lo() <= double(y) && double(y) <= R.hi();
        }
        expect(in, "the program's rounding beyond the exact gain is included (float)");
    }

    // a smoother of a nonnegative input stays nonnegative, below its gain
    const interval S = A.Iir(interval(0, 1, -24), {interval(0.99, 0.99, -24)});
    expect(S.lo() > -1e-6 && S.hi() < 100.01 && S.isValid(), "a smoother of [0, 1] by 0.99 is [0, 100]" + p);
    // no bound without constant coefficients, or for an unstable filter : anything
    const interval V = A.Iir(interval(0, 1, -24), {interval(0.5, 0.9, -24)});
    expect(!V.isValid() && !V.isBounded(), "a variable coefficient gives anything" + p);
    const interval U = A.Iir(interval(0, 1, -24), {interval(1, 1, -24)});
    expect(!U.isValid() && !U.isBounded(), "an integrator gives anything" + p);
    expect(!A.Iir(interval(0, 1, -24).withInvalid(true), {interval(0.5, 0.5, -24)}).isValid(), "iir keeps an alert" + p);
    // an integer recursion wraps : a congruential generator stays in int32, valid
    const interval G = A.Iir(interval(12345, 12345, 0), {interval(1103515245, 1103515245, 0)});
    expect(G.isValid() && G.lo() >= -2147483648.0 && G.hi() <= 2147483647.0, "an integer recursion stays in int32" + p);
}

void firs()
{
    programPrecision() = 2;
    const interval_algebra A;
    const interval         F = A.Fir(interval(0, 1, -24), {interval(1, 1, 0), interval(0, 0, 0), interval(-2, -2, 0)});
    expect(F.lo() <= -2 && F.lo() > -2.001 && F.hi() >= 1 && F.hi() < 1.001, "fir of [0, 1] by (1, 0, -2) is [-2, 1]");
    const interval Z = A.Fir(interval(0, 1, -24), {interval(0, 0, 0)});
    expect(Z.lo() == 0 && Z.hi() == 0, "an all-zero fir is 0");
    const interval S = A.Sum({interval(0, 1, -24), interval(2, 3, -24), interval(-1, 0, -24)});
    expect(S.lo() <= 1 && S.lo() > 0.999 && S.hi() >= 4 && S.hi() < 4.001, "sum of [0,1], [2,3], [-1,0] is [1, 4]");
    // a product from left to right : the bounds of the binary chain, bit for bit
    programPrecision() = 1;
    const interval a(0.1, 0.3, -24), b(-3.7, -2.1, -24), c(1.0 / 3, 0.7, -24);
    const interval P = A.Prod({a, b, c});
    const interval Q = A.Mul(A.Mul(a, b), c);
    expect(P.lo() == Q.lo() && P.hi() == Q.hi(), "prod(a, b, c) is (a * b) * c");
    expect(P.lo() <= -0.3 * 3.7 * 0.7 && P.hi() >= -0.1 * 2.1 / 3, "prod of [0.1,0.3], [-3.7,-2.1], [1/3,0.7] holds its exact range");
    const interval U = A.Prod({});
    expect(U.lo() == 1 && U.hi() == 1, "an empty product is 1");
    // a comb x + g y(n-N) is bounded as the first-order filter x + g y(n-1), whatever N
    programPrecision() = 2;
    std::vector<interval> comb(2000, interval(0, 0, 0));
    comb.back()                 = interval(0.7, 0.7, -52);
    const interval C1           = A.Iir(interval(-1, 1, -52), {interval(0.7, 0.7, -52)});
    const interval CN           = A.Iir(interval(-1, 1, -52), comb);
    expect(CN.lo() == C1.lo() && CN.hi() == C1.hi() && CN.isBounded(), "a comb of 2000 samples has the bound of its first-order filter");
    std::vector<interval> two(6, interval(0, 0, 0));
    two[2] = interval(0.5, 0.5, -52);
    two[5] = interval(-0.25, -0.25, -52);
    const interval T2 = A.Iir(interval(-1, 1, -52), two);
    const interval T1 = A.Iir(interval(-1, 1, -52), {interval(0.5, 0.5, -52), interval(-0.25, -0.25, -52)});
    expect(T2.lo() == T1.lo() && T2.hi() == T1.hi(), "taps at 3 and 6 : the filter of taps 1 and 2");
}

}  // namespace

int main()
{
    gains();
    inclusion<float>(1);
    inclusion<double>(2);
    firs();
    std::cout << (gFailed ? "FAILED" : "PASSED") << " (" << gFailed << " failure(s))\n";
    return gFailed ? 1 : 0;
}
