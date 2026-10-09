// Copyright 2026 Yann Orlarey, Stéphane Letz. Apache-2.0 (see LICENSE).
#include <array>
#include <cmath>
#include <string>

#include "interval/affine_ops.hh"
#include "interval/check.hh"

using namespace itv;

namespace {
// Equality of forms checks the identity itself, not merely its horizon hull.
// NaNs here are numeric-bottom sentinels and must compare like other coefficients.
bool sameCoefficient(double x, double y)
{
    return x == y || (std::isnan(x) && std::isnan(y));
}
bool sameForm(const AffItv& x, const AffItv& y)
{
    return sameCoefficient(x.a0, y.a0) && sameCoefficient(x.a1, y.a1) &&
           sameCoefficient(x.b0, y.b0) && sameCoefficient(x.b1, y.b1) &&
           x.lsb == y.lsb && x.mayBeInvalid == y.mayBeInvalid;
}

void identityTests(int precision)
{
    programPrecision() = precision;
    const std::string p = precision == 1 ? "float delay: " : "double delay: ";
    const interval_algebra ordinary;
    const affine_algebra affine(32);
    const std::array<interval, 9> signals = {
        interval(2, 8, 0), interval(-8, -2, 0), interval(4, 16, 2),
        interval(0.5, 2, -24), interval(-2, -0.5, -24),
        interval(-1, 1, -12), interval(HUGE_VAL, HUGE_VAL, -24),
        empty(-24), empty(2)
    };
    // The delay operand is implicitly converted to int32. Its whole valid image
    // must be zero; a fractional corridor may meet this condition too.
    const std::array<interval, 6> amounts = {
        ordinary.IntNum(0), ordinary.FloatNum(0), ordinary.FloatNum(-0.0),
        ordinary.FloatNum(0.75), ordinary.FloatNum(-0.75), interval(-0.75, 0.75, -24)
    };
    for (size_t i = 0; i < signals.size(); ++i) {
        for (size_t j = 0; j < amounts.size(); ++j) {
            for (bool signalAlert : {false, true}) {
                for (bool amountAlert : {false, true}) {
                    const interval x = signals[i].withInvalid(signalAlert);
                    const interval n = amounts[j].withInvalid(amountAlert);
                    const interval expected = x.withInvalid(amountAlert);
                    const std::string name = p + "identity " + std::to_string(i) + "/" +
                        std::to_string(j) + "/" + std::to_string(signalAlert) +
                        std::to_string(amountAlert);
                    checkExact(name + " ordinary", ordinary.Delay(x, n), expected);
                    check(name + " affine", true,
                          sameForm(affine.Delay(fromItv(x), fromItv(n)), fromItv(expected)));
                }
            }
        }
    }

    // Positive/negative rates exercise both branches that normally flatten a
    // temporal bound. Infinite rates also catch accidental infinity times zero.
    const std::array<AffItv, 5> forms = {
        AffItv{2, 0.125, 4, 0.25, -24}, AffItv{-4, -0.25, -2, -0.125, -24},
        AffItv{8, 1, 16, 2, 2}, AffItv{2, -0.125, 4, 0.125, -24},
        AffItv{-1, -HUGE_VAL, 1, HUGE_VAL, -24}
    };
    for (size_t i = 0; i < forms.size(); ++i) {
        for (bool signalAlert : {false, true}) {
            for (bool amountAlert : {false, true}) {
                const AffItv x = forms[i].withInvalid(signalAlert);
                const AffItv n = affine.IntNum(0).withInvalid(amountAlert);
                check(p + "zero delay preserves every affine coefficient " + std::to_string(i),
                      true, sameForm(affine.Delay(x, n), x.withInvalid(amountAlert)));
            }
        }
    }
}

void initializationTests(int precision)
{
    programPrecision() = precision;
    const std::string p = precision == 1 ? "float initialization: " : "double initialization: ";
    const interval_algebra ordinary;
    const affine_algebra affine(32);
    for (const interval& x : {interval(2, 8, 2), interval(-8, -2, 2),
                              interval(0.5, 2, -24), interval(-2, -0.5, -24)}) {
        for (const interval& n : {ordinary.IntNum(1), interval(0, 3, 0),
                                  interval(1, 3, 0), interval(0.75, 1.75, -24)}) {
            const interval expected(std::min(0.0, x.lo()), std::max(0.0, x.hi()), x.lsb());
            checkExact(p + "positive or possibly positive ordinary delay includes initial zero",
                       ordinary.Delay(x, n), expected);
            checkExact(p + "positive or possibly positive affine delay includes initial zero",
                       toItv(affine.Delay(fromItv(x), fromItv(n)), 32), expected);
        }
        checkExact(p + "Mem retains one-sample initialization", toItv(affine.Mem(fromItv(x)), 32),
                   interval(std::min(0.0, x.lo()), std::max(0.0, x.hi()), x.lsb()));
    }
    const AffItv moving{2, 0.125, 4, 0.25, -24};
    check(p + "minimum zero does not make a variable delay the identity", true,
          toItv(affine.Delay(moving, fromItv(interval(0, 3, 0))), 32).hasZero());
    check(p + "positive delay still shifts the growing upper bound", true,
          affine.Delay(moving, affine.IntNum(1)).hi(32) < moving.hi(32));
    // First truncate 1.75 to one sample: sliding by 1.75 would underestimate
    // the upper bound of a growing signal.
    check(p + "fractional delay uses the converted integer sample count", true,
          sameForm(affine.Delay(moving, affine.FloatNum(1.75)),
                   affine.Delay(moving, affine.IntNum(1))));
    check(p + "negative delay retains its domain alert", true,
          ordinary.Delay(interval(2, 8), ordinary.IntNum(-1)).mayBeInvalid() &&
          affine.Delay(moving, affine.IntNum(-1)).mayBeInvalid);
    check(p + "invalid conversion cannot certify a zero delay", true,
          ordinary.Delay(interval(2, 8), ordinary.FloatNum(HUGE_VAL)).mayBeInvalid() &&
          affine.Delay(moving, affine.FloatNum(HUGE_VAL)).mayBeInvalid);
}

void reciprocalWitness(int precision)
{
    programPrecision() = precision;
    const std::string p = precision == 1 ? "float witness 03: " : "double witness 03: ";
    const interval_algebra ordinary;
    const affine_algebra affine(32);
    // Isolate the library portion of the smoothed-frequency witness. The compiler
    // supplies a strictly positive bound; reading it at @0 must not invent zero.
    const interval frequency(0x1p-10, 1000, -24);
    const interval delayed = ordinary.Delay(frequency, ordinary.IntNum(0));
    const AffItv adelayed = affine.Delay(fromItv(frequency), affine.IntNum(0));
    const interval inverse = ordinary.Inv(delayed);
    const interval ainverse = toItv(affine.Inv(adelayed), 32);
    checkExact(p + "ordinary frequency range is unchanged", delayed, frequency);
    checkExact(p + "affine frequency range is unchanged", toItv(adelayed, 32), frequency);
    check(p + "both reciprocals stay positive and finite", true,
          inverse.lo() > 0 && inverse.isBounded() && !inverse.mayBeInvalid() &&
          ainverse.lo() > 0 && ainverse.isBounded() && !ainverse.mayBeInvalid());
    check(p + "both integer conversions remain defined", true,
          ordinary.IntCast(inverse).isValid() &&
          toItv(affine.IntCast(affine.Inv(adelayed)), 32).isValid());
    for (double sample : {frequency.lo(), 1.0, 50.0, 440.0, frequency.hi()}) {
        // These exact powers of two and ordinary positive values remain finite;
        // a runtime division is an independent inclusion witness in each precision.
        const double result = precision == 1 ? double(1.0f / float(sample)) : 1.0 / sample;
        check(p + "runtime reciprocal is included", true, inverse.has(result) && ainverse.has(result));
    }
}
}  // namespace

int main()
{
    for (int precision : {1, 2}) {
        identityTests(precision);
        initializationTests(precision);
        reciprocalWitness(precision);
    }
    return reportCheckResults();
}
