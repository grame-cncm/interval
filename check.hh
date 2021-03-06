#pragma once

#include <string>

#include "interval.hh"
#include "interval_algebra.hh"

void check(const std::string& testname, const interval& exp);

void check(const std::string& testname, const interval& exp, const interval& res);

void check(const std::string& testname, bool exp, bool res);

using ufun = double (*)(double);
using bfun = double (*)(double, double);

using umth = interval (interval_algebra::*)(const interval& x) const;
using bmth = interval (interval_algebra::*)(const interval& x, const interval& y) const;

interval testfun(int N, bfun f, const interval& x, const interval& y);
void     analyzemod(interval x, interval y);

void analyzeUnaryFunction(int E, int M, const char* title, const interval& D, ufun f);
void analyzeUnaryMethod(int E, int M, const char* title, const interval& D, ufun f, umth m);

void analyzeBinaryMethod(int E, int M, const char* title, const interval& C, const interval& D, bfun f, bmth m);
