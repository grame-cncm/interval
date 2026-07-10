/* Copyright 2023 Yann ORLAREY
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>

#include "interval/check.hh"
#include "interval/interval_algebra.hh"
#include "interval/interval_def.hh"

using namespace itv;

int main()
{
    // test interval representation
    // check("interval()", interval());
    check("interval(0,100,-24)", interval(100.0, 0.0));
    check("interval(0,0,0)", interval(0, 0));
    check("interval(-10,0,-24)", interval(0, -10));
    check("interval(-10,10,-5)", interval(-10, 10, -5));

    // test union intersection

    check("test2", intersection(interval(0, 100), interval(-10, 0)), interval(0));
    check("test3", intersection(interval(10, 100), interval(-10, 0)), interval());
    check("test4", reunion(interval(0, 100), interval(-100, 50)), interval(-100, 100));
    check("test5", reunion(interval(0, 100), interval(10, 500)), interval(0, 500));

    // test predicates
    interval a(1, 100), b(10, 20), c(-10, 0), n;

    std::cout << a << " == " << b << " = " << (a == b) << std::endl;
    std::cout << a << " <= " << b << " = " << (a <= b) << std::endl;
    std::cout << a << " < " << b << " = " << (a < b) << std::endl;

    std::cout << std::endl;

    std::cout << a << " != " << b << " = " << (a != b) << std::endl;
    std::cout << a << " > " << b << " = " << (a > b) << std::endl;
    std::cout << a << " >= " << b << " = " << (a >= b) << std::endl;

    std::cout << std::endl;

    std::cout << a << " == " << a << " = " << (a == a) << std::endl;
    std::cout << a << " <= " << a << " = " << (a <= a) << std::endl;
    std::cout << a << " < " << a << " = " << (a < a) << std::endl;

    std::cout << std::endl;

    std::cout << a << " != " << a << " = " << (a != a) << std::endl;
    std::cout << a << " > " << a << " = " << (a > a) << std::endl;
    std::cout << a << " >= " << a << " = " << (a >= a) << std::endl;

    std::cout << std::endl;

    std::cout << a << " != " << n << " = " << (a != n) << std::endl;
    std::cout << a << " > " << n << " = " << (a > n) << std::endl;
    std::cout << a << " >= " << n << " = " << (a >= n) << std::endl;

    std::cout << std::endl;

    std::cout << "must be true " << (reunion(a, n) == a) << std::endl;
    std::cout << "must be true " << (intersection(a, n) == n) << std::endl;

    std::cout << std::endl;

    // interval v = singleton(2);
    std::cout << singleton(2.5) << std::endl;
    std::cout << singleton(2.25) << std::endl;
    std::cout << singleton(2.1) << std::endl;
    std::cout << singleton(24) << std::endl;
    std::cout << singleton(44) << std::endl;
    std::cout << singleton(255) << std::endl;

    interval_algebra A;
    // A.testAll();
    /*A.testExp();
    A.testLog();*/
    /*A.testAcos();
    A.testAcosh();
    A.testAsin();*/
    // A.testAsinh();
    // A.testAtan();
    // A.testAtanh();
    // A.testCosh();
    // A.testLog10();
    // A.testLog();
    // A.testSin();
    // A.testCos();
    // A.testSin();
    // A.testSinh();
    // A.testSqrt();
    // A.testTan();
    // A.testTanh();
    A.testPow();
    // A.testAdd();
    // A.testMul();
    // A.testSub();
    // A.testInv();
    // A.testDiv();
    // A.testAtan2();

    // A.testAnd();
    // A.testNot();
    // A.testOr();
    // A.testXor();

    // A.testMod();

    /* A.testEq();
    A.testGe();
    A.testGt();
    A.testLe();
    A.testLt();*/

    // A.testAdd();
    // A.testMul();

    {
        double u = 0.0;
        double v = nextafter(u, -HUGE_VAL);
        double w = nextafter(v, -HUGE_VAL);

        if ((u != v) && (v != w) && (u > v) && (v > w)) {
            std::cout << "Order OK" << std::endl;
        } else {
            std::cout << "Order NOT OK" << std::endl;
        }
    }

    analyzemod(interval(9), interval(10));
    analyzemod(interval(8, 9), interval(10));
    analyzemod(interval(8, 9), interval(1, 10));
    analyzemod(interval(8, 11), interval(10));
    analyzemod(interval(9), interval(9, 10));
    analyzemod(interval(-9), interval(9, 10));
    analyzemod(interval(0, 9), interval(9, 10));
    analyzemod(interval(-9, 0), interval(9, 10));
    analyzemod(interval(0, 100), interval(1));
    analyzemod(interval(-100, 100), interval(1));

    analyzeUnaryFunction(10, 1000, "rint", interval(-10000, 10000), rint);
    analyzeUnaryFunction(10, 1000, "floor", interval(-10000, 10000), floor);
    analyzeUnaryFunction(10, 1000, "ceil", interval(-10000, 10000), ceil);

    /* interval X = interval(-10, 10, -10);
    propagateBackwardsUnaryMethod("exp", &interval_algebra::Exp, X, -24);
    std::cout << std::endl;
    analyzeUnaryMethod(10, 40000, "exp", X, exp, &interval_algebra::Exp);*/

    /* interval X = interval(-10, 10, -24);
    std::vector<const char*> titles{ "sin", "exp", "cos"};
    std::vector<umth> mps{&interval_algebra::Sin, &interval_algebra::Exp, &interval_algebra::Cos};
    std::cout << std::endl << "sin(exp(cos("<< X << "))) = " << A.Sin(A.Exp(A.Cos(X))) << std::endl << std::endl;
    propagateBackwardsComposition(titles, mps, X, -24);
    std::cout << std::endl << "sin(exp(cos("<< X << "))) = " << A.Sin(A.Exp(A.Cos(X))) << std::endl;
    std::cout << "----------------" << std::endl << std::endl;

    X = interval(1, 10, -24);
    std::vector<const char*> titles2{"log", "exp"};
    std::vector<umth> mps2{ &interval_algebra::Log, &interval_algebra::Exp};
    std::cout << std::endl << "log(exp("<< X << ")) = " << A.Log(A.Exp(X)) << std::endl << std::endl;
    propagateBackwardsComposition(titles2, mps2, X, -24);
    std::cout << std::endl << "log(exp("<< X << ")) = " << A.Log(A.Exp(X)) << std::endl;
    std::cout << "----------------" << std::endl << std::endl;

    X = interval(1, 10, -24);
    std::vector<const char*> titles3{"exp", "log"};
    std::vector<umth> mps3{ &interval_algebra::Exp, &interval_algebra::Log};
    std::cout << std::endl << "exp(log("<< X << ")) = " << A.Exp(A.Log(X)) << std::endl << std::endl;
    propagateBackwardsComposition(titles3, mps3, X, -24);
    std::cout << std::endl << "exp(log("<< X << ")) = " << A.Exp(A.Log(X)) << std::endl;
    std::cout << "----------------" << std::endl << std::endl; */

    /* interval X = interval(1, 2, -24);
    interval Y = interval(1, 10, -24);
    std::cout << "pow(" << X << ", " << Y << ") = " << A.Pow(X, Y) << std::endl << std::endl;
    propagateBackwardsBinaryMethod("pow", &interval_algebra::Pow, X, Y, -24);
    std::cout << std::endl;
    std::cout << "pow(" << X << ", " << Y << ") = " << A.Pow(X, Y) << std::endl;
    interval Z = A.Pow(X, Y);
    std::cout << "exp(" << Z << ") = " <<  A.Exp(Z) << std::endl;
    propagateBackwardsUnaryMethod("exp", &interval_algebra::Exp, Z, -24);
    std::cout << "exp(" << Z << ") = " <<  A.Exp(Z) << std::endl;
    std::cout << "----------------" << std::endl << std::endl;*/

    /* interval X = interval(0, 2, -6);
    std::cout << "sin(" << X << ") = "  << A.Sin(X) << std::endl << std::endl;
    X = itv::interval(X.lo(), X.hi(), 2*X.lsb());
    std::cout << "sin(" << X << ") = "  << A.Sin(X) << std::endl << std::endl;*/
    /* propagateBackwardsUnaryMethod("sin", &interval_algebra::Sin, X, -24);
    std::cout << std::endl;
    std::cout << "sin(" << X << ") = "  << A.Sin(X) << std::endl << std::endl;*/

    /* interval X = interval(1, 2, -24);
    interval Y = interval(1, 10, -24);
    std::cout << "rsh(" << X << ", " << Y << ") = " << A.Rsh(X, Y) << std::endl << std::endl;
    propagateBackwardsBinaryMethod("rsh", &interval_algebra::Rsh, X, Y, -20);
    std::cout << std::endl;
    std::cout << "rsh(" << X << ", " << Y << ") = " << A.Rsh(X, Y) << std::endl; */

    // Addition test (wrapping)
    /* interval X = interval(1, HUGE_VAL, 0);
    interval Y = interval(1, HUGE_VAL, 0);
    std::cout << X << " / " << Y << " = " << A.Div(X, Y) << std::endl;*/

    /* interval X = interval(-HUGE_VAL, HUGE_VAL, -780);
    interval Y = interval(-HUGE_VAL, HUGE_VAL, -804);
    std::cout << X << " / " << Y << " = " << A.Div(X, Y) << std::endl;*/

    /* interval X = interval(-HUGE_VAL, HUGE_VAL, -23);
    std::cout << X << std::endl;
    interval Y = A.Inv(X);
    std::cout << Y << std::endl;
    std::cout << A.Inv(Y) << std::endl;
    X = interval(-HUGE_VAL, HUGE_VAL, -20);
    std::cout << X << std::endl;
    std::cout << A.Inv(X) << std::endl;*/
    // propagateBackwardsUnaryMethod("inv", &interval_algebra::Inv, X, -24);
    // std::cout << X << std::endl;

    /* interval X = interval(10, 10, -29);
    interval Y = interval(-HUGE_VAL, HUGE_VAL, -61);
    std::cout << "Pow(" << X << ", " << Y << ") = " << A.Pow(X, Y) << std::endl;*/

    /* interval X = interval(0.5,96000,-65);
    interval Y = interval(23.4489,23.4489,-958);
    std::cout << X << "/" << Y << " = " << A.Div(X, Y) << std::endl;*/
}