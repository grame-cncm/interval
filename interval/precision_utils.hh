#pragma once

#include <cmath>
#include "check.hh"

/** 
 * truncates x at the precision induced by lsb
 * */
double truncate(double x, int lsb)
{
    double epsilon = pow(2, lsb);
    double res = epsilon*(double)floor(x/epsilon);
    return res;
}

/**
 * @brief compute the precision needed in the output of a function 
 * 
 * @param f the function to analyse
 * @param x the point at which the tightest precision is needed
 * @param epsilon the signed gap between the two consecutive numbers at which to compute the precision
*/
int exactPrecisionUnary(ufun f, double x, double epsilon)
{
    int res = floor((double)log2(std::abs(f(x + epsilon) - f(x))));
    return res;
}