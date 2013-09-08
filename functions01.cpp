#include "functions01.h"
#include <math.h>

const double h_sqrd = 0.076199682; // me ev nm^2

// exercise 2.1: cos(x) - x

double exercise_2_1(double x)
{
    return cos(x) - x;
}
double l8(double x)
{
    return 0.2734375 - 9.84375*pow(x,2.0) + 54.140625*pow(x,4.0) - 93.84375*pow(x,6.0) + 50.2734375*pow(x,8.0);
}
double l8deriv(double x)
{
    return - 9.84375*2*x + 54.140625*4*pow(x,3.0) - 93.84375*6*pow(x,5.0) + 50.2734375*8*pow(x,7.0);
}

double alpha(double mass, double Energy)
{
    return sqrt(2 * mass * Energy / h_sqrd);
}

double beta(double mass, double potential, double Energy)
{
    return sqrt(2 * mass * (potential - Energy) / h_sqrd);
}

// f(E) = beta * cos(alpha * a) - alpha * sin(alpha * a)

double bc_well_even(double E, double *params)
{
    double a = params[0];
    double m = params[1];
    double v_0 = params[2];
    double beta = sqrt(2 * m * (v_0 - E) / h_sqrd);
    double alpha = sqrt(2 * m * E / h_sqrd);
    return beta * cos(alpha * a) - alpha * sin(alpha * a);
}

double bc_well_odd(double E, double *params)
{
    double a = params[0];
    double m = params[1];
    double v_0 = params[2];
    double beta = sqrt(2 * m * (v_0 - E) / h_sqrd);
    double alpha = sqrt(2 * m * E / h_sqrd);
    return alpha * cos(alpha * a) + beta * sin(alpha * a);
}

// solution inside the box

double wavefunction_01(double x, double A, double B, double alpha, double beta)
{
    return A * sin(alpha * x) + B * cos(alpha * x);
}
// solution outside the box, to the right

double wavefunction_02(double x, double D, double beta)
{
    return D * exp(-beta * x);
}

// solution outside the box, to the left

double wavefunction_03(double x, double D, double beta)
{
    return D * exp(beta * x);
}
