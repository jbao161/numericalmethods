
#include "bisection.h"
#include "newton_raphson.h"
#include "functions01.h"
#include <stdio.h>
#include <math.h>

const double PI  =3.141592653589793238462;
int hw02_main()
{
    double solution;
    double lowerBound = 0;
    double upperBound = PI / 2;
    int max_iter = 1e4;
    double TOL = 0.5e-8;
    // exercise 2.1
    if (false)
    {
        solution = bisect(exercise_2_1, lowerBound, upperBound, max_iter, TOL);
    }
    // exercise 2.3
    if (false)
    {
        solution = bisect_r(exercise_2_1, lowerBound, upperBound, 1, max_iter, TOL);
    }
    // exercise 2.5
    // legendre polynomial, n = 8: 0.2734375 - 9.84375x^2.0 + 54.140625x^4.0 - 93.84375x^6.0 + 50.2734375x^8.0
    if (false)
    {
        lowerBound = 0;
        upperBound = 0.2;
        solution = bisect(l8, lowerBound, upperBound, max_iter, TOL);
        double guess = 0.2;
        solution = newtonsolve(l8, l8deriv, guess, max_iter, TOL); // solution: 0.1834346424796
    }

    // exercise 2.13 plot a function of energy
    if (true)
    {
        double a = 0.3; // nm (nanometer)
        double m = 1; // m_e (electron mass unit)
        double v_0 = 10; // eV (electron volt unit)
        double params[3] = {0.3, 1, 10};
        double start = 0;
        double stop = 10;
        double step = 0.5;
        for (double energy = start; energy < stop; energy += step)
        {
            solution = bc_well_even(energy, params);
            printf("%3.2f\r\n", solution);
        }
        for (double energy = start; energy < stop; energy += step)
        {
            solution = bc_well_odd(energy, params);
            printf("%3.2f\r\n", solution);
        }


        // exercise 2.14
        double energy_even = bisect_params(bc_well_even, params, lowerBound, upperBound, max_iter, TOL);
        lowerBound = 1;
        upperBound = 3;
        double energy_odd = bisect_params(bc_well_odd, params, lowerBound, upperBound, max_iter, TOL);
        double alpha1 = alpha(m, energy_even);
        double beta1 = beta(m, energy_even, v_0);
        double output;
        for (double i = 0.0; i < 10; i = i + 0.5)
        {
            output = wavefunction_01(i, 0, 1, alpha1, beta1);
            printf("%3.6f\r\n", output);
        }
    }
    return 0;
}
