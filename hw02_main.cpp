
#include "bisection.h"
#include "newton_raphson.h"
#include "functions01.h"
#include <stdio.h>
#include <math.h>

const double PI = 3.141592653589793238462;

int hw02_main() {
    double solution;
    double lowerBound = 0;
    double upperBound = PI / 2;
    int max_iter = 1e4;
    double TOL = 0.5e-8;
    // exercise 2.1
    if (false) {
        solution = bisect(exercise_2_1, lowerBound, upperBound, max_iter, TOL);
    }
    // exercise 2.3
    if (false) {
        solution = bisect_r(exercise_2_1, lowerBound, upperBound, 1, max_iter, TOL);
    }
    // exercise 2.5
    // legendre polynomial, n = 8: 0.2734375 - 9.84375x^2.0 + 54.140625x^4.0 - 93.84375x^6.0 + 50.2734375x^8.0
    if (false) {
        lowerBound = 0;
        upperBound = 0.2;
        solution = bisect(l8, lowerBound, upperBound, max_iter, TOL);
        double guess = 0.2;
        solution = newtonsolve(l8, l8deriv, guess, max_iter, TOL); // solution: 0.1834346424796
    }

    // exercise 2.13 plot a function of energy
    {
        double a = 0.3; // nm (nanometer)
        double m = 1; // m_e (electron mass unit)
        double v_0 = 10; // eV (electron volt unit)
        double params[3] = {0.3, 1, 10};
        double start = 0;
        double stop = 10;
        double step = 0.1;
        if (false) {
            for (double energy = start; energy < stop; energy += step) {
                solution = bc_well_even(energy, params);
                printf("%3.2f\r\n", solution);
            }
            for (double energy = start; energy < stop; energy += step) {
                solution = bc_well_odd(energy, params);
                printf("%3.2f\r\n", energy);
            }
        }
        if (false) {
            for (double energy = start; energy < stop; energy += step) {
                solution = bc_well_odd(energy, params);
                printf("%3.2f\r\n", solution);
            }
        }
        if (true) {
            // exercise 2.14
            // find the lowest even and odd energies, using the graphs as a guide
            double energy_even = bisect_params(bc_well_even, params, lowerBound, upperBound, max_iter, TOL);
            lowerBound = 2;
            upperBound = 3;
            double energy_odd = bisect_params(bc_well_odd, params, lowerBound, upperBound, max_iter, TOL);
            double alpha = get_alpha(m, energy_even);
            double beta = get_beta(m, v_0, energy_even);
            // see what the alpha and beta are
            printf("a = %3.2f\r\n", alpha);
            printf("b = %3.2f\r\n", beta);
            // III: C * exp(beta*x) | I: B * cos(alpha*x) | II: C * exp(-beta*x))
            // let's set C = 1 for simplicity, and then renormalize based on that
            // integral = sq anti III(infinity to -a) + sq anti I (-a to a) + sq anti II (a to infinity))
            double C = 1;
            double A = 0;
            // first use if C = 1, by matching the function values at x = -a
            // find B: exp(beta*-a))= B cos(alpha*-a) implies B = exp(beta*-a)/cos(alpha*-a)
            double B = exp(beta*-a) / cos(alpha*-a);
            printf("If C = 1, then B = %3.2f\r\n", B);
            // now we can sum that integral
            double integral = wavefunction_02sqint(a, C, beta);
            printf("If C = 1 and B = %3.2f, then the integral = %3.2f\r\n", B, integral);
            // now normalize the integral to one, i.e. divide the coefficients by the value of the integral
            C /= integral;
            B /= integral;
            // check that the integral now equals one
            integral = wavefunction_03sqint(-a, C, beta) + wavefunction_01sqint(a, A, B, alpha) - wavefunction_01sqint(-a, A, B, alpha) + wavefunction_02sqint(a, C, beta);
            printf("integral = %3.2f\r\n", integral);
            // see what the coefficients B and C are
            printf("C = %3.2f\r\n", C);
            printf("B = %3.2f\r\n", B);
            // plot wavefunction
            double output;
            if (true) {
                for (double i = -a; i <= a; i = i + 0.1) {
                    output = wavefunction_01(i, 0, 1, alpha);
                    printf("%3.6f\r\n", output);
                }
            }
        }

    }
    return 0;
}
