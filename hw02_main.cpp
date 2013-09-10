
#include "bisection.h"
#include "newton_raphson.h"
#include "functions01.h"
#include <stdio.h>
#include <math.h>

const double PI = 3.141592653589793238462;
const bool DEBUG = false;
void exercise_2_15(double, double);
double smallest_width(double);

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
    if (true) {
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
        double params[3] = {a, 1, 10};
        double start = 0;
        double stop = 10;
        double step = 0.1;
        if (false) {
            for (double energy = start; energy < stop; energy += step) {
                solution = bc_well_even(energy, params);
                printf("%3.2f\r\n", solution);
            }
        }
        if (false) {
            for (double energy = start; energy < stop; energy += step) {
                solution = bc_well_odd(energy, params);
                printf("%3.2f\r\n", solution);
            }
        }
        if (false) {
            // exercise 2.14
            // find the lowest even and odd energies, using the graphs as a guide
            double energy_even = bisect_params(bc_well_even, params, lowerBound, upperBound, max_iter, TOL);
            lowerBound = 2;
            upperBound = 3;
            double energy_odd = bisect_params(bc_well_odd, params, lowerBound, upperBound, max_iter, TOL);
            double alpha = get_alpha(m, energy_even);
            double beta = get_beta(m, v_0, energy_even);
            // see what the alpha and beta are
            printf("\r\nEven Solution\r\n");
            printf("a = %3.2f\r\n", alpha);
            printf("b = %3.2f\r\n", beta);
            // III: C * exp(beta*x) | I: B * cos(alpha*x) | II: C * exp(-beta*x))
            // let's set C = 1 for simplicity, and then renormalize based on that
            // integral = sq anti III(infinity to -a) + sq anti I (-a to a) + sq anti II (a to infinity))
            double C = 1;
            double A = 0;
            // first use if C = 1, by matching the function values at x = -a
            // find B: exp(beta*-a))= B cos(alpha*-a) implies B = exp(beta*-a)/cos(alpha*-a)
            double B = C * exp(beta*-a) / cos(alpha*-a);
            printf("If C = 1, then B = %3.8f\r\n", B);
            // now we can sum that integral
            double integral = wavefunction_03sqint(-a, C, beta) - wavefunction_02sqint(a, C, beta) + wavefunction_01sqint(a, A, B, alpha) - wavefunction_01sqint(-a, A, B, alpha);
            printf("If C = 1 and B = %3.8f, then the integral = %3.8f\r\n", B, integral);
            // now normalize the integral to one
            double params2[4] = {a, 1, 10, energy_even};
            C = bisect_params(get_integral_even, params2, 1, 1000, max_iter, TOL);
            B = C * exp(beta*-a) / cos(alpha*-a);
            // check that the integral now equals one
            integral = wavefunction_03sqint(-a, C, beta) + wavefunction_01sqint(a, A, B, alpha) - wavefunction_01sqint(-a, A, B, alpha) + wavefunction_02sqint(a, C, beta);
            printf("integral = %3.2f\r\n", integral);
            // see what the coefficients B and C are
            printf("C = %3.2f\r\n", C);
            printf("B = %3.2f\r\n", B);
            // plot wavefunction
            double output;
            if (false) {
                for (double i = -1; i < -a + 0.01; i = i + 0.01) {
                    output = wavefunction_03(i, C, beta);
                    printf("%3.6f\r\n", output);
                }

                printf("BREAK\r\n");
                for (double i = -a; i < a + 0.01; i = i + 0.01) {
                    output = wavefunction_01(i, 0, B, alpha);
                    printf("%3.6f\r\n", output);
                }
                printf("BREAK\r\n");
                for (double i = a; i < 1 + 0.01; i = i + 0.01) {
                    output = wavefunction_02(i, C, beta);
                    printf("%3.6f\r\n", output);
                }
            }

            // now for the odd solution
            printf("\r\nEven Solution\r\n");
            alpha = get_alpha(m, energy_odd);
            beta = get_beta(m, v_0, energy_odd);
            // see what the alpha and beta are
            printf("a = %3.2f\r\n", alpha);
            printf("b = %3.2f\r\n", beta);
            // III: C * exp(beta*x) | I: B * cos(alpha*x) | II: C * exp(-beta*x))
            // let's set C = 1 for simplicity, and then renormalize based on that
            // integral = sq anti III(infinity to -a) + sq anti I (-a to a) + sq anti II (a to infinity))
            C = 1;
            B = 0;
            // first use if C = 1, by matching the function values at x = -a
            // find B: exp(beta*-a))= B cos(alpha*-a) implies B = exp(beta*-a)/cos(alpha*-a)
            A = C * exp(beta*-a) / cos(alpha*-a);
            double D = -C;
            printf("If C = 1, then A = %3.8f\r\n", A);
            // now we can sum that integral
            integral = wavefunction_03sqint(-a, C, beta) - wavefunction_02sqint(a, D, beta) + wavefunction_01sqint(a, A, B, alpha) - wavefunction_01sqint(-a, A, B, alpha);
            printf("If C = 1 and A = %3.8f, then the integral = %3.8f\r\n", A, integral);
            // now normalize the integral to one
            double params3[4] = {a, 1, 10, energy_odd};
            C = bisect_params(get_integral_odd, params3, 1, 1000, max_iter, TOL);
            A = C * exp(beta*-a) / sin(alpha*-a);
            // check that the integral now equals one
            integral = wavefunction_03sqint(-a, C, beta) + wavefunction_01sqint(a, A, B, alpha) - wavefunction_01sqint(-a, A, B, alpha) + wavefunction_02sqint(a, D, beta);
            printf("integral = %3.2f\r\n", integral);
            // see what the coefficients B and C are
            printf("C = %3.2f\r\n", C);
            printf("A = %3.2f\r\n", A);
            // plot wavefunction
            output;
            if (false) {
                for (double i = -1; i < -a + 0.01; i = i + 0.01) {
                    output = wavefunction_03(i, C, beta);
                    printf("%3.6f\r\n", output);
                }

                printf("BREAK\r\n");
                for (double i = -a; i < a + 0.01; i = i + 0.01) {
                    output = wavefunction_01(i, A, 0, alpha);
                    printf("%3.6f\r\n", output);
                }
                printf("BREAK\r\n");
                for (double i = a; i < 1 + 0.01; i = i + 0.01) {
                    output = wavefunction_02(i, -C, beta);
                    printf("%3.2f\r\n", i);
                    printf("%3.6f\r\n", output);
                }
            }

        }

    }
    // test dependence on width
    if (false) {
        for (double i = 0.1; i < 0.15; i += 0.01) {
            exercise_2_15(i, 10);
        }
    }
    // test dependence on potential
      if (true) {
        for (double i = 1; i < 2; i += 0.5) {
            exercise_2_15(0.3, i);
        }
    }
    //bisect(smallest_width, 0.9, 0.11, max_iter, TOL);
    //exercise_2_15(0.12,10);
    
    return 0;
}

double smallest_width(double half_width) {
    double solution;
    double integral_start = -10;
    double integral_end = 100;
    int max_iter = 1e4;
    double TOL = 0.5e-8;
    double a = half_width; // nm (nanometer)
    double m = 1; // m_e (electron mass unit)
    double v_0 = 10; // eV (electron volt unit)
    double params[4];
    params[0] = a;
    params[1] = m;
    params[2] = v_0;
    double start = 0;
    double stop = 10;
    double step = 0.1;
    double lowerBound = start;
    double upperBound = stop;

    printf("\r\n using a = %3.2f and v_0 = %3.2f\r\n", a, v_0);

    // find the lowest even and odd energies, using the graphs as a guide
    double signvalue = 1;
    if (true) {
        //printf("\r\nodd spectral\r\n");
        for (double energy = start; energy < stop; energy += step) {
            solution = bc_well_odd(energy, params);
            signvalue *= solution;
            //printf("%3.2f\r\n", solution);
            // check at which energy the spectral function changes sign
            if (signvalue >= 0) {
                signvalue = 1;
            } else {
                lowerBound = energy - step;
                upperBound = energy;
                break;
            }
        }
    }
    double energy_odd = bisect_params(bc_well_odd, params, lowerBound, upperBound, max_iter, TOL);
    return energy_odd - 10;
}

void exercise_2_15(double half_width, double potential) {
    double solution;
    double integral_start = -10;
    double integral_end = 100;
    int max_iter = 1e4;
    double TOL = 0.5e-8;
    double a = half_width; // nm (nanometer)
    double m = 1; // m_e (electron mass unit)
    double v_0 = potential; // eV (electron volt unit)
    double params[4];
    params[0] = a;
    params[1] = m;
    params[2] = v_0;
    double start = 0;
    double stop = 10;
    double step = 0.1;
    double lowerBound = start;
    double upperBound = stop;

    //printf("\r\n using a = %3.2f and v_0 = %3.2f\r\n", a, v_0);

    // find the lowest even and odd energies, using the graphs as a guide
    double signvalue = 1;
    if (true) {
        //printf("\r\neven spectral\r\n");
        for (double energy = start; energy < stop; energy += step) {
            solution = bc_well_even(energy, params);
            signvalue *= solution;
            //printf("%3.2f\r\n", solution);
            // check at which energy the spectral function changes sign
            if (signvalue >= 0) {
                signvalue = 1;
            } else if (signvalue < 0){
                lowerBound = energy - step;
                upperBound = energy;
                break;
            }
        }
    }
    double energy_even = bisect_params(bc_well_even, params, lowerBound, upperBound, max_iter, TOL);
    signvalue = 1;
    if (true) {
        //printf("\r\nodd spectral\r\n");
        for (double energy = start; energy < stop; energy += step) {
            solution = bc_well_odd(energy, params);
            signvalue *= solution;
            //printf("%3.2f\r\n", solution);
            // check at which energy the spectral function changes sign
            if (signvalue >= 0) {
                signvalue = 1;
            } else if (signvalue < 0){
                lowerBound = energy - step;
                upperBound = energy;
                break;
            }
        }
    }
    lowerBound = 0.5; upperBound = 1.01;
    double energy_odd = bisect_params(bc_well_odd, params, lowerBound, upperBound, max_iter, TOL);
    printf("%3.6f\r\n%3.6f\r\n", energy_even, energy_odd);
    // exercise 2.14
    if (false) {
        double alpha = get_alpha(m, energy_even);
        double beta = get_beta(m, v_0, energy_even);
        // see what the alpha and beta are
        printf("\r\nEven Solution\r\n");
        printf("a = %3.2f\r\n", alpha);
        printf("b = %3.2f\r\n", beta);
        // III: C * exp(beta*x) | I: B * cos(alpha*x) | II: C * exp(-beta*x))
        // let's set C = 1 for simplicity, and then renormalize based on that
        // integral = sq anti III(infinity to -a) + sq anti I (-a to a) + sq anti II (a to infinity))
        double C = 1;
        double A = 0;
        // first use if C = 1, by matching the function values at x = -a
        // find B: exp(beta*-a))= B cos(alpha*-a) implies B = exp(beta*-a)/cos(alpha*-a)
        double B = C * exp(beta*-a) / cos(alpha*-a);
        printf("If C = 1, then B = %3.8f\r\n", B);
        // now we can sum that integral
        double integral = wavefunction_03sqint(-a, C, beta) - wavefunction_02sqint(a, C, beta) + wavefunction_01sqint(a, A, B, alpha) - wavefunction_01sqint(-a, A, B, alpha);
        printf("If C = 1 and B = %3.8f, then the integral = %3.8f\r\n", B, integral);
        // now normalize the integral to one
        double params2[4] = {a, 1, 10, energy_even};
        C = bisect_params(get_integral_even, params2, 1, 1000, max_iter, TOL);
        B = C * exp(beta*-a) / cos(alpha*-a);
        // check that the integral now equals one
        integral = wavefunction_03sqint(-a, C, beta) + wavefunction_01sqint(a, A, B, alpha) - wavefunction_01sqint(-a, A, B, alpha) + wavefunction_02sqint(a, C, beta);
        printf("integral = %3.2f\r\n", integral);
        // see what the coefficients B and C are
        printf("C = %3.2f\r\n", C);
        printf("B = %3.2f\r\n", B);
        // plot wavefunction
        double output;
        if (false) {
            printf("\r\nEven wavefunction\r\n");
            for (double i = -1; i < -a; i = i + 0.01) {
                output = wavefunction_03(i, C, beta);
                printf("%3.6f\r\n", output);
            }
            for (double i = -a; i < a; i = i + 0.01) {
                output = wavefunction_01(i, 0, B, alpha);
                printf("%3.6f\r\n", output);
            }
            for (double i = a; i < 1 + 0.01; i = i + 0.01) {
                output = wavefunction_02(i, C, beta);
                printf("%3.6f\r\n", output);
            }
        }

        // now for the odd solution
        printf("\r\nOdd Solution\r\n");
        alpha = get_alpha(m, energy_odd);
        beta = get_beta(m, v_0, energy_odd);
        params[4] = energy_odd;
        // see what the alpha and beta are
        printf("a = %3.2f\r\n", alpha);
        printf("b = %3.2f\r\n", beta);
        // III: C * exp(beta*x) | I: B * cos(alpha*x) | II: C * exp(-beta*x))
        // let's set C = 1 for simplicity, and then renormalize based on that
        // integral = sq anti III(infinity to -a) + sq anti I (-a to a) + sq anti II (a to infinity))
        C = 1;
        B = 0;
        // first use if C = 1, by matching the function values at x = -a
        // find B: exp(beta*-a))= B cos(alpha*-a) implies B = exp(beta*-a)/cos(alpha*-a)
        A = C * exp(beta*-a) / cos(alpha*-a);
        double D = -C;
        printf("If C = 1, then A = %3.8f\r\n", A);
        // now we can sum that integral
        integral = wavefunction_03sqint(-a, C, beta) - wavefunction_02sqint(a, D, beta) + wavefunction_01sqint(a, A, B, alpha) - wavefunction_01sqint(-a, A, B, alpha);
        printf("If C = 1 and A = %3.8f, then the integral = %3.8f\r\n", A, integral);
        // now normalize the integral to one
        double params3[4] = {a, 1, 10, energy_odd};
        C = bisect_params(get_integral_odd, params3, 1, 1000, max_iter, TOL);
        A = C * exp(beta*-a) / sin(alpha*-a);
        // check that the integral now equals one
        integral = wavefunction_03sqint(-a, C, beta) + wavefunction_01sqint(a, A, B, alpha) - wavefunction_01sqint(-a, A, B, alpha) + wavefunction_02sqint(a, D, beta);
        printf("integral = %3.2f\r\n", integral);
        // see what the coefficients B and C are
        printf("C = %3.2f\r\n", C);
        printf("A = %3.2f\r\n", A);
        // plot wavefunction
        if (false) {
            printf("\r\nOdd wavefunction\r\n");
            for (double i = -1; i < -a; i = i + 0.01) {
                output = wavefunction_03(i, C, beta);
                printf("%3.6f\r\n", output);
            }
            for (double i = -a; i < a; i = i + 0.01) {
                output = wavefunction_01(i, A, 0, alpha);
                printf("%3.6f\r\n", output);
            }
            for (double i = a; i < 1 + 0.01; i = i + 0.01) {
                output = wavefunction_02(i, -C, beta);
                printf("%3.6f\r\n", output);
            }
        }
    }
}