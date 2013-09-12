#include "newton_raphson.h"
#include "functions02.h"
#include <stdio.h>

const int max_iter = 1e4;
const double TOL = 5e-8;

void test01();
void test02(double*);

void hw03_main() {
    if (false) {
        test01();
    }
    if (true) {
        double in_params[3][2] = {
            {1, -27},
            {20, 0},
            {15, -13}
        };
        double out_params[3];
        convert_inputs(in_params, out_params);
        test02(out_params);
    }
    if (true) {

    }
}

// checks the hybrid newton method

void test01() {
    {
        double solution;
        double lowerBound = 0;
        double upperBound = 3;
        // hybrid method takes it more quickly to the solution
        solution = newtonhybrid(poly2test, poly2test_d, lowerBound, upperBound, max_iter, TOL);
        solution = newtonsolve(poly2test, poly2test_d, (lowerBound + upperBound) / 2, max_iter, TOL);
        // hybrid method brackets the intended solution, newton finds a solution far away
        solution = newtonsolve(poly2test, poly2test_d, lowerBound, max_iter, TOL);
    }
}

// finds the energy eigenvalues for ebb1, the half-infinite well

void test02(double params[3]) {
    double step_percent = 0.01;
    double step = params[2] * step_percent;
    for (double energy = 0; energy < params[2]+2; energy += step) {
        printf("energy: %3.6f ; function: %3.6f\r\n", energy, energy_function(energy, params));
    }
    newtonhybridp(energy_function, energy_function_d, params, 0.000000001, 0.8, max_iter, TOL);
}
