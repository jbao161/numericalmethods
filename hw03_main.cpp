#include "newton_raphson.h"
#include "functions02.h"
#include <stdio.h>
#include <vector>
#include <utility>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "bisection.h"

using namespace std;
const int max_iter = 1e4;
const double TOL = 5e-16;
const bool VERBOSE = false;
const double hbar_sqrd = 0.076199682; // me ev nm^2
void test01();
int test02(double*, vector<double>&);
void test03(double*, int);

double get_vvp(double potential, double *params, int num_energies) {
    // for a given L, start plotting the energy function until the specified number of energy roots are found
    params[2] = potential;
    double step_percent = 0.01;
    double step = potential * step_percent;
    int count = 0;
    double lowerBound = 0;
    int num_of_brackets = 0;
    double bracket_first, bracket_second;
    double root;
    double f_E;
    double f_E_prev = energy_function(0, params);
    double highest_energy;
    for (double energy = step; energy < potential; energy += step) {
        f_E = energy_function(energy, params);
        if (signbit(f_E) != signbit(f_E_prev)) {
            count++;
            if (count == num_energies) {
                lowerBound = energy - step;
                bracket_first = lowerBound;
                bracket_second = energy;
                num_of_brackets = num_of_brackets + 1;
                highest_energy = newtonhybridp(energy_function, energy_function_d, params, bracket_first, bracket_second, max_iter, TOL, VERBOSE);
            }
            // printf("num brackets: %d\r\n", num_of_brackets);
        }
        f_E_prev = f_E;
        if (false) {
            printf("energy: %3.8f ; function: %3.8f\r\n", energy, f_E);
        }
    }
    // printf("length: %3.8f\r\n", L);
    // printf("# roots: %d\r\n", count);
    double distance = -1;
    if (count == num_energies) {
        distance = (potential - highest_energy);
    }
    printf("%3.8f,%d,%3.8f\r\n", potential, count, distance);
    if (count > num_energies) {
        return -1;
    }
    if (count < num_energies) {
        return -2;
    }
    return (potential - highest_energy);
    // bisection search with each increment step that changes sign
}

void hw03_main() {
    if (false) {
        test01();
    }

    double in_params[3][2] = {
        {1, -27},
        {20, 0},
        {15, -13}
    };
    double out_params[6];
    convert_inputs(in_params, out_params);
    vector<double> energies;
    int num_eigenenergies = test02(out_params, energies); // put the eigenvalues in outparams
    double A; // wavefunction_01 coefficient
    double D; // wavefunction_02 coefficient
    double L; // length of the well
    // double v0 = out_params[2]; // potential of well
    double step = 0.01;
    for (int i = 0; i < num_eigenenergies; i++) {
        out_params[3] = energies.at(i);

        double m = out_params[0];
        double L = out_params[1];
        double v0 = out_params[2];
        double E = out_params[3];
        double alpha = sqrt(2 * m * E / hbar_sqrd);
        double beta = sqrt(2 * m * (v0 - E) / hbar_sqrd);
        A = bisect_params(get_sq_integral, out_params, 0, 1e3, max_iter, TOL);
        out_params[4] = A;
        D = get_D(A, out_params);
        //printf("D: %3.8f\r\n", D);

        out_params[5] = D;

        L = out_params[1];
        //printf("L: %3.20f\r\n", L);
        if (true) {
            ofstream datafile;
            stringstream stream;
            stream << "dataplot_" << i << ".txt";
            string fileName = stream.str();
            datafile.open(fileName.c_str(), std::ios_base::trunc);
            datafile << "m,L,v0,E,A,D: ";
            for (int i = 0; i < 6; i++) {
                datafile << out_params[i] << ",";
            }
            datafile << "\r\n";
            for (double position = 0; position < L; position += step) {
                datafile << std::fixed << std::setprecision(8) << position;
                datafile << ",";
                datafile << std::fixed << std::setprecision(8) << wavefunction_01(position, out_params);
                datafile << "\r\n";
            }
            for (double position = L; position < 2 * L; position += step) {
                datafile << std::fixed << std::setprecision(8) << position;
                datafile << ",";
                datafile << std::fixed << std::setprecision(8) << wavefunction_02(position, out_params);
                datafile << "\r\n";
            }
            datafile.close();
        }
    }

    L = out_params[1];
    //printf("L: %3.20f\r\n", L);
    // we need to solve for an energy eigenvalue = potential.
    // i.e. let E = v0,  then find L.
    //bisect_params(energy_function_L, out_params, 0, 10, max_iter, TOL);
    // populate a list of graphs and count for each graphs the number of roots.
    // if the number hits 4 we're going to bisect that interval using the return function as the discrete number of counts
    // we can't use tolerance to stop iteration because many solutions will be exactly equal to 4, so we just need to keep iterating for a fixed number of attempts
    // OR, we can find the location of the root close to v0 and calculate its distance from v0 and use that as the convergence criterion
    double step_L = 0.1;
    double distance;
    // distance search
    for (double length = 2.72; length < 2.73; length += .001) {
        distance = get_distance(length, out_params, 5);
        //printf("distance: %3.8f\r\n", distance);
    }

    // potential search
    out_params[1] = L;
    double v0 = out_params[2];
    for (double pt = 1.6; pt < 2.7; pt += 0.01) {
        distance = get_vvp(pt, out_params, 5);
        //printf("distance: %3.8f\r\n", distance);
    }


}


// checks the hybrid newton method

void test01() {
    {
        double solution;
        double lowerBound = 0;
        double upperBound = 3;
        // hybrid method takes it more quickly to the solution
        solution = newtonhybrid(poly2test, poly2test_d, lowerBound, upperBound, max_iter, TOL, VERBOSE);
        solution = newtonsolve(poly2test, poly2test_d, (lowerBound + upperBound) / 2, max_iter, TOL);
        // hybrid method brackets the intended solution, newton finds a solution far away
        solution = newtonsolve(poly2test, poly2test_d, lowerBound, max_iter, TOL);
    }
}

// finds the energy eigenvalues for ebb1, the half-infinite well

int test02(double params[6], vector<double> & energies) {
    int num_eigenenergies = 0;
    FILE *datafile;
    datafile = fopen("data01.txt", "wt");
    double step_percent = 0.01;
    double potential = params[2];
    double step = potential * step_percent;
    int count = 0;
    double lowerBound = 0;
    vector<pair<double, double> > root_bracket;
    int num_of_brackets = 0;
    double root;
    double f_E;
    double f_E_prev = energy_function(0, params);
    for (double energy = 0; energy < potential; energy += step) {
        f_E = energy_function(energy, params);
        if (signbit(f_E) != signbit(f_E_prev)) {
            count++;
            lowerBound = energy - step;
            pair<double, double> bracket;
            bracket.first = lowerBound;
            bracket.second = energy;
            root_bracket.push_back(bracket);
            num_of_brackets = num_of_brackets + 1;
            printf("lower bound: %3.8f, upper bound: %3.8f\r\n", lowerBound, energy);
        }
        f_E_prev = f_E;
        if (VERBOSE) {
            printf("energy: %3.8f ; function: %3.8f\r\n", energy, f_E);
        }
        fprintf(datafile, "%3.8f,%3.8f\r\n", energy, energy_function(energy, params));
    }
    fclose(datafile);

    for (int i = 0; i < num_of_brackets; i++) {
        root = newtonhybridp(energy_function, energy_function_d, params, root_bracket.at(i).first, root_bracket.at(i).second, max_iter, TOL, VERBOSE);
        energies.push_back(root);
        num_eigenenergies = num_eigenenergies + 1;
    }
    return num_eigenenergies;
}

void test03(double* roots, int num_of_roots) {

}