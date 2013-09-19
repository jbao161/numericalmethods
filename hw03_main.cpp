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
#include <float.h>

using namespace std;
const int max_iter = 1e4;
const double TOL = 5e-16;
const bool VERBOSE = false;
const double hbar_sqrd = 0.076199682; // me ev nm^2
void test01();
int test02(double*, vector<double>&);
void test03();
void test04(int);

void hw03_main() {
    if (false) {
        test01();
    }
    if (true) {
        test03();
    }
    if (false) {
        test04(4);
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
    fprintf(datafile, "length: %3.8f, potential: %3.8f", params[0], potential);
    for (double energy = step; energy < potential; energy += step) {
        f_E = energy_function(energy, params);
        if (signbit(f_E) != signbit(f_E_prev)) {
            count++;
            lowerBound = energy - step;
            pair<double, double> bracket;
            bracket.first = lowerBound;
            bracket.second = energy;
            root_bracket.push_back(bracket);
            num_of_brackets = num_of_brackets + 1;
            // printf("num brackets: %d\r\n", num_of_brackets);
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

// calls test02 to find the energy eigenvalues, then plots each eigenvalue's wavefunction to a file

void test03() {
    double in_params[3][2] = {
        {1, -27},
        {20, 0},
        {15, -13}
    };
    double out_params[6];
    convert_inputs(in_params, out_params);
    vector<double> energies;
    int num_eigenenergies = test02(out_params, energies); // put the eigenvalues in outparams
    double A = 1; // wavefunction_01 coefficient
    double D = 1; // wavefunction_02 coefficient
    double L; // length of the well
    // double v0 = out_params[2]; // potential of well
    double step = 0.01;
    for (int i = 0; i < num_eigenenergies; i++) {
        out_params[3] = energies.at(i);
        A = bisect_params(get_sq_integral, out_params, 0, 1e3, max_iter, TOL);
        double m = out_params[0];
        double L = out_params[1];
        double v0 = out_params[2];
        double E = out_params[3];
        double alpha = sqrt(2 * m * E / hbar_sqrd);
        double beta = sqrt(2 * m * (v0 - E) / hbar_sqrd);
        //A = sqrt(0.5 * L - 0.25 / alpha * sin(2 * alpha * L) + 0.5 / beta * pow(sin(alpha * L), 2)); vbcfg bvn

        D = get_D(A, out_params);
        out_params[4] = A;
        out_params[5] = D;
        L = out_params[1];
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
    if (true) {
        // we need to solve for an energy eigenvalue = potential.
        // i.e. let E = v0,  then find L.
        //bisect_params(energy_function_L, out_params, 0, 10, max_iter, TOL);
        // populate a list of graphs and count for each graphs the number of roots.
        // if the number hits 4 we're going to bisect that interval using the return function as the discrete number of counts
        // we can't use tolerance to stop iteration because many solutions will be exactly equal to 4, so we just need to keep iterating for a fixed number of attempts
        // OR, we can find the location of the root close to v0 and calculate its distance from v0 and use that as the convergence criterion
        double step_L = 0.1;
        for (double length = 0; length < L; length += step_L) {

        }
    }
}

void do_test03(double L_input) {
    double in_params[3][2] = {
        {1, -27},
        {20, 0},
        {15, -13}
    };
    double out_params[6];
    out_params[0] = L_input;
    convert_inputs(in_params, out_params);
    vector<double> energies;
    int num_eigenenergies = test02(out_params, energies); // put the eigenvalues in outparams
    double A = 1; // wavefunction_01 coefficient
    double D = 1; // wavefunction_02 coefficient
    double L = L_input; // length of the well
    // double v0 = out_params[2]; // potential of well
    double step = 0.01;
    for (int i = 0; i < num_eigenenergies; i++) {
        out_params[3] = energies.at(i);
        A = bisect_params(get_sq_integral, out_params, 0, 1e3, max_iter, TOL);
        double m = out_params[0];
        double L = out_params[1];
        double v0 = out_params[2];
        double E = out_params[3];
        double alpha = sqrt(2 * m * E / hbar_sqrd);
        double beta = sqrt(2 * m * (v0 - E) / hbar_sqrd);
        //A = sqrt(0.5 * L - 0.25 / alpha * sin(2 * alpha * L) + 0.5 / beta * pow(sin(alpha * L), 2)); vbcfg bvn

        D = get_D(A, out_params);
        out_params[4] = A;
        out_params[5] = D;
        L = out_params[1];
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
    if (true) {
        // we need to solve for an energy eigenvalue = potential.
        // i.e. let E = v0,  then find L.
        //bisect_params(energy_function_L, out_params, 0, 10, max_iter, TOL);
        // populate a list of graphs and count for each graphs the number of roots.
        // if the number hits 4 we're going to bisect that interval using the return function as the discrete number of counts
        // we can't use tolerance to stop iteration because many solutions will be exactly equal to 4, so we just need to keep iterating for a fixed number of attempts
        // OR, we can find the location of the root close to v0 and calculate its distance from v0 and use that as the convergence criterion
        double step_L = 0.1;
        for (double length = 0; length < L; length += step_L) {

        }
    }
}

int find_roots(double params[6], vector<double> & energies) {
    int num_eigenenergies = 0;
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
    for (double energy = step; energy < potential; energy += step) {
        f_E = energy_function(energy, params);
        if (signbit(f_E) != signbit(f_E_prev)) {
            count++;
            lowerBound = energy - step;
            pair<double, double> bracket;
            bracket.first = lowerBound;
            bracket.second = energy;
            root_bracket.push_back(bracket);
            num_of_brackets = num_of_brackets + 1;
            // printf("num brackets: %d\r\n", num_of_brackets);
        }
        f_E_prev = f_E;
        if (VERBOSE) {
            printf("energy: %3.8f ; function: %3.8f\r\n", energy, f_E);
        }
    }


    for (int i = 0; i < num_of_brackets; i++) {
        root = newtonhybridp(energy_function, energy_function_d, params, root_bracket.at(i).first, root_bracket.at(i).second, max_iter, TOL, VERBOSE);
        energies.push_back(root);
        num_eigenenergies = num_eigenenergies + 1;
    }
    return num_eigenenergies;
}
// finds value of L necessary to have 2 more bound states

void test04(int num_of_roots) {
    // convert the input mass, length, and potential into the right units
    double in_params[3][2] = {
        {1, -27},
        {20, 0},
        {15, -13}
    };
    double out_params[6];
    convert_inputs(in_params, out_params);

    // plot the energyfunction as a function of E, for each L in increments from L and up
    double L = out_params[1];
    double v_0 = out_params[2];
    ///double E = out_params[3];
    double increment = L * 1e-3;
    double increment_e = v_0 * 1e-3;
    double function_e;
    double function_e_prev;
    double root;
    int root_count;
    double lowerBound, upperBound;
    int exceeded = 0;
    do_test03(L);
    for (double length = L, i = 0; length < 2 * L || i < 2; length += increment, i++) {

    }
}
