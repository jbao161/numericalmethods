#include <math.h>
#include <stdio.h>
const double hbar_sqrd = 0.076199682; // me ev nm^2
const bool DEBUG = true;
// 6.62606957e-27 erg * second = g* cm^2 / s

//6.62606957 
// 2 nm

// hbar = 6.62606957 × 10-34 J*s
// 1 joule = 1e7 ergs
// erg = g·cm^2/s^2
// angstrom = 1e-8 cm;

double poly2test(double x) {
    return x * x - 2 * x - 2;
}

double poly2test_d(double x) {
    return 2 * x - 2;
}

double convert_mass(double gram, double exponent) {
    gram *= 1.0977693108;
    exponent += 27;
    double electron_mass = gram * pow(10, exponent);
    return electron_mass;
}

double convert_energy(double erg, double exponent) {
    erg *= 6.24150934;
    exponent += 11;
    double electron_volt = erg * pow(10, exponent);
    return electron_volt;
}

double convert_length(double angstrom, double exponent) {
    exponent -= 1;
    double nanometer = angstrom * pow(10, exponent);
    return nanometer;
}
// convert the units gram, angstrom, and ergs into electron mass, centimeter, and electron volt

void convert_inputs(double in_params[][2], double* out_params) {
    out_params[0] = convert_mass(in_params[0][0], in_params[0][1]);
    out_params[1] = convert_length(in_params[1][0], in_params[1][1]);
    out_params[2] = convert_energy(in_params[2][0], in_params[2][1]);
    if (DEBUG) {
        for (int i = 0; i < 3; i++) {
            printf("%3.8f\r\n", out_params[i]);
        }
    }

}

double energy_function(double E, double *params) {
    double m = params[0]; // in electron mass
    double L = params[1]; // in nanometer
    double v0 = params[2]; // in electron volt    
    double alpha = sqrt(2 * m * E / hbar_sqrd);
    double beta = sqrt(2 * m * (v0 - E) / hbar_sqrd);
    double alpha_L = alpha * L;
    return alpha * cos(alpha_L) + beta * sin(alpha_L);
    //sin(Lk) + sqrt(E / (v0 - E)) * cos(Lk);
}

double energy_function_L(double L, double *params) {
    double m = params[0]; // in electron mass
    double E = params[3]; // in nanometer
    double v0 = params[2]; // in electron volt    
    double alpha = sqrt(2 * m * E / hbar_sqrd); // this quantity is known
    double beta = sqrt(2 * m * (v0 - E) / hbar_sqrd);
    //double beta = 0; // since v0 = e
    double alpha_L = alpha * L;
    // need to solve 0 = cos (alpha * L) for L
    return alpha * cos(alpha_L) + beta * sin(alpha_L);
    //sin(Lk) + sqrt(E / (v0 - E)) * cos(Lk);
}

double energy_function_v0(double v0, double *params) {
    double m = params[0]; // in electron mass
    double E = params[3]; // in nanometer
    double L = params[1]; // in electron volt    
    double alpha = sqrt(2 * m * E / hbar_sqrd);
    double beta = sqrt(2 * m * (v0 - E) / hbar_sqrd);
    double alpha_L = alpha * L;
    return alpha * cos(alpha_L) + beta * sin(alpha_L);
    //sin(Lk) + sqrt(E / (v0 - E)) * cos(Lk);
}

double energy_function_d(double E, double *params) {
    double m = params[0];
    double L = params[1];
    double v0 = params[2];
    double alpha_L = L * sqrt(2 * m * E / hbar_sqrd);
    double dd_alpha = sqrt(2 * m / hbar_sqrd) / 2 * pow(E, -0.5);
    double dd_coef = 0.5 * pow((E / (v0 - E)), -0.5) * v0 / pow(v0 - E, 2);
    return dd_alpha * (L * cos(alpha_L) + sqrt(E / (v0 - E)) * -L * sin(alpha_L)) + dd_coef * cos(alpha_L);
}

double wavefunction_01(double x, double *params) {
    double m = params[0];
    double L = params[1];
    double E = params[3];
    double A = params[4];
    double alpha = sqrt(2 * m * E / hbar_sqrd);
    return A * sin(alpha * x);
}

double wavefunction_02(double x, double *params) {
    double m = params[0];
    double v0 = params[2];
    double E = params[3];
    double D = params[5];
    double beta = sqrt(2 * m * (v0 - E) / hbar_sqrd);
    return D * exp(-beta * x);
    ;
}

double get_sq_integral(double A, double *params) {
    double m = params[0];
    double L = params[1];
    double v0 = params[2];
    double E = params[3];
    double alpha = sqrt(2 * m * E / hbar_sqrd);
    double beta = sqrt(2 * m * (v0 - E) / hbar_sqrd);
    double D = A * exp(beta * L) * sin(alpha * L);
    double integral_01 = A * A * (0.5 * L + fabs(sin(2 * alpha * L))*0.25 / alpha);
    double integral_02 = fabs(D * D * 0.5 / beta * exp(-2 * beta * L));
    return -1 + integral_01 + integral_02;
}

double get_D(double A, double *params) {
    double m = params[0];
    double L = params[1];
    double v0 = params[2];
    double E = params[3];
    double alpha = sqrt(2 * m * E / hbar_sqrd);
    double beta = sqrt(2 * m * (v0 - E) / hbar_sqrd);
    double D = A * exp(beta * L) * sin(alpha * L);

    return D;
}

double get_distance(double L, double *params, int num_energies) {
    // for a given L, start plotting the energy function until the specified number of energy roots are found
    {
        int num_eigenenergies = 0;

        double step_percent = 0.01;
        double potential = params[2];
        double step = potential * step_percent;
        int count = 0;
        double lowerBound = 0;
        int num_of_brackets = 0;
        double bracket_first, bracket_second;
        double root;
        double f_E;
        double f_E_prev = energy_function(0, params);
        for (double energy = step; energy < potential; energy += step) {
            f_E = energy_function(energy, params);
            if (signbit(f_E) != signbit(f_E_prev)) {
                count++;
                lowerBound = energy - step;
                bracket_first = lowerBound;
                bracket_second = energy;
                num_of_brackets = num_of_brackets + 1;
                // printf("num brackets: %d\r\n", num_of_brackets);
            }
            f_E_prev = f_E;
            if (false) {
                printf("energy: %3.8f ; function: %3.8f\r\n", energy, f_E);
            }
        }
        // bisection search with each increment step that changes sign
    }
    // if there aren't at least that many found, return null = double.min
    // if there are more than that many, return double.max
    // return the distance of the energy root from v_0 the potential
}
