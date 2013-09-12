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
    return alpha * cos(alpha_L) - beta * sin(alpha_L);
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