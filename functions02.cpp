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
    exponent -= 8;
    double centimeter = angstrom * pow(10, exponent);
    return centimeter;
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

void convert_params_in(double *params) {
    params[0] *= 1e-23; // from 1e-27 g to 1e-4 in (1e-23 of a gram))
    params[1] *= 1e-2; // from 2e-8 cm to 1e-6 (in .01 of cm)
    params[2] *= 1e-13; // from 1e-13 to 1e0 (in 1e-13 of an erg)
}

void convert_params_out(double *params) {
    params[0] *= 1e23; // to 1e-27 g from 1e-4 in (1e-23 of a gram))
    params[1] *= 1e2; // to 2e-8 cm from 1e-6 (in .01 of cm)
    params[2] *= 1e13; // to 1e-13 from 1e0 (in 1e-13 of an erg)
}

double energy_function(double E, double *params) {
    double m = params[0]; // in (1e-23 of a gram)
    double L = params[1]; // in (.01 of a cm)
    double v0 = params[2]; // (in 1e-13 of an erg)
    double Lk = L * sqrt(2 * m * E / hbar_sqrd);
    double k = sqrt(2 * m * E / hbar_sqrd);
    double lambda = sqrt(2 * m * (v0 - E) / hbar_sqrd);
    return k * cos(Lk) - lambda * sin(Lk);
    //sin(Lk) + sqrt(E / (v0 - E)) * cos(Lk);
}

double energy_function_d(double E, double *params) {
    double m = params[0];
    double L = params[1];
    double v0 = params[2];
    double Lk = L * sqrt(2 * m * E / hbar_sqrd);
    double ddk = sqrt(2 * m / hbar_sqrd) / 2 * pow(E, -0.5);
    return L * cos(Lk) + sqrt(E / (v0 - E)) * cos(Lk);
}