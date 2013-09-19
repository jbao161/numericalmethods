/* 
 * File:   functions02.h
 * Author: jbao
 *
 * Created on September 12, 2013, 3:32 AM
 */

#ifndef FUNCTIONS02_H
#define	FUNCTIONS02_H

double poly2test(double);
double poly2test_d(double);
double energy_function(double, double*);
double energy_function_d(double, double*);
void convert_inputs(double[][2], double*);
double wavefunction_01(double, double*);
double wavefunction_02(double, double*);
double get_sq_integral(double, double*);
double get_D(double, double*);
double energy_function_L(double, double*);
double energy_function_v0(double, double*);
double get_distance(double, double*, int);

#endif	/* FUNCTIONS02_H */

