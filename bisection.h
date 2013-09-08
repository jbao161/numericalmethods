/* 
 * File:   bisection.h
 * Author: jbao
 *
 * Created on August 31, 2013, 7:58 PM
 */

#ifndef BISECTION_H
#define	BISECTION_H

double bisect(double (*function)(double), double, double, int, double);
double bisect_params(double (*function)(double, double*), double*, double, double, int, double);
double bisect_r(double (*function)(double), double, double, int,int, double);

#endif	/* BISECTION_H */

