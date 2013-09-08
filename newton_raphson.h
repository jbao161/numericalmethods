/* 
 * File:   newton_raphson.h
 * Author: jbao
 *
 * Created on August 31, 2013, 10:09 PM
 */

#ifndef NEWTON_RAPHSON_H
#define	NEWTON_RAPHSON_H

double newtonsolve(double (*function)(double), double (*fderivative)(double), double, int, double);


#endif	/* NEWTON_RAPHSON_H */

