/* 
 * File:   polynomial.h
 * Author: jbao
 *
 * Created on August 28, 2013, 5:01 PM
 */
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>
using namespace std;

class Polynomial {
public:
    map<double, double> polyCoefs;
    void set(double*, int);
    string toText(string, int);
    void print();
    double evaluate(double);
};
#endif


