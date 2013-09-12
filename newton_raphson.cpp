/*
 * File:   newton_raphson.cpp
 * Author: jbao
 *
 * Created on August 28, 2013, 4:10 PM
 */

#include <float.h>
#include <stdio.h>
#include <cmath>
#include "util_math.h"

using namespace std;
const bool DEBUG = true;

/**
 * Newton-Raphson method for finding real roots to a continuous function of a single variable.
 * @param function the equation we are trying to solve
 * @param fderivative the derivative of the equation we are trying to solve
 * @param initial_guess an approximation to the solution we are trying to find
 * @param max_iter maximum number of iterations allowed
 * @param TOL the maximum error term allowed for convergence
 * @return
 */
double newtonsolve(double (*function)(double), double (*fderivative)(double), double initial_guess, int max_iter, double TOL) {
    double f_x, derivative, next_guess;
    double convergence_criterion = DBL_MAX; // the relative distance between successive approximations
    if (DEBUG) {
        printf("\r\nNewton-Raphson search:\r\n");
        printf("initial guess: %3.2f\r\n", initial_guess);
    }
    // search for the solution using Newton-Raphson method
    for (int i = 0; i <= max_iter; i++) // i<= maxIter because results of nth iteration are not checked until the n+1th
    {
        f_x = function(initial_guess);
        // check if the function at the initial guess is sufficiently close to zero
        if (f_x == 0 || convergence_criterion < TOL) {
            if (DEBUG) {
                printf("after %i iterations\r\n", i); // if initial guess is the solution, i=0
                printf("solution: %3.12f\r\n", initial_guess);
                printf("f(%3.12f) = %3.12f\r\n", initial_guess, function(initial_guess));
            }
            return initial_guess; // found a solution!
        }
        // otherwise calculate the derivative at the initial guess.
        derivative = fderivative(initial_guess);
        if (derivative == 0) // verify that the derivative is nonzero
        {
            if (DEBUG) {
                printf("derivative at: %3.12f equals zero. Divide by zero error!\r\n", initial_guess);
            }
            return DBL_MAX; // cannot find a root due to divide by zero error!
        }
        // update the initial guess and convergence criterion
        next_guess = initial_guess - f_x / derivative;
        convergence_criterion = abs((next_guess - initial_guess) / next_guess);
        initial_guess = next_guess;
        if (DEBUG) {
            printf("%3.12f\r\n", next_guess);
        }
    }
    // unsuccessful search
    if (DEBUG) {
        printf("Method failed after %i iterations\r\n", max_iter);
        printf("inadequate solution: %3.12f with a tolerance of %3.12f\r\n", initial_guess, convergence_criterion);
        printf("f(%3.12f) = %3.12f\r\n", initial_guess, function(initial_guess));
    }
    return DBL_MAX; // cannot find a root because exceeded max iterations
}

/* comments:
 * 1. we solve f(x) = f(a) + (x-a) d/dx f(a) + error, for f(x) = 0
 * 2. the updated solution uses a subtraction, because f(a)/f'(a) +(x-a) = 0 gives x = a - f(a)/f'(a)
 * 3. linear approximation assumes higher order derivatives are small compared to f'(a)
 * 4. method fails if a derivative equals zero due to divide by zero error
 * 5. initial guess must be 'sufficiently close' to solution for convergence
 * 6. consecutive derivatives decrease to zero, so rate of convergence increases with each iteration
 * 7. convergence criterion is the change in x, divided by the approximation x
 * 8. approximations may oscillate between two values, and never converge!
 * 9. if the derivative dominates the function over the entire region, i.e. f(x)/df_dx(x) is close to zero, then convergence can be very slow!
 */

double newtonhybrid(double (*function)(double), double (*fderivative)(double), double lowerBound, double upperBound, int max_iter, double TOL, bool VERBOSE) {
    double precision; // half the length of the interval containing solution
    double functionBot = function(lowerBound);
    double functionTop = function(upperBound);
    double derivative, next_guess, midpoint;
    double initial_guess = (upperBound + lowerBound) / 2;
    double f_x = function(initial_guess);

    double convergence_criterion = DBL_MAX; // the relative distance between successive approximations

    // verify that inputs are in fact the lower and upper bound of the solution
    if (VERBOSE) {
        printf("\r\nHybrid search:\r\n");
        printf("initial guess: %3.2f\r\n", initial_guess);
        printf("f(%3.2f) = %3.2f\r\n", lowerBound, functionBot);
        printf("f(%3.2f) = %3.2f\r\n", upperBound, functionTop);
    }
    // function must be defined at the bounds, and its values must be of opposite sign
    if (isFiniteNumber(functionBot) == false || isFiniteNumber(functionTop) == false || signbit(functionBot) == signbit(functionTop)) {
        if (VERBOSE) {
            printf("invalid inputs: function must be defined at the bounds, and its values must be of opposite sign\r\n");
        }
        return DBL_MAX; // cannot find a root because of invalid inputs
    }

    // search for the solution using hybrid method
    for (int i = 0; i <= max_iter; i++) // i<= maxIter because results of nth iteration are not checked until the n+1th
    {
        precision = (upperBound - lowerBound) / 2;
        // check if the function at the initial guess is sufficiently close to zero
        if (abs(f_x) < TOL || convergence_criterion < TOL) {
            if (VERBOSE) {
                printf("after %i iterations\r\n", i); // if initial guess is the solution, i=0
                printf("solution: %3.12f with precision +/- %3.12f\r\n", initial_guess, precision);
                printf("f(%3.12f) = %3.12f\r\n", initial_guess, function(initial_guess));
            }
            return initial_guess; // found a solution!
        }
        // otherwise calculate the derivative at the initial guess.
        derivative = fderivative(initial_guess);
        if (derivative == 0) // verify that the derivative is nonzero
        {
            if (VERBOSE) {
                printf("derivative at: %3.12f equals zero. Divide by zero error!\r\n", initial_guess);
            }
            return DBL_MAX; // cannot find a root due to divide by zero error!
        }
        // determine the next approximation
        next_guess = initial_guess - f_x / derivative;

        // if the next approximation goes outside the bounds, take a bisection step
        if (!(next_guess > lowerBound & next_guess < upperBound)) {

            midpoint = (upperBound + lowerBound) / 2;
            f_x = function(midpoint);
            if (signbit(functionBot) == signbit(f_x)) {
                // solution is in the upper subinterval
                lowerBound = midpoint;
                functionBot = f_x;
            } else {
                // solution is in the lower subinterval
                upperBound = midpoint;
            }
            next_guess = (upperBound + lowerBound) / 2;
            f_x = function(next_guess);
            if (VERBOSE) {
                printf("bisection step: guess = %3.8f, low = %3.8f, high = %3.8f.\r\n", next_guess, lowerBound, upperBound);
            }
            if (VERBOSE) {
                printf("%3.12f (b)\r\n", next_guess);
            }
        }// if newton approximation is good, just update the bounds
        else {
            f_x = function(next_guess);
            if (signbit(functionBot) == signbit(f_x)) {
                // solution is in the upper subinterval
                lowerBound = next_guess;
                functionBot = f_x;
            } else {
                // solution is in the lower subinterval
                upperBound = next_guess;
            }
            if (VERBOSE) {
                printf("%3.12f (n)\r\n", next_guess);
            }
        }
        convergence_criterion = abs((next_guess - initial_guess) / next_guess);
        // set the initial guess for the next iteration
        initial_guess = next_guess;
    }
    // unsuccessful search
    if (VERBOSE) {
        printf("Method failed after %i iterations\r\n", max_iter);
        printf("inadequate solution: %3.12f with a tolerance of %3.12f\r\n", initial_guess, convergence_criterion);
        printf("f(%3.12f) = %3.12f\r\n", initial_guess, function(initial_guess));
    }
    return DBL_MAX; // cannot find a root because exceeded max iterations
}

double newtonhybridp(double (*function)(double, double*), double (*fderivative)(double, double*), double *params, double lowerBound, double upperBound, int max_iter, double TOL, bool VERBOSE) {
    double precision; // half the length of the interval containing solution
    double functionBot = function(lowerBound, params);
    double functionTop = function(upperBound, params);
    double derivative, next_guess, midpoint;
    double initial_guess = (upperBound + lowerBound) / 2;
    double f_x = function(initial_guess, params);

    double convergence_criterion = DBL_MAX; // the relative distance between successive approximations

    // verify that inputs are in fact the lower and upper bound of the solution
    if (VERBOSE) {
        printf("\r\nHybrid search:\r\n");
        printf("initial guess: %3.2f\r\n", initial_guess);
        printf("f(%3.2f) = %3.2f\r\n", lowerBound, functionBot);
        printf("f(%3.2f) = %3.2f\r\n", upperBound, functionTop);
    }
    // function must be defined at the bounds, and its values must be of opposite sign
    if (isFiniteNumber(functionBot) == false || isFiniteNumber(functionTop) == false || signbit(functionBot) == signbit(functionTop)) {
        if (VERBOSE) {
            printf("invalid inputs: function must be defined at the bounds, and its values must be of opposite sign\r\n");
        }
        return DBL_MAX; // cannot find a root because of invalid inputs
    }

    // search for the solution using hybrid method
    for (int i = 0; i <= max_iter; i++) // i<= maxIter because results of nth iteration are not checked until the n+1th
    {
        precision = (upperBound - lowerBound) / 2;
        // check if the function at the initial guess is sufficiently close to zero
        if (abs(f_x) < TOL || convergence_criterion < TOL) {
            if (VERBOSE) {
                printf("after %i iterations\r\n", i); // if initial guess is the solution, i=0
                printf("solution: %3.12f with precision +/- %3.12f\r\n", initial_guess, precision);
                printf("f(%3.12f) = %3.12f\r\n", initial_guess, function(initial_guess, params));
            }
            return initial_guess; // found a solution!
        }
        // otherwise calculate the derivative at the initial guess.
        derivative = fderivative(initial_guess, params);
        if (derivative == 0) // verify that the derivative is nonzero
        {
            if (VERBOSE) {
                printf("derivative at: %3.12f equals zero. Divide by zero error!\r\n", initial_guess);
            }
            return DBL_MAX; // cannot find a root due to divide by zero error!
        }
        // determine the next approximation
        next_guess = initial_guess - f_x / derivative;

        // if the next approximation goes outside the bounds, take a bisection step
        if (!(next_guess > lowerBound & next_guess < upperBound)) {

            midpoint = (upperBound + lowerBound) / 2;
            f_x = function(midpoint, params);
            if (signbit(functionBot) == signbit(f_x)) {
                // solution is in the upper subinterval
                lowerBound = midpoint;
                functionBot = f_x;
            } else {
                // solution is in the lower subinterval
                upperBound = midpoint;
            }
            next_guess = (upperBound + lowerBound) / 2;
            f_x = function(next_guess, params);
            if (VERBOSE) {
                printf("bisection step: guess = %3.8f, low = %3.8f, high = %3.8f.\r\n", next_guess, lowerBound, upperBound);
            }
            if (VERBOSE) {
                printf("%3.12f (b)\r\n", next_guess);
            }
        }// if newton approximation is good, just update the bounds
        else {
            f_x = function(next_guess, params);
            if (signbit(functionBot) == signbit(f_x)) {
                // solution is in the upper subinterval
                lowerBound = next_guess;
                functionBot = f_x;
            } else {
                // solution is in the lower subinterval
                upperBound = next_guess;
            }
            if (VERBOSE) {
                printf("%3.12f (n)\r\n", next_guess);
            }
        }
        convergence_criterion = abs((next_guess - initial_guess) / next_guess);
        // set the initial guess for the next iteration
        initial_guess = next_guess;
    }
    // unsuccessful search
    if (VERBOSE) {
        printf("Method failed after %i iterations\r\n", max_iter);
        printf("inadequate solution: %3.12f with a tolerance of %3.12f\r\n", initial_guess, convergence_criterion);
        printf("f(%3.12f) = %3.12f\r\n", initial_guess, function(initial_guess, params));
    }
    return DBL_MAX; // cannot find a root because exceeded max iterations
}

