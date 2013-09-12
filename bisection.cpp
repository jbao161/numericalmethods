
/*
 * File:   bisection.cpp
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
 * Bisection method for finding real roots to a continuous function of a single variable.
 * @param function the equation we are trying to solve
 * @param lowerBound a lower bound to the root we are trying to find
 * @param upperBound an upper bound to the root we are trying to find
 * @param max_iter maximum number of iterations allowed
 * @param TOL the maximum error term allowed for convergence
 * @return
 */
double bisect(double (*function)(double), double lowerBound, double upperBound, int max_iter, double TOL)
{
    double midpoint, prevMidpoint;
    double precision; // length of the interval containing solution
    double convergence_criterion; // size of interval divided by the midpoint
    double functionMidpoint;
    double functionBot = function(lowerBound);
    double functionTop = function(upperBound);

    // verify that inputs are in fact the lower and upper bound of the solution
    if (DEBUG)
    {
        printf("\r\nBisection search:\r\n");
        printf("f(%3.2f) = %3.2f\r\n", lowerBound, functionBot);
        printf("f(%3.2f) = %3.2f\r\n", upperBound, functionTop);
    }
    // function must be defined at the bounds, and its values must be of opposite sign
    if (isFiniteNumber(functionBot) == false || isFiniteNumber(functionTop) == false || signbit(functionBot) == signbit(functionTop))
    {
        if (DEBUG)
        {
            printf("invalid inputs: function must be defined at the bounds, and its values must be of opposite sign\r\n");
        }
        return DBL_MAX; // cannot find a root because of invalid inputs
    }
    // search for the solution using bisection method
    for (int i = 0; i < max_iter;)
    {
        // compute midpoint and f(midpoint)
        midpoint = (upperBound + lowerBound) / 2;
        functionMidpoint = function(midpoint);

        // determine how precisely the midpoint approximates the solution
        precision = abs(upperBound - lowerBound);
        convergence_criterion = abs(precision / midpoint);

        // check if the function at the midpoint is sufficiently close to zero
        if (functionMidpoint == 0 || convergence_criterion < TOL)
        {
            if (DEBUG)
            {
                printf("after %i iterations\r\n", i + 1);
                printf("solution: %3.12f, with precision +/- %3.12f\r\n", midpoint, precision / 2);
                printf("f(%3.12f) = %3.12f\r\n", midpoint, function(midpoint));
            }
            return midpoint; // found a suitable solution!
        }

        // otherwise continue searching
        i++;
        prevMidpoint = midpoint;
        // bisect the interval and determine which half contains the solution
        if (signbit(functionBot) == signbit(functionMidpoint))
        {
            // use the upper subinterval
            lowerBound = midpoint;
            functionBot = functionMidpoint;
        }
        else
        {
            // use the lower subinterval
            upperBound = midpoint;
        }
    }
    // unsuccessful search
    if (DEBUG)
    {
        printf("Method failed after %i iterations\r\n", max_iter);
        printf("inadequate solution: %3.12f within an interval of length %3.12f\r\n", midpoint, precision);
        printf("f(%3.12f) = %3.12f\r\n", midpoint, function(midpoint));
    }
    return DBL_MAX; // cannot find a root because exceeded max iterations
}

/*
 * comments:
 * 1. convergence criterion is size of interval containing solution divide by approximate solution
 * 2. this is basically an application of the mean value theorem using some continuous function and exploiting the fact that a solution to f = 0 will exist between a and b, provided that f(a) < 0 and f(b) > 0 both exist.
 * 3. we need to solve f=0, because the algorithm relies on positive/negative sign checking to determine which interval contains the solution
 * 4. we need to maintain the lower bound of the interval to map to a negative and the upper bound of the interval to map to a positive to ensure a zero exists within the interval
 * 5. it doesn't really matter where we split the interval, we use the midpoint in the hope that halving the interval will find a solution the fastest, but the algorithm works for any point in the interval
 * 6. if we wanted to allow f = constant instead of f = 0, we could simply use an input function that is subtracted by the constant.
 * 7. if there are multiple solutions within the interval, the algorithm still works. it will output whichever solution it stumbles upon and won't know that other solutions exist
 * 8. cannot find a solution if the function only touches the line y=0, but doesn't cross!
 */

double bisect_r(double (*function)(double), double lowerBound, double upperBound, int current_iter, int max_iter, double TOL)
{
    // unsuccessful search
    if (current_iter > max_iter)
    {
        if (DEBUG)
        {
            printf("Method failed after %i iterations\r\n", max_iter);
        }
        return DBL_MAX; // cannot find a root because exceeded max iterations
    }


    double midpoint;
    double precision; // length of the interval containing solution
    double convergence_criterion; // size of interval divided by the midpoint
    double functionMidpoint;
    double functionBot = function(lowerBound);
    double functionTop = function(upperBound);


    // function must be defined at the bounds, and its values must be of opposite sign
    if (isFiniteNumber(functionBot) == false || isFiniteNumber(functionTop) == false || signbit(functionBot) == signbit(functionTop))
    {
        if (DEBUG)
        {
            printf("invalid inputs: function must be defined at the bounds, and its values must be of opposite sign\r\n");
        }
        return DBL_MAX; // cannot find a root because of invalid inputs
    }
    // search for the solution using bisection method
    // compute midpoint and f(midpoint)
    midpoint = (upperBound + lowerBound) / 2;
    functionMidpoint = function(midpoint);

    // determine how precisely the midpoint approximates the solution
    precision = abs(upperBound - lowerBound);
    convergence_criterion = abs(precision / midpoint);

    // check if the function at the midpoint is sufficiently close to zero
    if (functionMidpoint == 0 || convergence_criterion < TOL)
    {
        if (DEBUG)
        {
            printf("after %i iterations\r\n", current_iter);
            printf("solution: %3.12f, with precision +/- %3.12f\r\n", midpoint, precision / 2);
            printf("f(%3.12f) = %3.12f\r\n", midpoint, function(midpoint));
        }
        return midpoint; // found a suitable solution!
    }

    // otherwise continue searching
    // bisect the interval and determine which half contains the solution
    if (signbit(functionBot) == signbit(functionMidpoint))
    {
        // use the upper subinterval
        return bisect_r(function, midpoint, upperBound, current_iter + 1, max_iter, TOL);
    }
    else
    {
        // use the lower subinterval
        return bisect_r(function, lowerBound, midpoint, current_iter + 1, max_iter, TOL);
    }
}

double bisect_params(double (*function)(double, double*), double *params, double lowerBound, double upperBound, int max_iter, double TOL)
{
    double midpoint, prevMidpoint;
    double precision; // length of the interval containing solution
    double convergence_criterion; // size of interval divided by the midpoint
    double functionMidpoint;
    double functionBot = function(lowerBound,params);
    double functionTop = function(upperBound,params);

    // verify that inputs are in fact the lower and upper bound of the solution
    if (DEBUG)
    {
        printf("\r\nBisection search:\r\n");
        printf("f(%3.2f) = %3.2f\r\n", lowerBound, functionBot);
        printf("f(%3.2f) = %3.2f\r\n", upperBound, functionTop);
    }
    // function must be defined at the bounds, and its values must be of opposite sign
    if (isFiniteNumber(functionBot) == false || isFiniteNumber(functionTop) == false || signbit(functionBot) == signbit(functionTop))
    {
        if (DEBUG)
        {
            printf("invalid inputs: function must be defined at the bounds, and its values must be of opposite sign\r\n");
        }
        return DBL_MAX; // cannot find a root because of invalid inputs
    }
    // search for the solution using bisection method
    for (int i = 0; i < max_iter;)
    {
        // compute midpoint and f(midpoint)
        midpoint = (upperBound + lowerBound) / 2;
        functionMidpoint = function(midpoint,params);

        // determine how precisely the midpoint approximates the solution
        precision = abs(upperBound - lowerBound);
        convergence_criterion = abs(precision / midpoint);

        // check if the function at the midpoint is sufficiently close to zero
        if (functionMidpoint == 0 || convergence_criterion < TOL)
        {
            if (DEBUG)
            {
                printf("after %i iterations\r\n", i + 1);
                printf("solution: %3.12f, with precision +/- %3.12f\r\n", midpoint, precision / 2);
                printf("f(%3.12f) = %3.12f\r\n", midpoint, function(midpoint,params));
            }
            return midpoint; // found a suitable solution!
        }

        // otherwise continue searching
        i++;
        prevMidpoint = midpoint;
        // bisect the interval and determine which half contains the solution
        if (signbit(functionBot) == signbit(functionMidpoint))
        {
            // use the upper subinterval
            lowerBound = midpoint;
            functionBot = functionMidpoint;
        }
        else
        {
            // use the lower subinterval
            upperBound = midpoint;
        }
    }
    // unsuccessful search
    if (DEBUG)
    {
        printf("Method failed after %i iterations\r\n", max_iter);
        printf("inadequate solution: %3.12f within an interval of length %3.12f\r\n", midpoint, precision);
        printf("f(%3.12f) = %3.12f\r\n", midpoint, function(midpoint,params));
    }
    return DBL_MAX; // cannot find a root because exceeded max iterations
}
