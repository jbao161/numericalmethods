
/*
 * File:   polynomial.cpp
 * Author: jbao
 *
 * Created on August 28, 2013, 5:01 PM
 */
#include "polynomial.h"

void Polynomial::set(double *coefs, int size)
{
    for (int i = 0; i < size; i++)
    {
        polyCoefs[i] = coefs[i];
    }
}

double Polynomial::evaluate(double x)
{
    double result = 0;
    typedef std::map<double, double>::iterator polyIter;
    for (polyIter i = polyCoefs.begin(); i != polyCoefs.end(); i++)
    {
        double coef = (*i).second;
        double exponent = (*i).first;
        result += coef * pow(x, exponent);
    }
    return result;
}

void Polynomial::print()
{
    cout << toText("x", 2) << endl;
}

string Polynomial::toText(string varName, int decimals)
{
    string result = "";
    string add = " + ";
    string subtract = " - ";
    string exp = "^";
    string operation = add;
    typedef std::map<double, double>::iterator polyIter;
    for (polyIter i = polyCoefs.begin(); i != polyCoefs.end(); i++)
    {
        double coef = (*i).second;
        double exponent = (*i).first;
        // if the coefficient is negative, and the term is not a constant, use subtract string and positive coefficient
        if (coef < 0 & exponent != 0)
        {
            operation = subtract;
            coef = -coef;
        }
        else
        {
            operation = add;
        }

        // don't print the term if its coefficient equals zero
        if (coef != 0)
        {
            result += operation;
            // print coefficient if term is constant, or coefficient does not equal one
            if (exponent == 0 || coef != 1)
            {

                ostringstream stringCoef; // stream used for the conversion
                stringCoef << std::fixed << std::setprecision(decimals) << coef;
                result += stringCoef.str();
            }
            // print the variable name if term is not a constant
            if (exponent != 0)
            {
                result += varName;
            }
            // if the term is not a constant and its exponent is not one, we need to show the exponent
            if (exponent != 0 & exponent != 1)
            {
                ostringstream stringExponent; // stream used for the conversion
                stringExponent << exponent;
                result += exp + stringExponent.str();
            }
        }
    }
    // formats the leading operation string
    string leadingOperator = result.substr(0, operation.length());
    if (leadingOperator == add)
    {
        result = result.substr(add.length());
    }
    else if (leadingOperator == subtract)
    {
        result = "-" + result.substr(subtract.length());
    }
    else if (result == "")
    {
        return "0";
    }
    return result;
}
