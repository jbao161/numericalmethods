
#include <float.h>

bool isFiniteNumber(double x) {
    return (x <= DBL_MAX && x >= -DBL_MAX);
}
