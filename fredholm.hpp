#ifndef fredholm_hpp
#define fredholm_hpp

#include <iostream>


void create_matrixes(double lambda, double (* ker)(double, double), double (* f)(double), double** a, double* f_i, unsigned n);

double solve(double** a, double* f, unsigned n);


#endif
