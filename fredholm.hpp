#ifndef fredholm_hpp
#define fredholm_hpp

#include <iostream>


void create_matrixes(double lambda, double (* ker)(double, double), double (* f)(double), double** a, double* f_i, unsigned n);

int solve(double** a, double* f, double* x, unsigned n, double eps, unsigned max_it, unsigned& it);


#endif
