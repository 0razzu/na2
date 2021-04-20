#include "fredholm.hpp"


void create_matrixes(double lambda, double (* ker)(double, double), double (* f)(double), double** a, double* f_i, unsigned n) {
    const double h = 1. / n;
    
    for (unsigned i = 0; i < n; i++) {
        double s_i = (double(i) + 0.5) * h;
        
        for (unsigned j = 0; j < n; j++) {
            a[i][j] = -h * lambda * ker(s_i, (double(j) + 0.5) * h);
            
            if (i == j)
                a[i][j]++;
        }
        
        f_i[i] = f(s_i);
    }
}


void revert_tr(double** u, unsigned n) {
    for (unsigned i = 0; i < n; i++)
        u[i][i] = 1. / u[i][i];
    
    for (unsigned i = 0; i < n; i++)
        for (unsigned j = 0; j < i; j++) {
            double sum = 0;
            for (unsigned k = j; k < i; k++)
                sum += u[i][k] * u[k][j];
            
            u[i][j] = -u[i][i] * sum;
    }
}


double solve(double** a, double* f, unsigned n) {
    return n;
}
