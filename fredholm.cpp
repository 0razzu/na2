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


double compute_accuracy(double* x_prev, double* x, unsigned n) {
    double x_prev_max = 0, x_max = 0;
    
    for (unsigned i = 0; i < n; i++) {
        double x_prev_cur = fabs(x_prev[i]);
        double x_cur = fabs(x[i]);
        
        if (x_prev_cur > x_prev_max)
            x_prev_max = x_prev_cur;
        
        if (x_cur > x_max)
            x_max = x_cur;
    }
    
    return fabs(x_max - x_prev_max);
}


int solve(double** a, double* f, double* x, unsigned n, double eps, unsigned max_it, unsigned& it, double& accuracy) {
    it = 0;
    
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < i; j++)
            a[i][j] /= a[i][i];
        for (unsigned j = i + 1; j < n; j++)
            a[i][j] /= a[i][i];
        
        f[i] /= a[i][i];
    }
    
    double* x_prev = new double[n];
    
    do {
        for (unsigned i = 0; i < n; i++)
            x_prev[i] = x[i];
        
        for (unsigned i = 0; i < n; i++) {
            x[i] = f[i];
            
            for (unsigned j = 0; j < i; j++)
                x[i] -= a[i][j] * x[j];
            for (unsigned j = i + 1; j < n; j++)
                x[i] -= a[i][j] * x_prev[j];
        }
        
        it++;
    } while (eps < (accuracy = compute_accuracy(x_prev, x, n)) && it < max_it);
    
    delete[] x_prev;
    
    if (eps < accuracy)
        return ACC_NOT_ACHVD;
    
    return OK;
}
