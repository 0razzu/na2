#include "fredholm.hpp"


void create_matrixes(double lambda, double (* ker)(double, double), double (* f)(double), double** a, double* f_n, unsigned n) {
    const double h = 1. / n;
    
    for (unsigned i = 0; i < n; i++) {
        double s_i = (double(i) + 0.5) * h;
        
        for (unsigned j = 0; j < n; j++) {
            a[i][j] = -h * lambda * ker(s_i, (double(j) + 0.5) * h);
            
            if (i == j)
                a[i][j]++;
        }
        
        f_n[i] = f(s_i);
    }
}


double compute_accuracy(double* x_prev, double* x, unsigned n) {
    double max_dif = 0;
    
    for (unsigned i = 0; i < n; i++) {
        double cur_dif = abs(x[i] - x_prev[i]);
        
        if (cur_dif > max_dif)
            max_dif = cur_dif;
    }
    
    return max_dif;
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
