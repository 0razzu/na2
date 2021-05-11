#include "fredholm.hpp"
#include "invert.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>


const unsigned WIDTH = 10;
const unsigned PRECISION = 7;


double ker(double t, double s) {
    return cos(t) + cos(s);
}


double f(double s) {
    return -sin(1) - cos(1) - .5 * cos(s) + 1 - s;
}


double phi_exact(double s) {
    return -s;
}


double matrix_norm(double** a, unsigned n) {
    double norm = 0;
    
    for (unsigned i = 0; i < n; i++) {
        double curr_str_len = 0;
        
        for (unsigned j = 0; j < n; j++)
            curr_str_len += fabs(a[i][j]);
        
        if (curr_str_len > norm)
            norm = curr_str_len;
    }
    
    return norm;
}


double inv_matrix_norm(double** a, unsigned n) {
    double** a_inv = new double*[n];
    for (unsigned i = 0; i < n; i++) {
        a_inv[i] = new double[n];
        
        for (unsigned j = 0; j < n; j++)
            a_inv[i][j] = a[i][j];
    }
    invert(a_inv, n);
    
    double norm = matrix_norm(a_inv, n);
    
    for (unsigned i = 0; i < n; i++)
        delete[] a_inv[i];
    delete[] a_inv;
    
    return norm;
}


void print(double** a, double* f_i, unsigned n) {
    std::cout.setf(std::ios::right | std::ios::fixed | std::ios::showpoint);
    std::cout.precision(PRECISION);
    
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++)
            std::cout << std::setw(WIDTH) << a[i][j] << '\t';
        std::cout << '\t' << std::setw(WIDTH) << f_i[i] << std::endl;
    }
    
    std::cout << std::endl;
}


void print(double* x, unsigned n) {
    std::cout.flags(std::ios::left | std::ios::fixed | std::ios::showpoint);
    std::cout << std::setw(WIDTH) << "Computed" << '\t';
    std::cout << std::setw(WIDTH) << "Exact" << '\t';
    std::cout << std::setw(WIDTH) << "Error" << std::endl;
    
    for (unsigned i = 0; i < n; i++) {
        double exact = phi_exact(1. / n * ((double)i + .5));
        
        std::cout << std::setw(WIDTH) << x[i] << '\t';
        std::cout << std::setw(WIDTH) << exact << '\t';
        std::cout << std::setw(WIDTH) << exact - x[i] << std::endl;
    }
    
    std::cout << std::endl;
}


double error_norm(double* x, unsigned n) {
    double z_norm = 0;
    
    for (unsigned i = 0; i < n; i++) {
        double exact = phi_exact(1. / n * ((double)i + .5));
        double z_i = exact - x[i];
        
        if (z_norm < fabs(z_i))
            z_norm = fabs(z_i);
    }
    
    return z_norm;
}


int main(int argc, const char* argv[]) {
    const unsigned N = 100;
    double** a = new double*[N];
    for (unsigned i = 0; i < N; i++)
        a[i] = new double[N];
    double* f_i = new double[N];
    double* x = new double[N];
    
    for (unsigned i = 0; i < N; i++)
        x[i] = 0;
    
    create_matrixes(-1, ker, f, a, f_i, N);
//    print(a, f_i, N);
    
    double a_norm = matrix_norm(a, N);
    double a_inv_norm = inv_matrix_norm(a, N);
    
    unsigned it;
    double accuracy;
    int ret = solve(a, f_i, x, N, 1E-10, 1000, it, accuracy);
//    print(x, N);
    
    std::cout.flags(std::ios::scientific);
    
    std::cout << "Status: " << (ret? "accuracy not achieved" : "OK") << std::endl;
    std::cout << "Quantity of iterations: " << it << std::endl;
    std::cout << "Accuracy: " << accuracy << std::endl;
    std::cout << "||A|| = " << a_norm << std::endl;
    std::cout << "||A^-1|| = " << a_inv_norm << std::endl;
    std::cout << "ν = " << a_norm * a_inv_norm << std::endl;
    std::cout << "||z|| = " << error_norm(x, N) << std::endl << std::endl;
    
    for (unsigned i = 0; i < N; i++)
        delete[] a[i];
    delete[] a;
    delete[] f_i;
    delete[] x;
}
