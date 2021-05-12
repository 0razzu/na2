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


double vector_norm(double* v, unsigned n) {
    double norm = 0;
    
    for (unsigned i = 0; i < n; i++)
        if (fabs(v[i]) > norm)
            norm = fabs(v[i]);
    
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


void print(double** a, double* f_n, unsigned n) {
    std::cout.setf(std::ios::right | std::ios::fixed | std::ios::showpoint);
    std::cout.precision(PRECISION);
    
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j < n; j++)
            std::cout << std::setw(WIDTH) << a[i][j] << '\t';
        std::cout << '\t' << std::setw(WIDTH) << f_n[i] << std::endl;
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
    double norm = 0;
    
    for (unsigned i = 0; i < n; i++) {
        double exact = phi_exact(1. / n * ((double)i + .5));
        double z_i = exact - x[i];
        
        if (norm < fabs(z_i))
            norm = fabs(z_i);
    }
    
    return norm;
}


double discrepancy_norm(double** a, double* x, double* f_n, unsigned n) {
    double norm = 0;
    
    for (unsigned i = 0; i < n; i++) {
        double r_i = 0;
        
        for (unsigned j = 0; j < n; j++)
            r_i += a[i][j] * x[j];
        
        r_i -= f_n[i];
        
        if (fabs(r_i) > norm)
            norm = fabs(r_i);
    }
    
    return norm;
}


int main(int argc, const char* argv[]) {
    const unsigned N = 100;
    double** a = new double*[N];
    for (unsigned i = 0; i < N; i++)
        a[i] = new double[N];
    double* f_n = new double[N];
    double* x = new double[N];
    for (unsigned i = 0; i < N; i++)
        x[i] = 0;
    
    create_matrixes(-1, ker, f, a, f_n, N);
//    print(a, f_n, N);
    
    double** a_copy = new double*[N];
    for (unsigned i = 0; i < N; i++) {
        a_copy[i] = new double[N];
        
        for (unsigned j = 0; j < N; j++)
            a_copy[i][j] = a[i][j];
    }
    double* f_copy = new double[N];
    for (unsigned i = 0; i < N; i++)
        f_copy[i] = f_n[i];
    
    unsigned it;
    double accuracy;
    int ret = solve(a, f_n, x, N, 1E-10, 1000, it, accuracy);
//    print(x, N);
    
    double a_norm = matrix_norm(a_copy, N);
    double a_inv_norm = inv_matrix_norm(a_copy, N);
    double x_norm = vector_norm(x, N);
    double z_norm = error_norm(x, N);
    double r_norm = discrepancy_norm(a_copy, x, f_copy, N);
    
    std::cout.flags(std::ios::scientific);
    
    std::cout << "Status: " << (ret? "accuracy not achieved" : "OK") << std::endl;
    std::cout << "Quantity of iterations: " << it << std::endl;
    std::cout << "Accuracy: " << accuracy << std::endl;
    std::cout << "||A|| = " << a_norm << std::endl;
    std::cout << "||A^-1|| = " << a_inv_norm << std::endl;
    std::cout << "ν(A) = " << a_norm * a_inv_norm << std::endl;
    std::cout << "||x|| = " << x_norm << std::endl;
    std::cout << "||z|| = " << z_norm << std::endl;
    std::cout << "ζ = " << z_norm / x_norm << std::endl;
    std::cout << "||r|| = " << r_norm << std::endl;
    std::cout << "ρ = " << r_norm / vector_norm(f_copy, N) << std::endl << std::endl;
    
    for (unsigned i = 0; i < N; i++) {
        delete[] a[i];
        delete[] a_copy[i];
    }
    delete[] a;
    delete[] a_copy;
    delete[] f_n;
    delete[] f_copy;
    delete[] x;
}
