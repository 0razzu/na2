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

    std::cout.setf(std::ios::right | std::ios::fixed | std::ios::showpoint);
    std::cout.precision(PRECISION);
    
    for (unsigned i = 0; i < N; i++) {
        for (unsigned j = 0; j < N; j++)
            std::cout << std::setw(WIDTH) << a[i][j] << '\t';
        std::cout << '\t' << f_i[i] << '\n';
    }
    
    unsigned it;
    double accuracy;
    int ret = solve(a, f_i, x, N, 1E-10, 1000, it, accuracy);
    
    std::cout << std::endl;
    std::cout.flags(std::ios::left | std::ios::fixed | std::ios::showpoint);
    std::cout << std::setw(WIDTH) << "Computed" << '\t';
    std::cout << std::setw(WIDTH) << "Exact" << '\t';
    std::cout << std::setw(WIDTH) << "Error" << std::endl;
    
    double z_norm = 0;
    for (unsigned i = 0; i < N; i++) {
        double exact = phi_exact(1. / N * ((double)i + .5));
        double z_i = exact - x[i];
        
        if (z_norm < fabs(z_i))
            z_norm = fabs(z_i);
        
        std::cout << std::setw(WIDTH) << x[i] << '\t';
        std::cout << std::setw(WIDTH) << exact << '\t';
        std::cout << std::setw(WIDTH) << z_i << std::endl;
    }
    
    std::cout.flags(std::ios::scientific);
    
    std::cout << std::endl;
    std::cout << "Status: " << (ret? "accuracy not achieved" : "OK") << std::endl;
    std::cout << "Quantity of iterations: " << it << std::endl;
    std::cout << "Accuracy: " << accuracy << std::endl;
    std::cout << "Error norm: " << z_norm << std::endl << std::endl;
    
    for (unsigned i = 0; i < N; i++)
        delete[] a[i];
    delete[] a;
    delete[] f_i;
    delete[] x;
}
