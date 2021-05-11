#include "invert.hpp"


int sign(double x) {
    return x > 0? 1 : -1;
}


void swap_str(double** a, unsigned n, unsigned str1, unsigned str2, unsigned from) {
    for (unsigned j = from; j < n; j++)
        std::swap(a[str1][j], a[str2][j]);
}


void swap_col(double** a, unsigned n, unsigned col1, unsigned col2) {
    for (unsigned i = 0; i < n; i++)
        std::swap(a[i][col1], a[i][col2]);
}


void up_triangularize(double** a, unsigned* col_transp, unsigned* str_transp, unsigned n) {
    for (unsigned j = 0; j < n - 1; j++) {
        if (j < n - 1) {
            double max_len = 0;
            unsigned max_len_col = j;
            for (unsigned i = 0; i < n; i++)
                max_len += a[i][j] * a[i][j];
            
            for (unsigned k = j + 1; k < n; k++) {
                double len = 0;
                for (unsigned i = 0; i < n; i++)
                    len += a[i][k] * a[i][k];
                
                if (len > max_len) {
                    max_len = len;
                    max_len_col = k;
                }
            }
            
            col_transp[j] = max_len_col;
            swap_col(a, n, j, max_len_col);
            
            double max_elem = fabs(a[j + 1][j]);
            unsigned max_elem_str = j + 1;
            for (unsigned i = j + 2; i < n; i++)
                if (fabs(a[i][j]) > max_elem) {
                    max_elem = fabs(a[i][j]);
                    max_elem_str = i;
                }
            
            str_transp[j + 1] = max_elem_str;
            swap_str(a, n, j + 1, max_elem_str, j);
        }
        
        for (unsigned i = j + 1; i < n; i++) {
            double a_i = a[i][j];
            
            if (a_i != 0) {
                double a_j = a[j][j];
                double z = fmax(fabs(a_j), fabs(a_i));
                double a_l = fmin(fabs(a_j), fabs(a_i)) / z;
                double den = sqrt(1 + a_l * a_l);
                double s = -a_i / z / den;
                double c = a_j / z / den;
                
                a[j][j] = c * a_j - s * a_i;
                a[i][j] = c + sign(s);
                
                for (unsigned k = j + 1; k < n; k++) {
                    double kj = a[j][k];
                    
                    a[j][k] = c * kj - s * a[i][k];
                    a[i][k] = s * kj + c * a[i][k];
                }
            }
        }
    }
}


void invert_tr(double** a, unsigned n) {
    for (unsigned i = n - 1; i >= 0 && i < UINT_MAX; i--) {
        for (unsigned j = n - 1; j >= i && j < UINT_MAX; j--) {
            if (i != j) {
                double ij = 0;
                for (unsigned k = i + 1; k <= j; k++)
                    ij += a[i][k] * a[k][j];
                a[i][j] = -1 / a[i][i] * ij;
            }
            
            else
                a[j][i] = 1 / a[i][i];
        }
    }
}


void de_up_triangularize(double** a, unsigned* col_transp, unsigned* str_transp, unsigned n) {
    for (unsigned j = n - 2; j >= 0 && j < UINT_MAX; j--) {
        double q[n];
        for (unsigned i = j + 1; i < n; i++) {
            q[i] = a[i][j];
            a[i][j] = 0;
        }
        
        for (unsigned i = n - 1; i > j; i--)
            if (q[i] != 0) {
                double s, c;
                
                if (q[i] < 0) {
                    c = q[i] + 1;
                    s = -sqrt(1 - c * c);
                }
                
                else {
                    c = q[i] - 1;
                    s = sqrt(1 - c * c);
                }
                
                for (unsigned k = 0; k < n; k++) {
                    double jk = a[k][j];
                    
                    a[k][j] = c * jk + s * a[k][i];
                    a[k][i] = -s * jk + c * a[k][i];
                }
            }
        
        swap_col(a, n, j + 1, str_transp[j + 1]);
        swap_str(a, n, j, col_transp[j], j);
    }
}


void invert(double** a, unsigned n) {
    unsigned* col_transp = new unsigned[n];
    unsigned* str_transp = new unsigned[n];
    
    up_triangularize(a, col_transp, str_transp, n);
    invert_tr(a, n);
    de_up_triangularize(a, col_transp, str_transp, n);
}
