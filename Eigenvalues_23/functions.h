#include <cmath>
#include <cstring>
#include <iostream>

void fill_matrix(double *matrix, int n, int k);
int read_matrix(char *filename, double *matrix, int n);
void print_matrix(double *matrix, int n, int r);
double scalar_product(double *vector1, double *vector2, int n);
void vector_subtract(double *a, double *b, double *c, int n);
void vector_addition(double *a, double *b, double *c, int n);
void product_vect_scal(double *a, double b, int n);
void division_vect_scal(double *a, double b, int n);
void product_subtraction(double *a, double *b, double scal, int n);
double norm_calculation(double *vector, int n);
int triangolization(double *matrix, double *vector, int n);
int run(double *matrix, double *vector, double *cos, double *sin,
        double *eigenvals, int n, int m, double &t1, double &t2, double &res1,
        double &res2, int &operation_count, double e);

static double matrix_norm = 1e-10;
static double matrix_trace = 0;
static double matrix_lengh = 0;
