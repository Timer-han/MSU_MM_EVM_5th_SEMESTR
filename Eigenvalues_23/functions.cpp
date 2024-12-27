#include "functions.h"

#include <cmath>
#include <cstdio>
#include <ctime>
#include <time.h>

double eps = 1e-10;

void fill_matrix(double *matrix, int n, int k)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            switch (k) {
                case 1:
                    matrix[i * n + j] = n - std::max(i, j);
                    break;
                case 2:
                    matrix[i * n + j] =
                        i == j ? 2 : std::abs(i - j) == 1 ? -1 : 0;
                    break;
                case 3:
                    matrix[i * n + j] =
                        ((i == j && j < n - 1)
                             ? 1
                             : (j == n - 1 ? i : (i == n - 1 ? j : 0)));
                    break;
                case 4:
                    matrix[i * n + j] = 1. / (i + j + 1);
                    break;
            }
        }
    }
}

int read_matrix(char *filename, double *matrix, int n)
{
    FILE *file = fopen(filename, "r");
    if (!file) return -1;

    for (int i = 0; i < n * n; i++) {
        if (fscanf(file, "%lf", &matrix[i]) != 1) return -2;
    }

    fclose(file);
    return 0;
}

void print_matrix(double *matrix, int n, int r)
{
    int rows = (r > n) ? n : r;
    int cols = (r > n) ? n : r;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf(" %10.3e", matrix[i * n + j]);
        }
        printf("\n");
    }
}

double scalar_product(double *vector1, double *vector2, int n)
{
    double result = 0;
    for (int i = 0; i < n; i++) result += vector1[i] * vector2[i];
    return result;
}

void vector_subtract(double *a, double *b, double *c, int n)
{  // a = b - c
    for (int i = 0; i < n; i++) {
        a[i] = b[i] - c[i];
    }
}

void vector_addition(double *a, double *b, double *c, int n)
{  // a = b + c
    for (int i = 0; i < n; i++) {
        a[i] = b[i] + c[i];
    }
}

void product_vect_scal(double *a, double b, int n)
{  // a = a * b
    for (int i = 0; i < n; i++) {
        a[i] *= b;
    }
}

void division_vect_scal(double *a, double b, int n)
{  // a = a / b
    for (int i = 0; i < n; i++) {
        a[i] /= b;
    }
}

void product_subtraction(double *a, double *b, double scal, int n)
{
    for (int i = 0; i < n; i++) a[i] -= b[i] * scal;
}


void print_vector(double *a, int n)
{
    for (int i = 0; i < n - 1; i++) printf("%8.3e ", a[i]);
    printf("%8.3e\n", a[n - 1]);
}

int triangolization(double *matrix, double *vector, int n)
{
    double norm, buf, scalar_prod;
    bool check;

    norm = 0;
    for (int j = 0; j < n; j++) norm += matrix[j * n + j];
    // printf("Trace = %8.3e\n", norm);

    for (int i = 0; i < n - 1; i++) {
        // printf("-------------- step %d --------------\n", i);
        // printf("Matrix:\n");
        // print_matrix(matrix, n, 4);

        norm = 0;
        for (int j = 0; j < n; j++) norm += matrix[j * n + j];
        // printf("Trace = %8.3e\n", norm);

        check = true;
        for (int j = i + 1; j < n; j++) {
            if (std::abs(matrix[i * n + j]) > eps) {
                check = false;
                break;
            }
        }
        if (check) continue;

        buf = 0;
        for (int j = i + 2; j < n; j++)
            buf += matrix[n * j + i] * matrix[n * j + i];
        norm = buf + matrix[n * (i + 1) + i] * matrix[n * (i + 1) + i];
        norm = sqrt(norm);
        // printf("Norm = %8.3e\n", norm);

        if (norm > eps) {
            for (int j = i + 1; j < n; j++) vector[j] = matrix[n * j + i];
            vector[i + 1] -= norm;
            // printf("vector is: ");
            // print_vector(vector + i + 1, n - i - 1);
            norm = vector[i + 1] * vector[i + 1] + buf;
            norm = sqrt(norm);

            if (norm < eps) {
                check = true;

                for (int row = i + 1; row < n && check; row++)
                    for (int col = 0; col < row && check; col++)
                        if (matrix[row * n + col] > eps) check = false;

                if (check) {
                    printf("[-] Division by zero!\n");
                    return -1;
                } else {
                    return 0;
                }
            }

            // x vector calculation
            for (int j = i + 1; j < n; j++) vector[j] /= norm;

            // printf("Yep!\n");
            // Multiplication x and matrix
            for (int col = i; col < n; col++) {
                scalar_prod = 0;
                for (int row = i + 1; row < n; row++)
                    scalar_prod += matrix[n * row + col] * vector[row];
                scalar_prod *= 2;
                for (int row = i + 1; row < n; row++)
                    matrix[n * row + col] -= scalar_prod * vector[row];
            }

            // Multiplication matrix and x
            for (int row = 0; row < n; row++) {
                scalar_prod = 0;
                for (int col = i + 1; col < n; col++)
                    scalar_prod += matrix[n * row + col] * vector[col];
                scalar_prod *= 2;
                for (int col = i + 1; col < n; col++)
                    matrix[n * row + col] -= scalar_prod * vector[col];
            }
        }
    }
    return 0;
}

int qr(double *matrix, double *cos, double *sin, int n, int r)
{
    double co, si, x, y, root;
    for (int i = 0; i < r; i++) {
        x = matrix[i * n + i];
        y = matrix[(i + 1) * n + i];

        if (std::abs(y) < eps) {
            cos[i] = 1;
            sin[i] = 0;
        } else {
            root = sqrt(x * x + y * y);
            if (root < eps) {
                return -1;
            }
            co = x / root;
            si = -y / root;

            cos[i] = co;
            sin[i] = si;

            matrix[i * n + i] = co * x - si * y;
            matrix[(i + 1) * n + i] = 0;

            for (int j = i + 1; j <= r; j++) {
                x = matrix[i * n + j];
                y = matrix[(i + 1) * n + j];
                matrix[i * n + j] = x * co - y * si;
                matrix[(i + 1) * n + j] = x * si + y * co;
            }
        }
    }
    return 0;
}

void rq(double *matrix, double *cos, double *sin, int n, int r)
{
    double si, co, x, y;
    for (int i = 0; i < r; i++) {
        si = sin[i];
        co = cos[i];
        for (int j = 0; j <= i; j++) {
            x = matrix[j * n + i];
            y = matrix[j * n + i + 1];
            matrix[j * n + i] = x * co - y * si;
            matrix[j * n + i + 1] = x * si + y * co;
        }
        y = matrix[(i + 1) * n + i + 1];
        matrix[(i + 1) * n + i] = -si * y;
        matrix[(i + 1) * n + i + 1] = co * y;
    }
}

double norm_calculation(double *matrix, int n) {
    double sum, norm = 0;
    for (int i = 0; i < n; i++) {
        sum = 0;
        for(int j = 0; j < n; j++) {
            sum += matrix[i * n + j] * matrix[i * n + j];
        }
        norm = std::max(norm, sum);
    }
    return sqrt(norm);
}

double lengh_calculation(double *matrix, int n) {
    double lengh = 0;
    for (int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            lengh += matrix[i * n + j] * matrix[j * n + i];
        }
    }
    return sqrt(lengh);
}

double res1_counter(double *eigenvals, int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++) sum += eigenvals[i];
    return std::abs(sum - matrix_trace) / matrix_norm;
}

double res2_counter(double *eigenvals, int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++) sum += eigenvals[i] * eigenvals[i];
    return std::abs(sqrt(sum) - matrix_lengh) / matrix_norm;
}

int run(double *matrix, double *vector, double *cos, double *sin,
        double *eigenvals, int n, int m, double &t1, double &t2, double &res1,
        double &res2, int &operation_count, double e)
{
    double shift;
    matrix_norm = norm_calculation(matrix, n);
    eps = matrix_norm * e;
    matrix_lengh = lengh_calculation(matrix, n);
    for (int j = 0; j < n; j++) {
        matrix_trace += matrix[j * n + j];
    }
    printf("Matrix: norm = %8.3e, trace\n", matrix_norm);

    t1 = clock();
    printf("\nTriangolized matrix:\n");
    if (triangolization(matrix, vector, n) != 0) {
        printf("[-] Can't triangolize!\n");
        return -2;
    }
    t1 = (clock() - t1) / CLOCKS_PER_SEC;

    print_matrix(matrix, n, m);
    printf("eps = %8.3e\n", eps);


    t2 = clock();
    for (int i = n - 1; i > 0; i--) {
        while (
            std::abs(matrix[i * n + i - 1]) > eps
        )
        {
            operation_count++;
            // printf("----------- step: %d ----------\n", operation_count);
            // printf("comparison with %8.3e\n", matrix[(i + 1) * n + i]);
            // printf("\nMatrix:\n");
            // print_matrix(matrix, n, m);

            shift = matrix[(i * n + i)] - 0.5 * matrix[i * n + i - 1];
            for (int j = 0; j <= i; j++) matrix[j * n + j] -= shift;

            if (qr(matrix, cos, sin, n, i) < 0) {
                printf("[-] Division by zero!");
                return -1;
            }
            rq(matrix, cos, sin, n, i);

            double trace = 0;
            for (int j = 0; j < n; j++) {
                trace += matrix[j * n + j];
            }
            for (int j = 0; j <= i; j++) matrix[j * n + j] += shift;
            // printf("Trace = %8.3e\n", trace);
        }
        eigenvals[i] = matrix[i * n + i];
    }
    eigenvals[0] = matrix[0];
    t2 = (clock() - t2) / CLOCKS_PER_SEC;

    res1 = res1_counter(eigenvals, n);
    res2 = res2_counter(eigenvals, n);

    printf("\nMatrix after:\n");
    print_matrix(matrix, n, m);

    return 0;
}