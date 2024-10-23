#include "functions.h"
#include <time.h>

int main(int argc, char *argv[])
{
    size_t i, j, k, column, row, min_norm_ind;
    double min_norm, norm, r1 = 0, r2 = 0;

    if (argc < 5) {
        fprintf(stderr, "Not enough arguments.\n");
        return 1;
    }

    size_t n = atoi(argv[1]);
    size_t m = atoi(argv[2]);
    size_t r = atoi(argv[3]);
    size_t s = atoi(argv[4]);
    size_t l = n % m;
    k = n / m;

    while (1. + EPS > 1.) {
        EPS /= 2.;
    }
    EPS *= 4;

    char *filename = nullptr;
    if (s == 0) {
        if (argc < 6) {
            fprintf(stderr, "File name do not defined.\n");
            return 1;
        }
        filename = argv[5];
    }

    double *matrix = new double[n*n];
    double *inversed_matrix =  new double[n*n];
    double *block_A = new double[m*m];
    double *block_B = new double[m*m];
    double *block_C = new double[m*m];

    clock_t t1 = clock();
    if (!block_C || !block_B || !block_C || !matrix || !inversed_matrix ||
        !run(matrix, inversed_matrix, block_A, block_B, block_C, n, m, k, l, s, filename)) {
        if (block_A) delete[] block_A;
        if (block_B) delete[] block_B;
        if (block_C) delete[] block_C;
        if (matrix) delete[] matrix;
        if (inversed_matrix) delete[] inversed_matrix;
        if (!block_C || !block_B || !block_C || !matrix || !inversed_matrix) {
            fprintf(stderr, "Can't allocate memory for matrix\n");
        }
        return 1;
    }

    print_matrix(inversed_matrix, n, r);
    t1 -= clock_t();

    clock_t t2 = clock();
    find_diff(matrix, inversed_matrix, n, m, s, r1, r2);
    t2 -= clock_t();


    printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %ld N = "
           "%ld M = %ld\n",
           argv[0], 18, r1, r2, double(t1 / CLOCKS_PER_SEC),
           double(t2 / CLOCKS_PER_SEC), s, n, m);

    delete[] block_A;
    delete[] block_B;
    delete[] block_C;
    delete[] matrix;
    delete[] inversed_matrix;
    return 0;
}


int find_diff(double *matrix, double *inversed_matrix, int n, int m, int s, double &r1, double &r2)
{
    if (n < 11000) {
        double *buf = new double[n * n];
        if (!buf) {
            fprintf(stderr, "Can't allocate memory for matrix\n");
            return 1;
        }
        if (s == 0) {
            if (read_matrix_from_file(buf, n, filename) != 0) {
                return 2;
            }
        } else {
            if (fill_matrix(buf, n, s) != 0) {
                return 2;
            }
        }
        mult(buf, inversed_matrix, matrix, n, m);
        unit_matrix(buf, n);
        matrix_subtr(buf, inversed_matrix, n, n);
        r1 = get_norm(buf, n);
        if (s == 0) {
            if (read_matrix_from_file(buf, n, filename) != 0) {
                return 2;
            }
        } else {
            if (fill_matrix(buf, n, s) != 0) {
                return 2;
            }
        }
        mult(inversed_matrix, buf, matrix, n, m);
        unit_matrix(buf, n);
        matrix_subtr(buf, inversed_matrix, n, n);
        r2 = get_norm(buf, n);

        delete[] buf;
    }
}


int run(double *matrix,
        double *inversed_matrix,
        double *block_A,
        double *block_B,
        double *block_C,
        size_t n,
        size_t m,
        size_t k,
        size_t n,
        size_t n,
        size_t s,
        char *filename
        )
{
    size_t i, j, k, column, row, min_norm_ind;
    double min_norm, norm, r1 = 0, r2 = 0;


    if (s == 0) {
        if (read_matrix_from_file(matrix, n, filename) != 0) {
            return 2;
        }
    } else {
        if (fill_matrix(matrix, n, s) != 0) {
            return 2;
        }
    }


    unit_matrix(inversed_matrix, n);

    for (column = 0; column < k + 1; column++) {
        // print_matrix(matrix, n, n);
        // printf("\n");
        min_norm = -1;
        min_norm_ind = column;
        for (row = column; row < k; row++) {
            get_block(matrix, block_A, n, m, k, l, column, row);
            print_matrix(block_A, m, m);
            printf("\n");
            if (column != k) {
                if (get_inverse_matrix(block_A, block_B, m) != 0) {
                    continue;
                }
                print_matrix(block_B, m, m);
            } else {
                if (get_inverse_matrix(block_A, block_B, l) != 0) {
                    continue;
                }
                print_matrix(block_B, l, l);
            }


            norm = get_norm(block_B, m);
            if (min_norm < 0 || norm < min_norm) {
                min_norm = norm;
                min_norm_ind = row;
            }
        }
        printf("min_norm: %lf\n", min_norm);
        if (min_norm < 0) {
            fprintf(stderr, "Matrix is invertable!\n");
            return -2;
        }

        rows_permutation(matrix, block_A, block_B, n, m, k, l, min_norm_ind, column,
                         column);
        rows_permutation(inversed_matrix, block_A, block_B, n, m, k, l,
                         min_norm_ind, column, 0);

        get_block(matrix, block_A, n, m, k, l, column, column);
        if (column != k) {
            unit_matrix(block_A, m);
            if (get_inverse_matrix(block_A, block_B, m) != 0) {
                fprintf(stderr, "Matrix is invertable!\n");
                return -1;
            }
        } else {
            unit_matrix(block_A, l);
            if (get_inverse_matrix(block_A, block_B, l) != 0) {
                fprintf(stderr, "Matrix is invertable!\n");
                return -1;
            }
        }

        put_block(matrix, block_A, n, m, k, l, column, column);

        for (i = column + 1; i < k + 1; i++) {
            get_block(matrix, block_A, n, m, k, l, i, column);
            if (column != k && i != k) {
                matrix_multiply(block_A, block_B, block_C, m, m, m);
            } else if (column != k && i == k) {
                matrix_multiply(block_B, block_A, block_C, m, m, l);
            }
            put_block(matrix, block_C, n, m, k, l, i, column);
        }

        for (i = 0; i < k + 1; i++) {
            get_block(inversed_matrix, block_A, n, m, k, l, i, column);
            if (column != k && i != k) {
                matrix_multiply(block_A, block_B, block_C, m, m, m);
            } else if (column != k && i == k) {
                matrix_multiply(block_B, block_A, block_C, m, m, l);
            } else if (column == k && i != k) {
                matrix_multiply(block_B, block_A, block_C, l, l, m);
            } else {
                matrix_multiply(block_B, block_A, block_C, l, l, l);
            }
            put_block(inversed_matrix, block_C, n, m, k, l, i, column);
        }

        for (i = 0; i < k + 1; i++) {
            if (i == column)
                continue;
            get_block(matrix, block_A, n, m, k, l, i, column);
            if (column == k) {
                zero_matrix(block_C, m, l);
            } else {
                zero_matrix(block_C, m, m);
            }
            put_block(matrix, block_C, n, m, k, l, i, column);
            for (j = column + 1; j < k + 1; k++) {
                get_block(matrix, block_C, n, m, k, l, i, j);
                if (i != k && j != k) {
                    matrix_multiply(block_C, block_A, block_B, m, m, m);
                } else if (i != k && j == k) {
                    matrix_multiply(block_C, block_A, block_B, m, m, l);
                } else if (i == k && j != k) {
                    matrix_multiply(block_A, block_C, block_B, l, m, m);
                } else {
                    matrix_multiply(block_A, block_C, block_B, l, m, l);
                }
                get_block(matrix, block_C, n, m, k, l, column, j);
                matrix_subtr(block_C, block_B, (i != k ? m : l),
                             (j != k ? m : l));
                put_block(matrix, block_C, n, m, k, l, column, j);
            }

            for (j = 0; j < k + 1; k++) {
                get_block(inversed_matrix, block_C, n, m, k, l, i, j);
                if (i != k && j != k) {
                    matrix_multiply(block_C, block_A, block_B, m, m, m);
                } else if (i != k && j == k) {
                    matrix_multiply(block_C, block_A, block_B, m, m, l);
                } else if (i == k && j != k) {
                    matrix_multiply(block_A, block_C, block_B, l, m, m);
                } else {
                    matrix_multiply(block_A, block_C, block_B, l, m, l);
                }
                get_block(inversed_matrix, block_C, n, m, k, l, column, j);
                matrix_subtr(block_C, block_B, (i != k ? m : l),
                             (j != k ? m : l));
                put_block(inversed_matrix, block_C, n, m, k, l, column, j);
            }
        }
    }

    print_matrix(inversed_matrix, n, r);
    t1 -= clock_t();

    clock_t t2 = clock();
    if (n < 11000) {
        double *buf = (double *)malloc(n * n * sizeof(double));
        if (!buf) {
            fprintf(stderr, "Can't allocate memory for matrix\n");
            return 1;
        }
        if (s == 0) {
            if (read_matrix_from_file(buf, n, filename) != 0) {
                return 2;
            }
        } else {
            if (fill_matrix(buf, n, s) != 0) {
                return 2;
            }
        }
        mult(buf, inversed_matrix, matrix, n, m);
        unit_matrix(buf, n);
        matrix_subtr(buf, inversed_matrix, n, n);
        r1 = get_norm(buf, n);
        if (s == 0) {
            if (read_matrix_from_file(buf, n, filename) != 0) {
                return 2;
            }
        } else {
            if (fill_matrix(buf, n, s) != 0) {
                return 2;
            }
        }
        mult(inversed_matrix, buf, matrix, n, m);
        unit_matrix(buf, n);
        matrix_subtr(buf, inversed_matrix, n, n);
        r2 = get_norm(buf, n);

        free(buf);
    }
    t2 -= clock_t();


    printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %ld N = "
           "%ld M = %ld\n",
           argv[0], 18, r1, r2, double(t1 / CLOCKS_PER_SEC),
           double(t2 / CLOCKS_PER_SEC), s, n, m);

    delete[] block_A;
    delete[] block_B;
    delete[] block_C;
    delete[] matrix;
    delete[] inversed_matrix;
    return 0;
}