#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ABS(a) ((a) > 0 ? (a) : (-a))

double EPS = std::numeric_limits<double>::epsilon();

int run(double *matrix,
        double *inversed_matrix,
        double *block_A,
        double *block_B,
        double *block_C,
        size_t n,
        size_t m,
        size_t k,
        size_t l,
        size_t s,
        char *filename);
int find_diff(double *matrix, double *inversed_matrix, char* filename, int n, int m, int s, double &r1, double &r2);
int fill_matrix(double *matrix, size_t n, size_t s);
int read_matrix_from_file(double *matrix, size_t n, const char *filename);
void unit_matrix(double *matrix, size_t n);
void rows_permutation(double *A, double *block1, double *block2, size_t n,
                      size_t m, size_t k, size_t l, size_t i1, size_t i2,
                      size_t begin);
void print_matrix(double *matrix, size_t n, size_t r);
void print_matrix_l_x_n(double *matrix, size_t l, size_t n);
void get_block(double *a, double *block, size_t n, size_t m, size_t k, size_t l,
               size_t i, size_t j);
void put_block(double *a, double *block, size_t n, size_t m, size_t k, size_t l,
               size_t i, size_t j);
void matrix_sum(double *A, double *B, double *C, size_t n, size_t m);
void matrix_subtr(double *A, double *B, size_t n, size_t m);
double get_norm(double *matrix, size_t m);
int get_inverse_matrix(double *A, double *B, size_t m);
void mult(double *a, double *b, double *c, size_t n, size_t m);
void matrix_multiply(const double *A, const double *B, double *C, size_t p,
                     size_t q, size_t r);
void zero_matrix(double *matrix, size_t n, size_t m);


int find_diff(double *matrix, double *inversed_matrix, char* filename, int n, int m, int s, double &r1, double &r2)
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
    } else {
        r1 = 0;
        r2 = 0;
    }
    return 0;
}


int run(
        double *matrix,
        double *inversed_matrix,
        double *block_A,
        double *block_B,
        double *block_C,
        size_t n,
        size_t m,
        size_t k,
        size_t l,
        size_t s,
        char *filename
    )
{
    size_t i, j, column, row, min_norm_ind;
    double min_norm, norm;


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
            for (j = column + 1; j < k + 1; j++) {
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

            for (j = 0; j < k + 1; j++) {
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

    delete[] block_A;
    delete[] block_B;
    delete[] block_C;
    delete[] matrix;
    delete[] inversed_matrix;
    return 0;
}


int fill_matrix(double *matrix, size_t n, size_t s)
{
    size_t i, j;
    switch (s) {
    case 1:
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                matrix[i * n + j] = n - MAX(i, j) + 2;
        break;
    case 2:
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                matrix[i * n + j] = MAX(i, j) - 1;
        break;
    case 3:
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                matrix[i * n + j] = ABS(i - j);
        break;
    case 4:
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                matrix[i * n + j] = 1. / (i + j + 1);
        break;
    default:
        fprintf(stderr, "Unknown formula %ld\n", s);
        return -1;
    }
    return 0;
}

int read_matrix_from_file(double *matrix, size_t n, const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Can't open file %s\n", filename);
        return -1;
    }

    for (size_t i = 0; i < n * n; i++) {
        if (fscanf(file, "%lf", &matrix[i]) != 1) {
            fprintf(stderr, "Can't read file.\n");
            fclose(file);
            return -1;
        }
    }

    fclose(file);
    return 0;
}

void unit_matrix(double *matrix, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            matrix[i * n + j] = (i == j ? 1 : 0);
        }
    }
}

void zero_matrix(double *matrix, size_t n, size_t m)
{
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            matrix[i * n + j] = 0;
        }
    }
}

void matrix_multiply(const double *A, const double *B, double *C, size_t p,
                     size_t q, size_t r)
{
    size_t i, j, k;
    double temp;

    for (i = 0; i < p * r; ++i) {
        C[i] = 0.0;
    }

    for (i = 0; i < p; ++i) {
        for (k = 0; k < q; ++k) {
            temp = A[i * q + k];

            for (j = 0; j + 4 <= r; j += 4) {
                C[i * r + j] += temp * B[k * r + j];
                C[i * r + j + 1] += temp * B[k * r + j + 1];
                C[i * r + j + 2] += temp * B[k * r + j + 2];
                C[i * r + j + 3] += temp * B[k * r + j + 3];
            }

            for (; j < r; ++j) {
                C[i * r + j] += temp * B[k * r + j];
            }
        }
    }
}

void rows_permutation(double *A, double *block1, double *block2, size_t n,
                      size_t m, size_t k, size_t l, size_t i1, size_t i2,
                      size_t begin)
{
    size_t i;
    for (i = begin; i < k + 1; i++) {
        get_block(A, block1, n, m, k, l, i1, i);
        get_block(A, block2, n, m, k, l, i2, i);
        put_block(A, block2, n, m, k, l, i1, i);
        put_block(A, block1, n, m, k, l, i2, i);
    }
}

void matrix_sum(double *A, double *B, double *C, size_t n,
                size_t m) // C = A + B
{
    size_t i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            C[j + i * m] = A[j + i * m] + B[j + i * m];
        }
    }
}

void matrix_subtr(double *A, double *B, size_t n, size_t m) // A -= B
{
    size_t i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            A[j + i * m] -= B[j + i * m];
        }
    }
}

double get_norm(double *matrix, size_t m)
{
    size_t i, j;
    double max = 0, sum;
    for (j = 0; j < m; j++) {
        sum = 0;
        for (i = 0; i < m; i++) {
            sum += ABS(matrix[j + i * m]);
        }
        max = MAX(max, sum);
    }
    return max;
}

int get_inverse_matrix(double *A, double *B, size_t m)
{
    size_t i, j, k, max_ind;
    double max, buf;
    unit_matrix(B, m);
    for (j = 0; j < m; j++) {
        max = 0;
        max_ind = 0;
        for (i = j; i < m; i++) {
            if (ABS(A[j + i * m]) > max) {
                max = ABS(A[j + i * m]);
                max_ind = i;
            }
        }
        printf("max: %lf, max_ind: %ld\n", max, max_ind);

        if (max <= EPS) {
            return -1;
        }
        for (i = j; i < m; i++) {
            buf = A[i + j * m];
            A[i + j * m] = A[i + max_ind * m];
            A[i + max_ind * m] = buf;
        }

        for (i = 0; i < m; i++) {
            buf = B[i + j * m];
            B[i + j * m] = B[i + max_ind * m];
            B[i + max_ind * m] = buf;
        }

        A[j + j * m] = 1;
        for (i = j + 1; i < m; i++) {
            A[i + j * m] /= max;
        }

        for (i = 0; i < m; i++) {
            B[i + j * m] /= max;
        }

        for (i = 0; i < m; i++) {
            if (i == j)
                continue;
            max = A[j + i * m];
            A[j + i * m] = 0;
            for (k = j + 1; k < m; k++) {
                A[k + i * m] -= max * A[k + j * m];
            }
            for (k = 0; k < m; k++) {
                B[k + i * m] -= max * B[k + j * m];
            }
        }
    }
    return 0;
}

void mult(double *a, double *b, double *c, size_t n, size_t m)
{
    size_t i, j, s, r, t, q;
    size_t k = n / m;
    size_t l = n - k * m; // n = k * m + l
    size_t bl = (l != 0 ? k + 1 : k); // Общее количество блоков
    size_t v, h, ah;

    // Проходим по всем блокам в матрице C
    for (i = 0; i < bl; i++) {
        for (j = 0; j < bl; j++) {
            // Определяем размер текущего блока C[v x h]
            v = (i < k ? m : l); // вертикальный размер блока
            h = (j < k ? m : l); // горизонтальный размер блока

            // Указатель на начало текущего блока C
            double *pc = c + (i * m) * n + j * m;

            // Инициализируем блок C нулями
            for (r = 0; r < v; r++) {
                for (t = 0; t < h; t++) {
                    pc[r * n + t] = 0.0;
                }
            }

            // Перемножаем соответствующие блоки A и B и добавляем к C
            for (s = 0; s < bl; s++) {
                // Определяем размер внутреннего блока
                ah = (s < k ? m : l);

                // Указатели на текущие блоки A и B
                double *pa = a + (i * m) * n + s * m; // блок A[i][s]
                double *pb = b + (s * m) * n + j * m; // блок B[s][j]

                // Основные циклы с разверткой для блоков 3x3
                size_t r_end = (v / 3) * 3;
                size_t t_end = (h / 3) * 3;

                // Обработка блоков размером 3x3
                for (r = 0; r < r_end; r += 3) {
                    for (t = 0; t < t_end; t += 3) {
                        double s00 = 0.0, s01 = 0.0, s02 = 0.0;
                        double s10 = 0.0, s11 = 0.0, s12 = 0.0;
                        double s20 = 0.0, s21 = 0.0, s22 = 0.0;

                        for (q = 0; q < ah; q++) {
                            double a0q = pa[(r + 0) * n + q];
                            double a1q = pa[(r + 1) * n + q];
                            double a2q = pa[(r + 2) * n + q];
                            double bq0 = pb[q * n + (t + 0)];
                            double bq1 = pb[q * n + (t + 1)];
                            double bq2 = pb[q * n + (t + 2)];

                            s00 += a0q * bq0;
                            s01 += a0q * bq1;
                            s02 += a0q * bq2;

                            s10 += a1q * bq0;
                            s11 += a1q * bq1;
                            s12 += a1q * bq2;

                            s20 += a2q * bq0;
                            s21 += a2q * bq1;
                            s22 += a2q * bq2;
                        }

                        pc[(r + 0) * n + (t + 0)] += s00;
                        pc[(r + 0) * n + (t + 1)] += s01;
                        pc[(r + 0) * n + (t + 2)] += s02;

                        pc[(r + 1) * n + (t + 0)] += s10;
                        pc[(r + 1) * n + (t + 1)] += s11;
                        pc[(r + 1) * n + (t + 2)] += s12;

                        pc[(r + 2) * n + (t + 0)] += s20;
                        pc[(r + 2) * n + (t + 1)] += s21;
                        pc[(r + 2) * n + (t + 2)] += s22;
                    }

                    // Обработка оставшихся столбцов в блоке
                    for (t = t_end; t < h; t++) {
                        double s0 = 0.0, s1 = 0.0, s2 = 0.0;

                        for (q = 0; q < ah; q++) {
                            double a0q = pa[(r + 0) * n + q];
                            double a1q = pa[(r + 1) * n + q];
                            double a2q = pa[(r + 2) * n + q];
                            double bqt = pb[q * n + t];

                            s0 += a0q * bqt;
                            s1 += a1q * bqt;
                            s2 += a2q * bqt;
                        }

                        pc[(r + 0) * n + t] += s0;
                        pc[(r + 1) * n + t] += s1;
                        pc[(r + 2) * n + t] += s2;
                    }
                }

                // Обработка оставшихся строк в блоке
                for (r = r_end; r < v; r++) {
                    for (t = 0; t < t_end; t += 3) {
                        double s0 = 0.0, s1 = 0.0, s2 = 0.0;

                        for (q = 0; q < ah; q++) {
                            double a0q = pa[r * n + q];
                            double bq0 = pb[q * n + (t + 0)];
                            double bq1 = pb[q * n + (t + 1)];
                            double bq2 = pb[q * n + (t + 2)];

                            s0 += a0q * bq0;
                            s1 += a0q * bq1;
                            s2 += a0q * bq2;
                        }

                        pc[r * n + (t + 0)] += s0;
                        pc[r * n + (t + 1)] += s1;
                        pc[r * n + (t + 2)] += s2;
                    }

                    // Обработка оставшихся элементов
                    for (t = t_end; t < h; t++) {
                        double sum = 0.0;

                        for (q = 0; q < ah; q++) {
                            sum += pa[r * n + q] * pb[q * n + t];
                        }

                        pc[r * n + t] += sum;
                    }
                }
            }
        }
    }
}

void print_matrix(double *matrix, size_t n, size_t r)
{
    size_t rows = (r > n) ? n : r;
    size_t cols = (r > n) ? n : r;
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            printf(" %10.3e", matrix[i * n + j]);
        }
        printf("\n");
    }
}

void print_matrix_l_x_n(double *matrix, size_t l, size_t n)
{
    for (size_t i = 0; i < l; i++) {
        for (size_t j = 0; j < n; j++) {
            printf(" %10.3e", matrix[i * n + j]);
        }
        printf("\n");
    }
}

void get_block(double *a, double *block, size_t n, size_t m, size_t k, size_t l,
               size_t i, size_t j)
{
    size_t block_index = 0;
    size_t count_skip_elements = 0;
    double *block_beginning = nullptr;
    if (i < k) {
        block_beginning = a + m * (m * k * i + i * l + m * j);
        if (j < k) {
            for (block_index = 0; block_index < m * m; block_index++) {
                block[block_index] =
                    block_beginning[block_index + count_skip_elements];
                if ((block_index + 1) % m == 0) {
                    count_skip_elements += n - m;
                }
            }
        } else {
            for (block_index = 0; block_index < m * l; block_index++) {
                block[block_index] =
                    block_beginning[block_index + count_skip_elements];
                if ((block_index + 1) % l == 0) {
                    count_skip_elements += n - l;
                }
            }
        }
    } else {
        block_beginning = a + m * (m * k * k + l * k + l * j);
        if (j < k) {
            for (block_index = 0; block_index < m * l; block_index++) {
                block[block_index] =
                    block_beginning[block_index + count_skip_elements];
                if ((block_index + 1) % m == 0) {
                    count_skip_elements += n - m;
                }
            }
        } else {
            for (block_index = 0; block_index < l * l; block_index++) {
                block[block_index] =
                    block_beginning[block_index + count_skip_elements];
                if ((block_index + 1) % l == 0) {
                    count_skip_elements += n - l;
                }
            }
        }
    }
}

void put_block(double *a, double *block, size_t n, size_t m, size_t k, size_t l,
               size_t i, size_t j)
{
    size_t block_index = 0;
    size_t count_skip_elements = 0;
    double *block_beginning = nullptr;
    if (i < k) {
        block_beginning = a + m * (m * k * i + i * l + m * j);
        if (j < k) {
            for (block_index = 0; block_index < m * m; block_index++) {
                block_beginning[block_index + count_skip_elements] =
                    block[block_index];
                if ((block_index + 1) % m == 0) {
                    count_skip_elements += n - m;
                }
            }
        } else {
            for (block_index = 0; block_index < m * l; block_index++) {
                block_beginning[block_index + count_skip_elements] =
                    block[block_index];
                if ((block_index + 1) % l == 0) {
                    count_skip_elements += n - l;
                }
            }
        }
    } else {
        block_beginning = a + m * (m * k * k + l * k + l * j);
        if (j < k) {
            for (block_index = 0; block_index < m * l; block_index++) {
                block_beginning[block_index + count_skip_elements] =
                    block[block_index];
                if ((block_index + 1) % l == 0) {
                    count_skip_elements += n - l;
                }
            }
        } else {
            for (block_index = 0; block_index < l * l; block_index++) {
                block_beginning[block_index + count_skip_elements] =
                    block[block_index];
                if ((block_index + 1) % l == 0) {
                    count_skip_elements += n - l;
                }
            }
        }
    }
}