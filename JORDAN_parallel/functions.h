#include <iostream>
#include <cfloat>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <time.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define ABS(a) ((a) > 0 ? (a) : (-a))

static double EPS = 1e-16;


enum class io_status
{
    bad_allocation,
	error_open,
	error_read,
    irreversible,
	success,
	undef,
	unknown_formula
};


enum class reduce
{
    abs_max,
    abs_min,
    abs_max_first,
    abs_min_first,
    max,
    min,
    mult,
	sum,
    synchronize
};


class Args
{
public:
    // Данные для отладки
	size_t p = 0; // количество потоков
	size_t pi = 0; // номер потока
	pthread_t tid = 1;
	const char * name = nullptr;
	io_status error_type = io_status::undef;
	double error_flag = 0;

    // Данные для вычислений
	size_t n = 0;  // размер матрицы
	size_t m = 0;  // размер блока
	size_t l = 0;  // размер последнего блока
	size_t k = 0;  // количество блоков размера m в строке
	size_t s = 0;  // номер формулы
	size_t r = 0;  // размер выводимой
    double r1 = 0; // невязка 1
    double r2 = 0; // невязка 2
    double t1 = 0; // time 1
    double t2 = 0; // time 2
    FILE *file;
    double *block = nullptr;
    double *matrix = nullptr;
    double *inversed_matrix = nullptr;
    double *thread_part = nullptr;
    double *norm = nullptr;
};

void synchronize(int p, double *a, int n, reduce reduce_type);
void *thread_func(void *args);
// int run(double *matrix,
//         double *inversed_matrix,
//         double *block_A,
//         double *block_B,
//         double *block_C,
//         size_t n,
//         size_t m,
//         size_t k,
//         size_t l,
//         size_t s,
//         size_t r,
//         char *filename);
int find_diff(double *matrix, double *inversed_matrix, double *block, double *norm, FILE*file, int n, int m, int s, double &r1, double &r2, size_t p, size_t pi);
int fill_matrix(double *matrix, size_t n, size_t s);
int read_matrix_from_file(double *matrix, size_t n, FILE *file);
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
// void mult(double *a, double *b, double *c, size_t n, size_t m);
double mult_sub_norm_p(double *a, double *b, double *pc, double *norm, size_t n, size_t m, size_t p, size_t pi);
void matrix_multiply(const double *A, const double *B, double *C, size_t p,
                     size_t q, size_t r);
void zero_matrix(double *matrix, size_t n, size_t m);

void zero_matrix_p(double *matrix, size_t n, size_t m, size_t p, size_t pi);
void unit_matrix_p(double *matrix, size_t n, size_t m, size_t p, size_t pi);
int fill_matrix_p(double *matrix, size_t n, size_t s, size_t m, size_t p, size_t pi);
double get_norm_p(double *matrix, size_t n, size_t m, size_t p, size_t pi);
int process_args(Args *a);
void rows_permutation_p(double *A, double *block1, double *block2, size_t n,
                      size_t m, size_t k, size_t l, size_t i1, size_t i2,
                      size_t begin, size_t p, size_t pi);
void *thread_func(void *args);
double get_time();


void synchronize(int p, double *a=nullptr, int n = 0, reduce reduce_type = reduce::abs_max) {
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
	static int t_in = 0;//количество вошедших потоков
	static int t_out = 0;//количество вышедших потоков
	static double *r = nullptr;
	// static int ind = -1;
    // static double buf = 0;
	if (p <= 1) return;
	pthread_mutex_lock(&m);
	if (r == nullptr) r = a;
	else {
        if (reduce_type == reduce::sum) for(int i = 0; i < n; i++) r[i] += a[i];
        if (reduce_type == reduce::max)
            for(int i = 0; i < n; i++)
                if (a[i] > r[i])
                    r[i] = a[i];
        if (reduce_type == reduce::min)
            for(int i = 0; i < n; i++)
                if (a[i] < r[i])
                    r[i] = a[i];
        if (reduce_type == reduce::abs_max)
            for(int i = 0; i < n; i++)
                if (std::abs(a[i]) > std::abs(r[i]))
                    r[i] = a[i];
        if (reduce_type == reduce::abs_min)
            for(int i = 0; i < n; i++)
                if (std::abs(a[i]) < std::abs(r[i]))
                    r[i] = a[i];
        if (reduce_type == reduce::abs_min_first) {
            // printf("p=%d:  %8.3e, %8.3e, %8.3e, %8.3e\n", p, a[0], a[1], r[0], r[1]);
            if ((std::abs(a[0]) < std::abs(r[0]) || r[1] < -0.5) && a[1] > -0.5) {
                // printf("Yep!\n");
                for(int i = 0; i < n; i++) {
                    r[i] = a[i];
                }
            }
        }
        if (reduce_type == reduce::abs_max_first) {
            // printf("p=%d:  %8.3e, %8.3e, %8.3e, %8.3e\n", p, a[0], a[1], r[0], r[1]);
            if ((std::abs(a[0]) > std::abs(r[0]) || r[1] < -0.5) && a[1] > -0.5) {
                // printf("Yep!\n");
                for(int i = 0; i < n; i++) {
                    r[i] = a[i];
                }
            }
        }
    }
	t_in++;
	if (t_in >= p) {
		t_out = 0;
		pthread_cond_broadcast(&c_in);
	}
	else {
		while (t_in < p) {
			pthread_cond_wait(&c_in, &m);
		}
	}
	if (r != a){
		for (int i = 0; i < n; i++) {
			a[i] = r[i];
		}
	}
	t_out++;
	if (t_out >= p) {
		t_in = 0;
		r = nullptr;
		pthread_cond_broadcast(&c_out);
	}
	else {
		while (t_out < p){
			pthread_cond_wait(&c_out, &m);
		}
	}
	pthread_mutex_unlock(&m);

}

void get_block(double *a, double *block, size_t n, size_t m, size_t k, size_t l,
               size_t i, size_t j)
{
    size_t row_count, col_count;
    size_t row_offset, col_offset;
    size_t block_row, block_col;

    // Определяем размеры блока
    if (i < k)
        row_count = m;
    else
        row_count = l;

    if (j < k)
        col_count = m;
    else
        col_count = l;

    // Вычисляем смещение в исходной матрице
    row_offset = i * m;
    col_offset = j * m;

    // Копируем элементы блока
    for (block_row = 0; block_row < row_count; block_row++)
    {
        for (block_col = 0; block_col < col_count; block_col++)
        {
            block[block_row * col_count + block_col] =
                a[(row_offset + block_row) * n + (col_offset + block_col)];
        }
    }
}

void put_block(double *a, double *block, size_t n, size_t m, size_t k, size_t l,
               size_t i, size_t j)
{
    size_t row_size, col_size;
    size_t row_offset, col_offset;
    size_t block_row, block_col;

    // Определяем размеры блока
    if (i < k)
        row_size = m;
    else
        row_size = l;

    if (j < k)
        col_size = m;
    else
        col_size = l;

    // Вычисляем смещение в матрице
    row_offset = i * m;
    col_offset = j * m;

    // Копируем элементы из блока обратно в матрицу
    for (block_row = 0; block_row < row_size; block_row++)
    {
        for (block_col = 0; block_col < col_size; block_col++)
        {
            a[(row_offset + block_row) * n + (col_offset + block_col)] =
                block[block_row * col_size + block_col];
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

int fill_matrix(double *matrix, size_t n, size_t s)
{
    size_t i, j;
    switch (s) {
    case 1:
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                matrix[i * n + j] = n - MAX(i, j);
        break;
    case 2:
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                matrix[i * n + j] = MAX(i, j) + 1;
        break;
    case 3:
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                matrix[i * n + j] = (i > j ? i - j : j - i);
        break;
    case 4:
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                matrix[i * n + j] = 1. / (i + j + 1);
        break;
    default:
        fprintf(stderr, "[-] Unknown formula %ld\n", s);
        return -1;
    }
    return 0;
}

int read_matrix_from_file(double *matrix, size_t n, FILE* file)
{
    rewind(file);
    // printf("size is %ld\n", n * n);
    for (size_t i = 0; i < n * n; i++) {
        if (fscanf(file, "%lf", &matrix[i]) != 1) {
            fprintf(stderr, "[-] Can't read file.\n");
            return -1;
        }
    }
    return 0;
}

void unit_matrix(double *matrix, size_t n)
{
    memset(matrix, 0, n * n * sizeof(double));
    for (size_t i = 0; i < n; i++) {
        matrix[i * n + i] = 1;
    }
}

void zero_matrix(double *matrix, size_t n, size_t m)
{
    memset(matrix, 0, n * m * sizeof(double));
}

void matrix_multiply(const double *A, const double *B, double *C, size_t p,
                     size_t q, size_t r) {
    size_t i, j, k;
    for (i = 0; i < p; ++i) {
        for (k = 0; k < r; ++k) {
            C[i * r + k] = 0.0;
        }
    }


    for (i = 0; i < p; ++i) {
        for (j = 0; j < q; ++j) {

            double a_ij = A[i * q + j];
            k = 0;
            for (; k + 3 < r; k += 4) {
                C[i * r + k]     += a_ij * B[j * r + k];
                C[i * r + k + 1] += a_ij * B[j * r + k + 1];
                C[i * r + k + 2] += a_ij * B[j * r + k + 2];
                C[i * r + k + 3] += a_ij * B[j * r + k + 3];
            }

            for (; k < r; ++k) {
                C[i * r + k] += a_ij * B[j * r + k];
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
    size_t i;
    for (i = 0; i < n * m; i++) {
        C[i] = A[i] + B[i];
    }
}

void matrix_subtr(double *A, double *B, size_t n, size_t m) // A -= B
{
    size_t i;
    for (i = 0; i < n * m; i++) {
        A[i] -= B[i];
    }
}

double get_norm(double *matrix, size_t m)
{
    size_t i, j;
    double max = 0, sum;
    for (j = 0; j < m; j++) {
        sum = 0;
        for (i = 0; i < m; i++) {
            sum += std::abs(matrix[j + i * m]);
        }
        // printf("sum is %8.3e\n", sum);
        max = std::max(max, sum);
    }
    return max;
}

int get_inverse_matrix(double *A, double *B, size_t m)
{
    size_t i, j, k, max_ind;
    double max, buf;
    unit_matrix(B, m);
    for (j = 0; j < m; j++) {
        // print_matrix(A, m, m);
        // printf("\n");
        max = std::abs(A[j + j * m]);
        max_ind = j;
        for (i = j + 1; i < m; i++) {
            if (std::abs(A[j + i * m]) > max) {
                max = std::abs(A[j + i * m]);
                max_ind = i;
            }
        }
        // printf("max: %lf, max_ind: %ld\n", max, max_ind);

        if (max < EPS) {
            return -1;
        }
        // print_matrix(A, m, m);
        // printf("^^^A0\n");
        // print_matrix(B, m, m);
        // printf("^^^B0\n");
        // printf("max_ind is %ld\n", max_ind);

        // Меняем местами строки, если нашли нужную
        if (max_ind != j) {
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
        }
        
        max = A[j + j * m];
        A[j + j * m] = 1;
        for (i = j + 1; i < m; i++) {
            A[i + j * m] /= max;
        }

        for (i = 0; i < m; i++) {
            B[i + j * m] /= max;
        }

        // print_matrix(A, m, m);
        // printf("^^^A\n");
        // print_matrix(B, m, m);
        // printf("^^^B\n");
        for (i = 0; i < m; i++) {
            if (i != j) {
                max = A[j + i * m];
                A[j + i * m] = 0;
                for (k = j + 1; k < m; k++) {
                    A[k + i * m] -= max * A[k + j * m];
                }
                for (k = 0; k < m; k++) {
                    B[k + i * m] -= max * B[k + j * m];
                }
                // if (std::abs(max) > EPS) {
                // }
            }
        }
        // print_matrix(A, m, m);
        // printf("^^^A2\n\n");
        // print_matrix(B, m, m);
        // printf("^^^B2\n\n");

    }
    return 0;
}

double mult_sub_norm_p(double *a, double *b, double *pc, double *norm, size_t n, size_t m, size_t p, size_t pi)
{
    size_t i, j, s, r, t, q;
    size_t k = n / m;
    size_t l = n - k * m; // n = k * m + l
    size_t bl = (l != 0 ? k + 1 : k); // Общее количество блоков
    size_t v, h, ah;
    double max_norm = 0;

    for (i = 0; i < n; i++) norm[i] = 0;

    // Проходим по всем блокам в матрице C
    for (i = 0; i < bl; i++) {
        for (j = pi; j < bl; j += p) {
            // Определяем размер текущего блока C[v x h]
            v = (i < k ? m : l); // вертикальный размер блока
            h = (j < k ? m : l); // горизонтальный размер блока

            // Указатель на начало текущего блока C
            // double *pc = c + (i * m) * n + j * m;

            // Инициализируем блок C нулями
            for (r = 0; r < v; r++) {
                for (t = 0; t < h; t++) {
                    pc[r * h + t] = 0.0;
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

                        pc[(r + 0) * h + (t + 0)] += s00;
                        pc[(r + 0) * h + (t + 1)] += s01;
                        pc[(r + 0) * h + (t + 2)] += s02;

                        pc[(r + 1) * h + (t + 0)] += s10;
                        pc[(r + 1) * h + (t + 1)] += s11;
                        pc[(r + 1) * h + (t + 2)] += s12;

                        pc[(r + 2) * h + (t + 0)] += s20;
                        pc[(r + 2) * h + (t + 1)] += s21;
                        pc[(r + 2) * h + (t + 2)] += s22;

                        // if (r == t && i == j) {
                        //     max_norm = MAX(std::abs(1. - s00), max_norm);
                        //     max_norm = MAX(std::abs(1. - s11), max_norm);
                        //     max_norm = MAX(std::abs(1. - s22), max_norm);
                        // } else {
                        //     max_norm = MAX(std::abs(s00), max_norm);
                        //     max_norm = MAX(std::abs(s11), max_norm);
                        //     max_norm = MAX(std::abs(s22), max_norm);
                        // }

                        // max_norm = MAX(std::abs(s01), max_norm);
                        // max_norm = MAX(std::abs(s02), max_norm);
                        // max_norm = MAX(std::abs(s10), max_norm);
                        // max_norm = MAX(std::abs(s12), max_norm);
                        // max_norm = MAX(std::abs(s20), max_norm);
                        // max_norm = MAX(std::abs(s21), max_norm);

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

                        pc[(r + 0) * h + t] += s0;
                        pc[(r + 1) * h + t] += s1;
                        pc[(r + 2) * h + t] += s2;
                        // max_norm = MAX(std::abs(s0), max_norm);
                        // max_norm = MAX(std::abs(s1), max_norm);
                        // max_norm = MAX(std::abs(s2), max_norm);
                        
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

                        pc[r * h + (t + 0)] += s0;
                        pc[r * h + (t + 1)] += s1;
                        pc[r * h + (t + 2)] += s2;

                        // max_norm = MAX(std::abs(s0), max_norm);
                        // max_norm = MAX(std::abs(s1), max_norm);
                        // max_norm = MAX(std::abs(s2), max_norm);

                    }

                    // Обработка оставшихся элементов
                    for (t = t_end; t < h; t++) {
                        double sum = 0.0;

                        for (q = 0; q < ah; q++) {
                            sum += pa[r * n + q] * pb[q * n + t];
                        }

                        pc[r * h + t] += sum;
                        // if (r == t && i == j) max_norm = MAX(std::abs(1 - sum), max_norm);
                        // else max_norm = MAX(std::abs(sum), max_norm);
                    }
                }
            }
            for (r = 0; r < v; r++) {
                for (t = 0; t < h; t++) {
                    if (i == j && r == t) {
                        norm[t + m * j] += std::abs(1 - pc[r * h + t]);
                    } else {
                        norm[t + m * j] += std::abs(pc[r * h + t]);
                    }
                }
            }
        }
    }

    // synchronize(p, &max_norm, 1, reduce::abs_max);
    synchronize(p);
    for (i = 0; i < n; i++) max_norm = MAX(max_norm, norm[i]);
    synchronize(p);


    return max_norm;
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

int find_diff(double *matrix, double *inversed_matrix, double* block, double *norm, FILE *file, int n, int m, int s, double &r1, double &r2, size_t p, size_t pi)
{
    if (n < 11000) {
        if (s == 0 && pi == 0) {
            if (read_matrix_from_file(matrix, n, file) != 0) {
                return 2;
            }
        } else if (s != 0) {
            if (fill_matrix(matrix, n, s) != 0) {
                return 2;
            }
        }
        synchronize(p);

        r1 = mult_sub_norm_p(matrix, inversed_matrix, block, norm, n, m, p, pi);
        r2 = mult_sub_norm_p(inversed_matrix, matrix, block, norm, n, m, p, pi);
        // printf("F pi: %ld r1 = %8.3e, r2 = %8.3e\n", pi, r1, r2);

    } else {
        r1 = 0;
        r2 = 0;
    }
    return 0;
}


void zero_matrix_p(double *matrix, size_t n, size_t m, size_t p, size_t pi)
{
    size_t l = n % m, k = n / m, bl = (l > 0 ? k + 1 : k);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = pi; j < bl; j += p) {
            memset(matrix + i * n + j * m, 0, (j == k ? l : m) * sizeof(double));
        }
    }
}


void unit_matrix_p(double *matrix, size_t n, size_t m, size_t p, size_t pi)
{
    zero_matrix_p(matrix, n, m, p, pi);
    // printf("pi: %ld\n", pi);
    synchronize(p);

    // if (pi == 0) {
    //     for (size_t i = 0; i < n; i++) {
    //         matrix[i * n + i] = 1;
    //     }
    // }
    for (size_t i = pi * m; i < n; i += p * m) {
        for (size_t j = i; j < i + m && j < n; j++) {
            matrix[j * n + j] = 1;
        }
    }
}


int fill_matrix_p(double *matrix, size_t n, size_t s, size_t m, size_t p, size_t pi)
{
    size_t i, j, k;
    switch (s) {
    case 1:
        for (i = 0; i < n; i++)
            for (j = pi * m; j < n; j+= p * m)
                for (k = j; k < j + m && k < n; k++)
                    matrix[i * n + k] = n - MAX(i, k);
        break;
    case 2:
        for (i = 0; i < n; i++)
            for (j = pi * m; j < n; j+= p * m)
                for (k = j; k < j + m && k < n; k++)
                    matrix[i * n + k] = MAX(i, k) + 1;
        break;
    case 3:
        for (i = 0; i < n; i++)
            for (j = pi * m; j < n; j+= p * m)
                for (k = j; k < j + m && k < n; k++)
                    matrix[i * n + k] = (i > k ? i - k : k - i);
        break;
    case 4:
        for (i = 0; i < n; i++)
            for (j = pi * m; j < n; j+= p * m)
                for (k = j; k < j + m && k < n; k++)
                    matrix[i * n + k] = 1. / (i + k + 1);
        break;
    default:
        fprintf(stderr, "[-] Unknown formula %ld\n", s);
        return -1;
    }
    return 0;
}


double get_norm_p(double *matrix, size_t n, size_t m, size_t p, size_t pi)
{
    size_t i, j, k;
    double max = 0, sum;
    for (j = pi * m; j < n; j += p * m) {
        for (k = j; k < j + m && k < n; k++) {
            sum = 0;
            for (i = 0; i < n; i++) {
                sum += std::abs(matrix[k + i * n]);
            }
            max = std::max(max, sum);
            // printf("sum is %8.3e\n", sum);
        }
    }
    return max;
}


int process_args(Args *a)
{
    int error = 0;
	for (size_t i = 0; i < a->p; i++) {
		if (a[i].error_type != io_status::success) {
			printf("[-] Error in thread %s\n    ", a[i].name);
            if (a[i].error_type == io_status::bad_allocation) printf("Can't allocate memory!\n");
            if (a[i].error_type == io_status::error_open) printf("Can't open the file!\n");
            if (a[i].error_type == io_status::error_read) printf("Bad file!\n");
            if (a[i].error_type == io_status::irreversible) printf("Irreversible matrix!\n");
            if (a[i].error_type == io_status::unknown_formula) printf("Unknown formula!\n");
			error++;
		}
	}
	return error;
}


void rows_permutation_p(double *A, double *block1, double *block2, size_t n,
                      size_t m, size_t k, size_t l, size_t i1, size_t i2,
                      size_t begin, size_t p, size_t pi)
{
    size_t i;
    if (begin % p <= pi) begin += pi - begin % p;
    else begin += p + pi - begin % p;

    for (i = begin; i < k + 1; i += p) {
        get_block(A, block1, n, m, k, l, i1, i);
        get_block(A, block2, n, m, k, l, i2, i);
        put_block(A, block2, n, m, k, l, i1, i);
        put_block(A, block1, n, m, k, l, i2, i);
    }
}

double get_time() {
    struct timeval buf;
    gettimeofday(&buf, 0);
    return buf.tv_sec + buf.tv_usec / 1e6;
}


void *thread_func(void *args)
{
	Args *a = (Args*) args;
	
	size_t p = a->p, pi = a->pi, n = a->n, m = a->m, s = a->s, k = a->k, l = a->l, r = a->r;
    size_t diag, i, j, bl = (l == 0) ? k : k + 1, min_norm_ind, row, x, y, z, begin;
    double *matrix = a->matrix, *inversed_matrix = a->inversed_matrix, *block = a->block, *norm_arr = a->norm;
    // double r1 = a->r1, r2 = a->r2;
    double norm = 0, min_norm, buf_array[2];
    FILE *file = a->file;
    pthread_t tid = a->tid;
    cpu_set_t cpu;
	int n_cpus = get_nprocs(); // число процессоров
	int cpu_id = n_cpus - 1 - (pi % n_cpus);

	CPU_ZERO(&cpu); // вместо конструктора
    CPU_SET(cpu_id, &cpu);

    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);

    double *block_A = new double[m*m];
    double *block_B = new double[m*m];
    double *block_C = new double[m*m];

    if (!block_A || !block_B || !block_C) {
        if (block_A) delete[] block_A;
        if (block_B) delete[] block_B;
        if (block_C) delete[] block_C;
        a->error_flag = 1;
        a->error_type = io_status::bad_allocation;
    }
	synchronize(p, & a->error_flag, 1);
	if (a -> error_flag > 0) {
        delete[] block_A;
        delete[] block_B;
        delete[] block_C;
		return nullptr;
	}

    a->t1 = get_time();

    zero_matrix_p(matrix, n, m, p, pi);
    unit_matrix_p(inversed_matrix, n, m, p, pi);
    synchronize(p);


    if (s == 0) {
        if (pi == 0) {
            if (read_matrix_from_file(matrix, n, file) != 0) {
                a->error_flag = 1;
                a->error_type = io_status::error_read;
            }
        }
    } else {
        if (fill_matrix_p(matrix, n, s, m, p, pi) != 0) {
            a->error_flag = 1;
            a->error_type = io_status::unknown_formula;
        }
    }

    if (pi == 0) {
        printf("[+] Given matrix:\n");
        print_matrix(matrix, n, 4);
    }
	// printf("Yep!\n");
	
	synchronize(a->p, & a->error_flag, 1);
	if (a -> error_flag > 0) {
        delete[] block_A;
        delete[] block_B;
        delete[] block_C;
		return nullptr;
	}

    norm = get_norm_p(matrix, n, m, p, pi);
	synchronize(a->p, &norm, 1);
	if (pi == 0) EPS *= norm;
    // printf("norm = %8.3e, eps = %8.3e\n", norm, EPS);
	synchronize(a->p);

    // Пробегаюсь по диагональным элементам
    for (diag = 0; diag < bl; diag++) {
        // print_matrix(matrix, n, n);
        // printf("\n");
        min_norm = -1;
        min_norm_ind = -1;
        // Ищу в столбце матрицу у кооторой обратная имеет наименьшую норму
        if (diag < k) {
            for (row = diag + pi; row < k; row += p) {
                get_block(matrix, block_A, n, m, k, l, diag, row);
                if (get_inverse_matrix(block_A, block_B, m) == 0) {
                    // print_matrix(block_B, m, print_size);
                    norm = get_norm(block_B, m);
                    // printf("norm is %8.3e\n", norm);
                    if (min_norm < 0 || norm < min_norm) {
                        min_norm = norm;
                        min_norm_ind = row;
                    }
                }
            }
        } else {
            if (pi == 0) {
                get_block(matrix, block_A, n, m, k, l, k, k);
                if (get_inverse_matrix(block_A, block_B, l) == 0) {
                    norm = get_norm(block_B, l);
                    min_norm = norm;
                    min_norm_ind = k;
                }
            }
        }
        buf_array[0] = min_norm;
        buf_array[1] = min_norm_ind * 1. > 1e+13 ? -1. : min_norm_ind;
        // buf_array[1] = buf_array[1] * 1. > 1e+13 ? -1. : buf_array[1];

        
        // printf("%ld min_norm: %8.3e, ind: %.0e ---buf\n", pi, buf_array[0], buf_array[1]);
        // printf("%ld min_norm: %8.3e, ind: %lu\n", pi, min_norm, min_norm_ind);
        synchronize(p, buf_array, 2, reduce::abs_min_first);
        // if (pi == 0) printf("%ld min_norm: %8.3e, ind: %.0e\n", pi, buf_array[0], buf_array[1]);
        // printf("Yep!\n");
        if (buf_array[1] < 0) {

            fprintf(stderr, "[-] Matrix is irreversible, place: 923!\n");
            a->error_flag = 1;
            a->error_type = io_status::irreversible;
            delete[] block_A;
            delete[] block_B;
            delete[] block_C;
            return nullptr;
        }
        // min_norm = buf_array[0];
        min_norm_ind = buf_array[1];
        // printf("Step: %ld, thread: %ld\n", diag, pi);

        // Переставляю строки так, чтобы на текущем диагональном элементе была
        // матрица с наименьшей нормой её обратной
        rows_permutation_p(
            matrix, block_A, block_B, n, m, k, l, min_norm_ind, diag, diag, p, pi
        );
        rows_permutation_p(
            inversed_matrix, block_A, block_B, n, m, k, l, min_norm_ind, diag, 0, p, pi
        );
        
        // Нахожу обратную в потоке pi
        if (diag % p == pi) {
            get_block(matrix, block_A, n, m, k, l, diag, diag);
            if (diag != k) {
                if (get_inverse_matrix(block_A, block, m) != 0) {
                    fprintf(stderr, "[-] Matrix is irreversible, place 949!\n");
                    a->error_flag = 1;
                    a->error_type = io_status::irreversible;
                    delete[] block_A;
                    delete[] block_B;
                    delete[] block_C;
                    return nullptr;
                }
            } else {
                if (get_inverse_matrix(block_A, block, l) != 0) {
                    fprintf(stderr, "[-] Matrix is irreversible!\n");
                    a->error_flag = 1;
                    a->error_type = io_status::irreversible;
                    delete[] block_A;
                    delete[] block_B;
                    delete[] block_C;
                    return nullptr;
                }
            }
        }
        synchronize(p);

        memcpy(block_B, block, m * m * sizeof(double));
        
        // Вставляю единичную на место текущего диагонального элемента
        // put_block(matrix, block_A, n, m, k, l, diag, diag);

        // Начиная со следующего каждый элемент в строке домножаю на обратную к диагональному
        // block_B - обратная к диагональному
        begin = diag + 1;
        if (begin % p <= pi) begin += pi - begin % p;
        else begin += pi + p - begin % p;

        for (i = begin; i < bl; i+=p) {
            get_block(matrix, block_A, n, m, k, l, diag, i);
            x = (diag == k ? l : m);
            y = (i == k    ? l : m);
            matrix_multiply(block_B, block_A, block_C, x, x, y);
            put_block(matrix, block_C, n, m, k, l, diag, i);
        }

        // Каждый элемент той же строки матрицы В домножаю на эту же обратную
        for (i = pi; i < bl; i+=p) {
            get_block(inversed_matrix, block_A, n, m, k, l, diag, i);
            // printf("------matrix------\n");
            // printf("Diag = %ld, i = %ld\n", diag, i);
            // print_matrix(block_A, (i < k ? m : l), 4);
            x = (diag == k ? l : m);
            y = (i == k    ? l : m);
            // printf("x = %ld, y = %ld\n", x, y);
            // printf("------------------\n");
            matrix_multiply(block_B, block_A, block_C, x, x, y);
            put_block(inversed_matrix, block_C, n, m, k, l, diag, i);
        }

        synchronize(p);

        // Для каждой строки, кроме той, на которой есть текущий диагональный
        for (i = 0; i < bl; i++) {
            if (i != diag) {
                // Запоминаем блок диагонального столбца и текущей строки

                get_block(matrix, block_A, n, m, k, l, i, diag);
                // Для каждого следующего элемента текущей строки, вычитаем из него А х С
                for (j = begin; j < bl; j+=p) { // j - номер столбца

                    get_block(matrix, block_C, n, m, k, l, diag, j);

                    x = (i == k    ? l : m);
                    y = (diag == k ? l : m);
                    z = (j == k    ? l : m);
                    matrix_multiply(block_A, block_C, block_B, x, y, z);

                    get_block(matrix, block_C, n, m, k, l, i, j);
                    matrix_subtr(block_C, block_B, x, z);
                    put_block(matrix, block_C, n, m, k, l, i, j);
                }

                // Аналогично для матрицы В
                for (j = pi; j < bl; j+=p) {
                    get_block(inversed_matrix, block_C, n, m, k, l, diag, j);
                    x = (i == k    ? l : m);
                    y = (diag == k ? l : m);
                    z = (j == k    ? l : m);
                    matrix_multiply(block_A, block_C, block_B, x, y, z);

                    get_block(inversed_matrix, block_C, n, m, k, l, i, j);
                    matrix_subtr(block_C, block_B, x, z);
                    put_block(inversed_matrix, block_C, n, m, k, l, i, j);
                }
            }
        }
        synchronize(p);
    }

	synchronize(p, & a-> error_flag, 1);
	if (a -> error_flag > 0) {
		delete[] block_A;
        delete[] block_B;
        delete[] block_C;
		return nullptr;
	}

    a->t1 -= get_time();

    if (pi == 0) {
        printf("\n[+] Inversed matrix:\n");
        print_matrix(inversed_matrix, n, r);
    }


    a->t2 = get_time();

    if (find_diff(matrix, inversed_matrix, block_A, norm_arr, file, n, m, s, a->r1, a->r2, p, pi) != 0) {
        a->error_type = io_status::error_read;
        a->error_flag = 1;
    }
	synchronize(p, & a-> error_flag, 1);
    if (a -> error_flag > 0) {
		delete[] block_A;
        delete[] block_B;
        delete[] block_C;
		return nullptr;
	}

    a->t2 -= get_time();


	a -> error_type = io_status::success;

    delete[] block_A;
    delete[] block_B;
    delete[] block_C;
	return nullptr;
}
