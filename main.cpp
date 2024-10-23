#include "functions.h"

int main(int argc, char *argv[])
{
    double r1 = 0, r2 = 0;

    if (argc < 5) {
        fprintf(stderr, "Not enough arguments.\n");
        return 1;
    }

    size_t n = atoi(argv[1]);
    size_t m = atoi(argv[2]);
    size_t r = atoi(argv[3]);
    size_t s = atoi(argv[4]);
    size_t l = n % m;
    size_t k = n / m;

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
    find_diff(matrix, inversed_matrix, filename, n, m, s, r1, r2);
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
