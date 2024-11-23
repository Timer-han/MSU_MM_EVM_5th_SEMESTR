#include "functions.h"

int main(int argc, char *argv[])
{
    double r1 = 0, r2 = 0;
    int z = 0;

    if (argc < 5) {
        fprintf(stderr, "[-] Not enough arguments.\n");
        return 1;
    }

    size_t n = 0, m = 0, r = 0, s = 0;
    if (
        sscanf(argv[1], "%lu", &n) != 1 ||
        sscanf(argv[2], "%lu", &m) != 1 ||
        sscanf(argv[3], "%lu", &r) != 1 ||
        sscanf(argv[4], "%lu", &s) != 1
    )
    {
        printf("[-] Mistake in args!\n");
        return 2;
    }
    size_t l = n % m;
    size_t k = n / m;

    char *filename = nullptr;
    if (s == 0) {
        if (argc < 6) {
            fprintf(stderr, "[-] File name do not defined.\n");
            return 1;
        }
        filename = argv[5];
    }

    double *matrix = new double[n*n];
    double *inversed_matrix =  new double[n*n];
    double *block_A = new double[m*m];
    double *block_B = new double[m*m];
    double *block_C = new double[m*m];
    double *norm = new double[n];


    double t1 = clock();
    if (!block_C || !block_B || !block_C || !matrix || !inversed_matrix || !norm ||
        (z = run(matrix, inversed_matrix, block_A, block_B, block_C, n, m, k, l, s, r, filename))) {
        if (block_A) delete[] block_A;
        if (block_B) delete[] block_B;
        if (block_C) delete[] block_C;
        if (matrix) delete[] matrix;
        if (inversed_matrix) delete[] inversed_matrix;
        if (norm) delete[] norm;
        if (!block_C || !block_B || !block_C || !matrix || !inversed_matrix) {
            fprintf(stderr, "[-] Can't allocate memory for matrix\n");
        }
        return z;
    }

    t1 = clock() - t1;


    double t2 = clock();
    if ((z = find_diff(matrix, inversed_matrix, block_A, norm, filename, n, m, s, r1, r2)) != 0) {
        return z;
    } else {
        t2 = clock() - t2;
        printf("[+] Inversed matrix:\n");
        print_matrix(inversed_matrix, n, r);
        printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %ld N = "
            "%ld M = %ld\n",
            argv[0], 18, r1, r2, t1 / CLOCKS_PER_SEC,
            t2 / CLOCKS_PER_SEC, s, n, m);
    }



    delete[] block_A;
    delete[] block_B;
    delete[] block_C;
    delete[] matrix;
    delete[] inversed_matrix;
    delete[] norm;
    return 0;
}
