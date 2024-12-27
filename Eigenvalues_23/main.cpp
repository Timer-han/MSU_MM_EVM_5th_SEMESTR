#include <cstdio>

#include "functions.h"

int main(int argc, char *argv[])
{
    if (argc < 5) {
        printf("[-] Not enough args!\n");
        return -1;
    }

    int n, m, k, its = 0;
    double t1, t2, res1, res2, eps;

    if (sscanf(argv[1], "%d", &n) != 1 || sscanf(argv[2], "%d", &m) != 1 ||
        sscanf(argv[3], "%le", &eps) != 1 || sscanf(argv[4], "%d", &k) != 1) {
        printf("[-] Bad args!\n");
        return -2;
    }

    if (k == 0 && argc < 6) {
        printf("[-] Can't find filename!\n");
        return -3;
    }

    char *filename = argv[5];
    matrix_norm = matrix_norm;
    matrix_trace = matrix_trace;
    matrix_lengh = matrix_lengh;

    double *matrix = new double[n * n];
    double *vector1 = new double[n];
    double *vector2 = new double[n];
    double *vector3 = new double[n];
    double *eigenvals = new double[n];

    if (!matrix || !vector1 || !vector2 || !vector3) {
        printf("[-] Can't allocate memory!\n");
        if (matrix) delete[] matrix;
        if (vector1) delete[] vector1;
        if (vector2) delete[] vector2;
        if (vector3) delete[] vector3;
        if (eigenvals) delete[] eigenvals;
        return -4;
    }

    if (k == 0 && read_matrix(filename, matrix, n) != 0) {
        printf("[-] Bad file!\n");
        delete[] matrix;
        delete[] vector1;
        delete[] vector2;
        delete[] vector3;
        delete[] eigenvals;
        return -5;
    } else if (k > 0 && k < 5) {
        fill_matrix(matrix, n, k);
    } else if (k < 0 || n > 4) {
        printf("[-] Unknown formula %d!\n", k);
        return -6;
    }

    printf("Given matrix:\n");
    print_matrix(matrix, n, 4);

    if (
        run(matrix, vector1, vector2, vector3, eigenvals, n, m, t1, t2, res1, res2, its, eps) < 0
    )
    {
        delete[] matrix;
        delete[] vector1;
        delete[] vector2;
        delete[] vector3;
        delete[] eigenvals;
    }
    printf("RESULT:");
    for (int i = 0; i < m; i++) printf("% 8.3e", eigenvals[i]);
    printf("\n");

    printf(
        "%s : Residual1 = %e Residual2 = %e Iterations = %d \
Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n",
        argv[0], res1, res2, its, its / n, t1, t2
    );
    delete[] matrix;
    delete[] vector1;
    delete[] vector2;
    delete[] vector3;
    delete[] eigenvals;
    return 0;
}