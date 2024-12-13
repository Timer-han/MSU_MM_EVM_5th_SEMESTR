#include "functions.h"

int main(int argc, char *argv[])
{
    double r1 = 0, r2 = 0;
    FILE *file = nullptr;

    if (argc < 6) {
        fprintf(stderr, "[-] Not enough arguments.\n");
        return 1;
    }

    size_t n = 0, m = 0, p = 0, r = 0, s = 0;
    if (
        sscanf(argv[1], "%lu", &n) != 1 ||
        sscanf(argv[2], "%lu", &m) != 1 ||
        sscanf(argv[3], "%lu", &p) != 1 ||
        sscanf(argv[4], "%lu", &r) != 1 ||
        sscanf(argv[5], "%lu", &s) != 1
    ) {
        printf("[-] Mistake in args!\n");
        return 2;
    }

    size_t l = n % m;
    size_t k = n / m;

    if (l > 0 && k + 1 < p) p = k + 1;
    if (l == 0 && k < p) p = k;

    char *filename = nullptr;
    if (s == 0) {
        if (argc < 7) {
            fprintf(stderr, "[-] File name do not defined.\n");
            return 1;
        }
        filename = argv[6];

        file = fopen(filename, "r");
        if (!file) {
            printf("[-] Can't open file \"%s\"\n", filename);
            return -2;
        }
    }


    double *matrix  = new double[n*n];
    double *inversed_matrix =  new double[n*n];
    double *block_A =  new double[m*m];
    double *norm    = new double[n];
    Args   *a       = new Args[p];

	if (!matrix || !inversed_matrix || !norm || !a) {
        printf("[-] Not enough memory!\n");
        if (matrix)  delete[] matrix;
        if (inversed_matrix) delete[] inversed_matrix;
        if (block_A)    delete[] block_A;
        if (norm)    delete[] norm;
        if (a)       delete[] a;
        if (file) fclose(file);
        return -1;
    }
	int error;
    size_t i;

	for (i = 0; i < p; i++) {
        a[i].n = n;
        a[i].k = k;
        a[i].l = l;
        a[i].m = m;
        a[i].s = s;
		a[i].p = p;
		a[i].pi = i;
		a[i].file = file;
        a[i].block = block_A;
        a[i].matrix = matrix;
        a[i].inversed_matrix = inversed_matrix;

	}

    double t1 = clock();
	for (i = 1; i < p; i++) {
		if ((error = pthread_create(&a[i].tid, nullptr, thread_func, a + i))) {
            delete[] matrix;
            delete[] inversed_matrix;
            delete[] block_A;
            delete[] norm;
            delete[] a;
            if (file) fclose(file);
			return error;
		}
	}
    a[0].tid = pthread_self();
    thread_func(a);

	for (i = 1; i < p; i++) {
		pthread_join(a[i].tid, nullptr);
	}
	if ((error = process_args(a))) {
		delete[] matrix;
        delete[] inversed_matrix;
        delete[] block_A;
        delete[] norm;
        delete[] a;
        if (file) fclose(file);
        return error;
	}
    // printf("Yep!\n");
    t1 = clock() - t1;
    printf("\n[+] Inversed matrix:\n");
    print_matrix(inversed_matrix, n, r);


    // if ((z = run(matrix, inversed_matrix, block_A, block_B, block_C, n, m, k, l, s, r, filename))) {
    //     if (block_A) delete[] block_A;
    //     if (block_B) delete[] block_B;
    //     if (block_C) delete[] block_C;
    //     if (matrix) delete[] matrix;
    //     if (inversed_matrix) delete[] inversed_matrix;
    //     if (norm) delete[] norm;
    //     if (!block_A || !block_B || !block_C || !matrix || !inversed_matrix) {
    //         fprintf(stderr, "[-] Can't allocate memory for matrix\n");
    //     }
    //     return z;
    // }


    
    double t2 = clock();
    if ((error = find_diff(matrix, inversed_matrix, block_A, norm, file, n, m, s, r1, r2)) != 0) {
        delete[] matrix;
        delete[] inversed_matrix;
        delete[] block_A;
        delete[] norm;
        delete[] a;
        if (file) fclose(file);
        return error;
    } else {
        t2 = clock() - t2;
        // printf("[+] Inversed matrix:\n");
        // print_matrix(inversed_matrix, n, r);
        printf("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %ld N = "
            "%ld M = %ld P = %ld\n",
            argv[0], 18, r1, r2, t1 / CLOCKS_PER_SEC,
            t2 / CLOCKS_PER_SEC, s, n, m, p);
    }



    delete[] matrix;
    delete[] inversed_matrix;
    delete[] block_A;
    delete[] norm;
    delete[] a;
    if (file) fclose(file);
    return 0;
}
