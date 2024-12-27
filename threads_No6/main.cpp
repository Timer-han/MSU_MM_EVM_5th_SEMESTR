#include <cstdio>
#include <iostream>
#include <pthread.h>

enum class io_status
{
	undef,
	error_open,
	error_read,
	success
};

class Args
{
public:
	int p = 0; // количество потоков
	int pi = 0; // номер потока
	pthread_t tid = 1;
	const char * name = nullptr;

	int n = 0;
	int m = 0;
	int r = 0; // количество строк для потока
	int prev_r = 0;
	double count = 0;
	double sum = 0;
	double *matrix = nullptr;
	double *part = nullptr;
	double *prev = nullptr;
	double *next = nullptr;
	double res = 0; // количество элементов, больших предыдущего
	io_status error_type = io_status::undef;
	double error_flag = 0;
};

void synchronize(int p, double *a, int n);
io_status nums_counter(FILE * f, double* count);
void *thread_func(void *args);
int process_args(Args *a);

void synchronize(int p, double *a=nullptr, int n = 0) {
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
	static int t_in = 0;//количество вошедших потоков
	static int t_out = 0;//количество вышедших потоков
	static double *r = nullptr;
	int i;
	if (p <= 1) return;
	pthread_mutex_lock(&m);
	if (r == nullptr) r = a;
	else for(i = 0; i < n; i++) r[i] += a[i];
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
		for (i = 0; i < n; i++) {
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

io_status nums_counter(FILE * f, double* count)
{
    int a, b, c;
    *count = 0;
    
    if (fscanf(f, "%d", &b) != 1 && !feof(f)) {
        return io_status::error_read;
    }
    if (fscanf(f, "%d", &c) != 1 && !feof(f)) {
        return io_status::error_read;
    }
    for (;;) {
        a = b;
        b = c;
        if (fscanf(f, "%d", &c) != 1) {
            break;
        }
        if (a < b && a < c) {
            count[0]++;
        }
    }
    if (feof(f)) {
        return io_status::success;
    } else {
        return io_status::error_read;
    }
}

void *thread_func(void *args)
{
	Args *a = (Args*) args;
	int r = a->r;
	int prev_r = a->prev_r - 1;
	int m = a->m;
	int count = 0;
	double sum = 0;
	double *part = a->part;
	double *prev = a->prev;
	double *next = a->next;
	
	a->error_flag = 0;
	synchronize(a->p, & a->error_flag, 1);

	if (r == 1) {
		if (
			part[0] > part[1] && 
			(!prev || part[0] > prev[0]) && 
			(!next || part[0] > next[0])
		)
		{
			printf("1 x = %d, y = %d\n", 0, (int) (part - a->matrix) / m);
			sum += part[0];
			count++;
		}
		if (
			part[m - 1] > part[m - 2] && 
			(!prev || part[m - 1] > prev[prev_r * m + m - 1]) && 
			(!next || part[m - 1] > next[m - 1])
		)
		{
			printf("2 x = %d, y = %d\n", m - 1, (int) (part - a->matrix) / m);
			sum += part[m - 1];
			count++;
		}

		for (int i = 1; i < m - 1; i++) {
			if (
				part[i] > part[i + 1] &&
				part[i] > part[i - 1] &&
				(!prev || part[i] > prev[prev_r * m + i]) &&
				(!next || part[i] > next[i])
			)
			{
				printf("3 x = %d, y = %d\n", i, (int) (part - a->matrix) / m);
				sum += part[i];
				count++;
			}
		}
	} else {
		if (
			part[0] > part[1] && 
			part[0] > part[m] && 
			(!prev || part[0] > prev[prev_r * m + 0])
		)
		{
			printf("4 x = %d, y = %d\n", 0, (int) (part - a->matrix) / m);
			sum += part[0];
			count++;
		}
		
		if (
			part[m - 1] > part[m - 2] && 
			part[m - 1] > part[2 * m - 1] && 
			(!prev || part[m - 1] > prev[prev_r * m + m - 1])
		)
		{
			printf("5 x = %d, y = %d\n", m - 1, (int) (part - a->matrix) / m);
			sum += part[m - 1];
			count++;
		}

		if (
			part[(r - 1) * m] > part[(r - 1) * m + 1] && 
			part[(r - 1) * m] > part[(r - 2) * m] && 
			(!next || part[(r - 1) * m] > next[0])
		)
		{
			printf("6 x = %d, y = %d\n", 0, (int) (part - a->matrix) / m + r - 1);
			sum += part[(r - 1) * m];
			count++;
		}

		if (
			part[r * m - 1] > part[r * m - 2] && 
			part[r * m - 1] > part[(r - 1) * m - 1] && 
			(!next || part[r * m - 1] > next[m - 1])
		)
		{
			printf("7 x = %d, y = %d\n", m - 1, (int) (part - a->matrix) / m + r - 1);
			sum += part[r * m - 1];
			count++;
		}

		for (int i = 1; i < m - 1; i++) {
			if (
				part[i] > part[i + 1] &&
				part[i] > part[i - 1] &&
				part[i] > part[m + i] &&
				(!prev || part[i] > prev[prev_r * m + i])
			)
			{
				printf("8 x = %d, y = %d\n", i, (int) (part - a->matrix) / m);
				sum += part[i];
				count++;
			}
		}

		for (int i = 1; i < m - 1; i++) {
			if (
				part[(r - 1) * m + i] > part[(r - 1) * m + i + 1] &&
				part[(r - 1) * m + i] > part[(r - 1) * m + i - 1] &&
				part[(r - 1) * m + i] > part[(r - 2) * m + i] &&
				(!next || part[(r - 1) * m + i] > next[i])
			)
			{
				printf("9 x = %d, y = %d\n", i, (int) (part - a->matrix) / m + r - 1);
				sum += part[(r - 1) * m + i];
				count++;
			}
		}

		for (int i = 1; i < r - 1; i++) {
			if (
				part[i * m] > part[i * m + 1] &&
				part[i * m] > part[(i - 1) * m] &&
				part[i * m] > part[(i + 1) * m]
			)
			{
				printf("10 x = %d, y = %d\n", 0, (int) (part - a->matrix) / m + i);
				sum += part[i * m];
				count++;
			}
		}

		for (int i = 2; i < r; i++) {
			if (
				part[i * m - 1] > part[i * m - 2] &&
				part[i * m - 1] > part[(i - 1) * m - 1] &&
				part[i * m - 1] > part[(i + 1) * m - 1]
			)
			{
				printf("11 x = %d, y = %d\n", m - 1, (int) (part - a->matrix) / m + i - 1);
				sum += part[i * m - 1];
				count++;
			}
		}

		for (int i = 1; i < r - 1; i++) {
			for (int j = 1; j < m - 1; j++) {
				if (
					part[i * m + j] > part[i * m + j + 1] &&
					part[i * m + j] > part[i * m + j - 1] &&
					part[i * m + j] > part[(i - 1) * m + j] &&
					part[i * m + j] > part[(i + 1) * m + j]
				)
				{
					printf("12 x = %d, y = %d\n", j, (int) (part - a->matrix) / m + i);
					sum += part[i * m + j];
					count++;
				}
			}
		}
	}
	
	a->sum = sum;
	a->count += count;
	a -> error_type = io_status::success;
	synchronize(a->p, & a->sum, 1);
	synchronize(a->p, & a->count, 1);


	double average = a->sum / a->count;
	if (r == 1) {
		if (
			part[0] > part[1] && 
			(!prev || part[0] > prev[0]) && 
			(!next || part[0] > next[0])
		)
		{
			part[0] = average;
		}
		if (
			part[m - 1] > part[m - 2] && 
			(!prev || part[m - 1] > prev[m - 1]) && 
			(!next || part[m - 1] > next[m - 1])
		)
		{
			part[m - 1] = average;
		}

		for (int i = 1; i < m - 1; i++) {
			if (
				part[i] > part[i + 1] &&
				part[i] > part[i - 1] &&
				(!prev || part[i] > prev[i]) &&
				(!next || part[i] > next[i])
			)
			{
				part[i] = average;
			}
		}
	} else {

		if (
			part[0] > part[1] && 
			part[0] > part[m] && 
			(!prev || part[0] > prev[0])
		)
		{
			part[0] = average;
		}
		
		if (
			part[m - 1] > part[m - 2] && 
			part[m - 1] > part[2 * m - 1] && 
			(!prev || part[m - 1] > prev[m - 1])
		)
		{
			part[m - 1] = average;
		}

		if (
			part[(r - 1) * m] > part[(r - 1) * m + 1] && 
			part[(r - 1) * m] > part[(r - 2) * m] && 
			(!next || part[(r - 1) * m] > next[0])
		)
		{
			part[(r - 1) * m] = average;
		}

		if (
			part[r * m - 1] > part[r * m - 2] && 
			part[r * m - 1] > part[(r - 1) * m - 1] && 
			(!next || part[r * m - 1] > next[m - 1])
		)
		{
			part[r * m - 1] = average;
		}

		for (int i = 1; i < m - 1; i++) {
			if (
				part[i] > part[i + 1] &&
				part[i] > part[i - 1] &&
				part[i] > part[m + i] &&
				(!prev || part[i] > prev[i])
			)
			{
				part[i] = average;
			}
		}

		for (int i = 1; i < m - 1; i++) {
			if (
				part[(r - 1) * m + i] > part[(r - 1) * m + i + 1] &&
				part[(r - 1) * m + i] > part[(r - 1) * m + i - 1] &&
				part[(r - 1) * m + i] > part[(r - 2) * m + i] &&
				(!next || part[(r - 1) * m + i] > next[i])
			)
			{
				part[(r - 1) * m + i] = average;
			}
		}

		for (int i = 1; i < r - 1; i++) {
			if (
				part[i * m] > part[i * m + 1] &&
				part[i * m] > part[(i - 1) * m] &&
				part[i * m] > part[(i + 1) * m]
			)
			{
				part[i * m] = average;
			}
		}

		for (int i = 2; i < r; i++) {
			if (
				part[i * m - 1] > part[i * m - 2] &&
				part[i * m - 1] > part[(i - 1) * m - 1] &&
				part[i * m - 1] > part[(i + 1) * m - 1]
			)
			{
				part[i * m - 1] = average;
			}
		}

		for (int i = 1; i < r - 1; i++) {
			for (int j = 1; j < m - 1; j++) {
				if (
					part[i * m + j] > part[i * m + j + 1] &&
					part[i * m + j] > part[i * m + j - 1] &&
					part[i * m + j] > part[(i - 1) * m + j] &&
					part[i * m + j] > part[(i + 1) * m + j]
				)
				{
					part[i * m + j] = average;
				}
			}
		}
	}

	return nullptr;
}


int process_args(Args *a)
{
	for (int i = 0; i < a->p; i++) {
		if (a[i].error_type != io_status::success) {
			printf("Error in file %s\n", a[i].name);
			return 1;
		}
	}
	return 0;
}


int read_file(char *filename, double *matrix, int n, int m)
{
	FILE *file = fopen(filename, "r");

	if (!file) {
		return -2;
	}

	for (int i = 0; i < n * m; i++) {
		if (fscanf(file, "%lf", matrix + i) != 1) {
			fclose(file);
			return -1;
		}
	}
	fclose(file);
	return 0;
}

void print_matrix(double *matrix, int n, int m)
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			printf(" %8.3e", matrix[i * m + j]);
		}
		printf("\n");
	}
	printf("\n");
}

int main(int argc, char *argv[])
{
	if (argc < 5) {
		printf("[-] Usage: %s files!\n", argv[0]);
		return 1;
	}

	int p, error, n, m;

	if (
		sscanf(argv[1], "%d", &p) != 1 ||
		sscanf(argv[2], "%d", &n) != 1 ||
		sscanf(argv[3], "%d", &m) != 1
	)
	{
		printf("[-] Bad argv!\n");
		return -1;
	}
	int pp = p;
	if (p > n) p = n;

    Args *a = new Args[p];
	double *matrix = new double[n * m]; 
	if (!a || !matrix) {
		if (a) delete[] a;
		if (matrix) delete[] matrix;
        printf("[-] Not enough memory!\n");
        return -2;
    }

	if (read_file(argv[4], matrix, n, m) != 0) {
		printf("[-] Bad file!\n");
		delete[] a;
		delete[] matrix;
		return -3; 
	}

	printf("---------------------- p = %2d ----------------------\n", pp);

	for (int k = 0; k < p; k++) {
		a[k].pi = k;
		a[k].p = p;
		a[k].n = n;
		a[k].m = m;
		a[k].r = n / p + (k < n % p ? 1 : 0);
		a[k].matrix = matrix;
		a[k].part = (k > 0 ? a[k - 1].part + a[k - 1].r * m : matrix);
		a[k].prev_r = (k > 0 ? a[k - 1].r : 0);

		if (k > 0) {
			a[k].prev = new double[m];
			if (!a[k].prev) {
				printf("[-] Not enough memory 2!\n");
				for (int j = 0; j <= k; j++) {
					if (a[j].prev) delete [] a[j].prev;
					if (a[j].next) delete [] a[j].next;
				}
				delete[] a;
				delete[] matrix;
				return -4;
			}
			for (int j = 0; j < m; j++) {
				a[k].prev[j] = matrix[(k - 1) * m + j];
			}
		}
		if (k < p - 1) {
			a[k].next = new double[m];
			if (!a[k].next) {
				printf("[-] Not enough memory 3!\n");
				for (int j = 0; j <= k; j++) {
					if (a[j].prev) delete [] a[j].prev;
					if (a[j].next) delete [] a[j].next;
				}
				delete[] a;
				delete[] matrix;
				return -4;
			}
			for (int j = 0; j < m; j++) {
				a[k].next[j] = a[k].part[a[k].r * m + j];
			}
		}
	}

	for (int k = 1; k < p; k++) {
		if ((error = pthread_create(&a[k].tid, nullptr, thread_func, a + k))) {
			for (int i = 0; i < p; i++) {
				if (a[i].prev) delete[] a[i].prev;
				if (a[i].next) delete[] a[i].next;
			}
			delete[] matrix;
			delete[] a;

			return error;
		}
	}

	a[0].tid = pthread_self();
	thread_func(a + 0);

	for (int k = 1; k < p; k++) {
		pthread_join(a[k].tid, nullptr);
	}
	// printf("Yep!\n");
	// if (!process_args(a)) {
	// 	printf("[+] Result = %lld\n", (long long) a[0].res);
	// }

	print_matrix(matrix, n, m);
	printf("Count = %d\n", (int) a->count);

	for (int i = 0; i < p; i++) {
		if (a[i].prev) delete[] a[i].prev;
		if (a[i].next) delete[] a[i].next;
	}
	delete[] matrix;
	delete[] a;
	return 0;
}