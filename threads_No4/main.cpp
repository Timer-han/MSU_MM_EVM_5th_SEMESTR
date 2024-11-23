#include <iostream>
#include <pthread.h>
#include <cfloat>
#include <chrono>


class Args
{
public:
	pthread_t tid = 1;
	int p = 0; // количество потоков
	int k = 0; // номер потока
	int len = 0; // количество элементов внутри куска массива
	int *count = nullptr;
	double *array = nullptr;
	double prev_last;
	double next_first;
	double error_flag = 0;
	const char *name = nullptr;
};

void synchronize(int p, double *a, int n);
void *thread_func(void *args);

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


void array_editor(
	double* a,
	int len,
	double prev_last,
	double next_first,
	int *count,
	int k,
	int p
)
{
	double buf, prev;
	prev = a[0];
	if (len == 1 && k != 0 && k != p - 1) {
		buf = (prev_last + next_first) / 2;
		if (std::abs(buf - a[0]) > DBL_EPSILON) count[0]++;
		a[0] = buf;
		return;
	}
	else if (len == 1) return;

	if (k != 0) {
		buf = (prev_last + a[1]) / 2;
		if (std::abs(buf - a[0]) > DBL_EPSILON) count[0]++;
		a[0] = buf;
	}

	for (int i = 1; i < len - 1; i++) {
		buf = (prev + a[i + 1]) / 2;
		if (std::abs(buf - a[i]) > DBL_EPSILON) count[0]++;
		prev = a[i];
		a[i] = buf;
	}

	if (k != p - 1) {
		buf = (next_first + prev) / 2;
		if (std::abs(buf - a[len - 1]) > DBL_EPSILON) count[0]++;
		a[len - 1] = buf;
	}

	return;
}

void *thread_func(void *args)
{
	auto start = std::chrono::steady_clock::now();
	Args *a = (Args*) args;

	array_editor(
		a->array,
		a->len,
		a->prev_last,
		a->next_first,
		a->count,
		a->k,
		a->p
	);
	synchronize(a->p);

	auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = end - start;
	printf("Thread %3ld: Time taken: %8.2f seconds\n", a->tid, duration.count());


	return nullptr;
}


int read_file(const char *file, double *a, int n) {
	FILE *f = fopen(file, "r");
	if (!f) {
		printf("[-] Can't open the file!\n");
		return -1;
	}

	for (int i = 0; i < n; i++) {
		if (fscanf(f, "%lf", &a[i]) != 1) {
			printf("[-] Incorrect file!\n");
			fclose(f);
			return -2;
		}
	}

	fclose(f);
	return 0;
}


int main(int argc, char *argv[])
{
	if (argc == 1) {
		printf("[-] Usage: %s files!\n", argv[0]);
		return 1;
	}

	int p = std::stoi(argv[1]), n = std::stoi(argv[2]), error, k, count = 0, pp = p;
	if (p > n) p = n;
    Args *a = new Args[p];
	double *array = new double[n];

	if (!a || !array) {
		if (a) delete[] a;
		if (array) delete[] array;
        printf("[-] Not enough memory!\n");
        return -1;
    }

	if ((error = read_file(argv[3], array, n) < 0)) {
		delete[] a;
		delete[] array;
		return error;
	}

	for (k = 0; k < p; k++) {
		a[k].k = k;
		a[k].p = p;
		a[k].len = n / p + (k < n % p ? 1 : 0);
		a[k].array = k > 0 ? a[k - 1].array + a[k - 1].len : array;
		if (k > 0) a[k - 1].next_first = a[k].array[0];
		if (k > 0)     a[k].prev_last  = a[k - 1].array[a[k - 1].len - 1];
		// printf("Yep!\n");
		a[k].count = &count;
	}

	for (k = 1; k < p; k++) {
		if ((error = pthread_create(&a[k].tid, nullptr, thread_func, a + k))) {
			return error;
		}
	}

	a[0].tid = pthread_self();
	thread_func(a + 0);

	for (k = 1; k < p; k++) {
		pthread_join(a[k].tid, nullptr);
	}
	
	printf("[+] RESULT %3d:", pp);
	for (k = 0; k < n; k++) {
		printf(" %8.2e", array[k]);
	}
	printf("\n");

	delete[] a;
	delete[] array;
	return 0;
}
