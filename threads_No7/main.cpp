#include <cstddef>
#include <cstdio>
#include <iostream>
#include <pthread.h>
#include <sys/time.h>
#include <time.h>

const int chunk = 1e4;

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

	int prev_count = 0;
	int all = 0;

	unsigned long long int n = 0;
	unsigned long long int begin = 3;
	unsigned long long int res = 0;
	unsigned long long int *count = nullptr;
	unsigned long long int *sum = nullptr;
	unsigned long long int prev_sum = 0;
	double t1 = 0;
	io_status error_type = io_status::undef;
	double error_flag = 0;
};

void synchronize(int p, double *a, int n);
io_status nums_counter(FILE * f, double* count);
void *thread_func(void *args);
int process_args(Args *a);

void synchronize(int p, unsigned long long int *a=nullptr, int n = 0) {
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
	static int t_in = 0;//количество вошедших потоков
	static int t_out = 0;//количество вышедших потоков
	static unsigned long long int *r = nullptr;
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

double get_time() {
    struct timeval buf;
    gettimeofday(&buf, 0);
    return buf.tv_sec + buf.tv_usec / 1e6;
}

void *thread_func(void *args)
{
	Args *a = (Args*) args;
	unsigned long long int n = a->n;
	int p = a->p;
	int pi = a->pi;
	unsigned long long int begin = 3 + pi * chunk;
	bool check;
	unsigned long long int cur_sum;
	unsigned long long int cur_count;

	a->t1 = get_time();

	while (true) {
		cur_sum = 0;
		cur_count = 0;
		for (size_t i = begin; i < begin + chunk; i += 2) {
			check = true;
			for (size_t j = 3; j * j <= i && check; j += 2) {
				if (i % j == 0) {
					check = false;
				}
			}
			if (check) {
				cur_sum += i;
				cur_count++;
			}
		}
		synchronize(p, &cur_sum, 1);
		synchronize(p, &cur_count, 1);

		// printf("cur = %llu, n = %llu, count = %lld\n", cur_count, n, *a->count);
		if (*a->count + cur_count == n) {
			if (pi == 0) *a->sum += cur_sum;
			synchronize(p);
			a->t1 -= get_time();
			return nullptr;
		} else if (*a->count + cur_count > n) {
			// printf("Yep!\n");
			if (pi == 0) {
				cur_sum = *a->sum;
				cur_count = *a->count;
				for (size_t i = begin; i < begin + chunk * p + 2; i += 2) {
					check = true;
					for (size_t j = 3; j * j <= i + 2 && check; j += 2) {
						if (i % j == 0) {
							check = false;
						}
					}
					if (check) {
						cur_sum += i;
						cur_count++;
						// printf("curr count = %lld, n = %lld\n", cur_count, n);
						if (cur_count == n) {
							*a->sum = cur_sum;
							*a->count = cur_count;
							synchronize(p);
							a->t1 -= get_time();
							return nullptr;
						}
					}
				}
			}
			synchronize(p);
			a->t1 -= get_time();
			return nullptr;
		}

		synchronize(p);
		if (pi == 0) {
			*a->sum += cur_sum;
			*a->count += cur_count;
		}
		begin += p * chunk;
	}


	a -> error_type = io_status::success;
	synchronize(a->p, & a->res, 1);

	a->t1 -= get_time();

	return nullptr;
}




int main(int argc, char *argv[])
{
	if (argc < 3) {
		printf("[-] Usage: %s files!\n", argv[0]);
		return 1;
	}

	int p, error, n;
	double t = get_time();
	
	if (
		sscanf(argv[1], "%d", &p) != 1 ||
		sscanf(argv[2], "%d", &n) != 1
	)
	{
		printf("[-] Bad argv!\n");
		return -1;
	}
	if (p <= 0 || n <= 0) {
		printf("[-] Bad args!\n");
		return -2;
	}

	if (n == 1) {
		printf("[+] Result = %11d, p = %3d, n = %5u\n", 2, p, n);
		return 0;
	} else if (n == 2) {
		printf("[+] Result = %11d, p = %3d, n = %5u\n", 5, p, n);
		return 0;
	}
    Args *a = new Args[p];
	if (a == nullptr) {
        printf("[-] Not enough memory!\n");
        return -1;
    }

	int k;
	unsigned long long int count = 1, sum = 2, begin = 3;
	for (k = 0; k < p; k++) {
		a[k].pi = k;
		a[k].p = p;
		a[k].count = &count;
		a[k].sum = &sum;
		a[k].begin = begin;
		a[k].n = n;

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

	t = get_time() - t;

	for (int k = 0; k < p; k++) {
		printf("[T] Time in thread %d = %lf\n", k, -a[k].t1);
	}
	printf("[T] All time = %lf\n", t);
	printf("[+] Result = %11llu, p = %2d, n = %5u\n", sum, p, n);

	delete[] a;
	return 0;
}