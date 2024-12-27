#include <cstddef>
#include <iostream>
#include <pthread.h>
#include <cfloat>
#include <chrono>
#include <limits>

double eps = std::numeric_limits<double>::epsilon() * 32;


class Args
{
public:
	pthread_t tid = 1;
	double error_flag = 0;
	const char *name = nullptr;

	double *array = nullptr;
	double *curr_array = nullptr;
	double *prev_array = nullptr;
	double *next_array = nullptr;
	double l1, l2, r1, r2; // -2, -1, len, len + 1 элементы
	double ll = 0; // левый  элемент левой  подпоследовательности
	double rr = 0; // правый элемент правой подпоследовательности
	double lr = 0; // правый элемент левой подпоследовательности
	double rl = 0; // левый элемент правой подпоследовательности
	int *count = nullptr;
	int p = 0;    // количество потоков
	int pi = 0;   // номер потока
	int n = 0;	  // длина всего массива
	int len = 0;  // длина куска массива
	int LBound = 0; // левая  граница куска
	int RBound = 0; // правая граница куска
	int lli = -3; // левая  граница левой  подпоследовательности
	int lri = -3; // правая граница левой  подпоследовательности
	int rli = -3; // левая  граница правой подпоследовательности
	int rri = -3; // правая граница правой подпоследовательности

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
	int n,
	int len,

	double *array,
	double *curr_array,
	int *count,
	int &lli,
	int &lri,
	int &rli,
	int &rri
)
{
	int li = -3, ri = -3, pos = curr_array - array;
	if (len == 1) {
		if (pos > 1) {
			if (std::abs(curr_array[-2] - curr_array[-1] - curr_array[0]) < eps) {
				li = -2;
				ri = 0;
			}
		}
		if (pos > 0 && n - (pos + len) > 0) {
			if (std::abs(curr_array[-1] - curr_array[0] - curr_array[1]) < eps) {
				if (li < -2) li = -1;
				ri = 1;
			}
		}
		if (n - (pos + len) > 1) {
			if (std::abs(curr_array[0] - curr_array[1] - curr_array[2]) < eps) {
				if (li < -2) li = 0;
				ri = 2;
			}
		}
		lli = li;
		rli = li;
		lri = ri;
		rri = ri;
		return;
	} else if (len == 2) {
		if (pos > 1) {
			if (std::abs(curr_array[-2] - curr_array[-1] - curr_array[0]) < eps) {
				li = -2;
				ri = 0;
			}
		}
		if (pos > 0) {
			if (std::abs(curr_array[-1] - curr_array[0] - curr_array[1]) < eps) {
				if (li < -2) li = -1;
				ri = 1;
			}
		}
		if (n - pos - len > 0) {
			if (std::abs(curr_array[0] - curr_array[1] - curr_array[2]) < eps) {
				if (li < -2) li = 0;
				ri = 2;
			}
		}
		if (n - pos - len > 1) {
			if (std::abs(curr_array[1] - curr_array[2] - curr_array[3]) < eps) {
				if (li < -2) li = 1;
				ri = 3;
			}
		}
		lli = li;
		rli = li;
		lri = ri;
		rri = ri;
		return;
	}


	lli = li;
	rli = li;
	lri = ri;
	rri = ri;

	if (pos > 1) {
		if (std::abs(curr_array[-2] - curr_array[-1] - curr_array[0]) < eps) {
			li = -2;
			ri = 0;
		}
	}
	if (pos > 0) {
		if (std::abs(curr_array[-1] - curr_array[0] - curr_array[1]) < eps) {
			if (li < -2) li = -1;
			ri = 1;
		}
	}

	double average;
	for (int i = 0; i < len - 2; i++) {
		if (std::abs(curr_array[i] - curr_array[i + 1] - curr_array[i + 2]) < eps) {
			if (ri == -3) {
				li = i;
				ri = i + 2;
			} else if (i - ri > 1) {
				if (li < 1) {
					lli = li;
					lri = ri;
				} else {
					average = (curr_array[li] + curr_array[ri]) / 2;
					for (int j = li; j <= ri; j++) curr_array[j] = average;
					*count += ri - li + 1;
				}
				li = i;
				ri = i + 2;
			} else {
				ri = i + 2;
			}
		}
	}

	if (n - pos - len > 0) {
		if (std::abs(curr_array[len - 2] - curr_array[len - 1] - curr_array[len]) < eps) {
			if (ri == -3) {
				li = len - 2;
				ri = len;
			} else if (len - 2 - ri > 1) {
				if (li < 1) {
					lli = li;
					lri = ri;
				} else {
					average = (curr_array[li] + curr_array[ri]) / 2;
					for (int j = li; j <= ri; j++) curr_array[j] = average;
					*count += ri - li + 1;
				}
				li = len - 2;
				ri = len;
			} else {
				ri = len;
			}
		}
	}
	if (n - pos - len > 1) {
		if (std::abs(curr_array[len - 1] - curr_array[len] - curr_array[len + 1]) < eps) {
			if (ri == -3) {
				li = len - 1;
				ri = len + 1;
			} else if (len - 1 - ri > 1) {
				if (li < 1) {
					lli = li;
					lri = ri;
				} else {
					average = (curr_array[li] + curr_array[ri]) / 2;
					for (int j = li; j <= ri; j++) curr_array[j] = average;
					*count += ri - li + 1;
				}
				li = len - 1;
				ri = len + 1;
			} else {
				ri = len + 1;
			}
		}
	}

	if (li < 1) {
		lli = li;
		lri = ri;
	} else if (ri <= len - 2) {
		average = (curr_array[li] + curr_array[ri]) / 2;
		for (int j = li; j <= ri; j++) curr_array[j] = average;
		*count += ri - li + 1;
	}

	if (lli < -2) {
		lli = li;
		lri = ri;
	}

	rli = li;
	rri = ri;

	return;
}

void connection(Args *a){
	double average;
	double *curr_array = a[0].curr_array;
	double *next_array = nullptr;
	double *array = a->array;
	int *count = a->count;
	int p = a->p;
	int lli, lri, rli, ri, li, npos; // li, ri - global positions!

	if (p == 1) {
		if (a[0].lli >= 0) {
			li = a[0].lli;
			ri = a[0].lri;
			printf("li = %d, ri = %d\n", li, ri);

			average = (curr_array[li] + curr_array[ri]) / 2;
			for (int j = li; j <= ri; j++) curr_array[j] = average;
			*count += ri - li + 1;
		} if (a[0].lli != a[0].rli && a[0].rri == a[0].n - 1) {
			li = a[0].rli;
			ri = a[0].rri;

			average = (curr_array[li] + curr_array[ri]) / 2;
			for (int j = li; j <= ri; j++) curr_array[j] = average;
			*count += ri - li + 1;
		}
		return;
	}
	else npos = 0;

	lli = a[0].lli;
	lri = a[0].lri;
	rli = a[0].rli;


	if (lli > -3 && lli != rli) {
		average = (curr_array[lli] + curr_array[lri]) / 2;
		for (int j = lli; j <= lri; j++) curr_array[j] = average;
		*count += lri - lli + 1;
	}

	li = a[0].rli;
	ri = a[0].rri;

	for (int i = 0; i < p - 1; i++) {
		printf("li = %d, ri = %d\n", li, ri);
		curr_array = a[i].curr_array;
		next_array = a[i].next_array;
		npos = next_array - array;
		if (ri < -2) {
			if (a[i + 1].lli > -3) {
				li = a[i + 1].lli + npos;
				ri = a[i + 1].lri + npos;
			}
		} else if (a[i + 1].lli < -2) {
			average = (array[li] + array[ri]) / 2;
			for (int j = li; j <= ri; j++) array[j] = average;
			*count += ri - li + 1;

			li = -3;
			ri = -3;
		} else if (ri + 1 < a[i + 1].lli + npos) {
			average = (array[li] + array[ri]) / 2;
			for (int j = li; j <= ri; j++) array[j] = average;
			*count += ri - li + 1;

			if (a[i + 1].lli != a[i + 1].rli && a[i + 1].lli == 0) {
				li = a[i + 1].lli + npos;
				ri = a[i + 1].lri + npos;
				average = (array[li] + array[ri]) / 2;
				for (int j = li; j <= ri; j++) array[j] = average;
				*count += ri - li + 1;
			}

			li = a[i + 1].rli + npos;
			ri = a[i + 1].rri + npos;
		} else {
			if (a[i + 1].rri != a[i + 1].lri) {
				ri = a[i + 1].lri + npos;
				average = (array[li] + array[ri]) / 2;
				for (int j = li; j <= ri; j++) array[j] = average;
				*count += ri - li + 1;

				li = a[i + 1].rli + npos;
				ri = a[i + 1].rri + npos;
			} else {
				ri = a[i + 1].rri + npos;
			}
		}
	}
	printf("li = %d, ri = %d\n", li, ri);

	if (li >= -2) {
		if (ri == a[p - 1].rri + npos) {
			average = (array[li] + array[ri]) / 2;
			for (int j = li; j <= ri; j++) array[j] = average;
			*count += ri - li + 1;
		} else {
			li = a[p - 1].rli + npos;
			ri = a[p - 1].rri + npos;

			average = (array[li] + array[ri]) / 2;
			for (int j = li; j <= ri; j++) array[j] = average;
			*count += ri - li + 1;
		}
	}

	return;
}

void *thread_func(void *args)
{
	auto start = std::chrono::steady_clock::now();
	Args *a = (Args*) args;

	array_editor(
		a->n,
		a->len,
		a->array,
		a->curr_array,
		a->count,
		a->lli,
		a->lri,
		a->rli,
		a->rri
	);
	synchronize(a->p);

	if (a->pi == 0) {
		connection(a);
	}
	synchronize(a->p);



	auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = end - start;
	printf("Thread %3d: Time taken: %8.2f seconds\n", a->pi, duration.count());


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

	int p = 0, n = 0, error, pi, count = 0, pp;

	if (argc < 4 || sscanf(argv[1], "%d", &p) != 1 || sscanf(argv[2], "%d", &n) != 1) {
		printf("[-] Not enough memory!\n");
		return 2;
	}
	if (n <= 0 || p <= 0) {
		printf("[-] Bad args!\n");
		return 3;
	}
	pp = p;
	if (p > n) p = n;

    Args *args = new Args[p];
	double *array = new double[n];

	if (!args || !array) {
		if (args) delete[] args;
		if (array) delete[] array;
        printf("[-] Not enough memory!\n");
        return -1;
    }

	if ((error = read_file(argv[3], array, n) < 0)) {
		delete[] args;
		delete[] array;
		return error;
	}

	for (pi = 0; pi < p; pi++) {
		args[pi].pi = pi;
		args[pi].p = p;
		args[pi].n = n;
		args[pi].array = array;
		args[pi].len = n / p + (pi < n % p ? 1 : 0);
		args[pi].curr_array = pi > 0 ? args[pi - 1].curr_array + args[pi - 1].len : array;
		if (pi > 0) args[pi].prev_array = args[pi - 1].curr_array;
		if (pi > 0) args[pi - 1].next_array = args[pi].curr_array;



		args[pi].count = &count;
	}

	for (pi = 1; pi < p; pi++) {
		if ((error = pthread_create(&args[pi].tid, nullptr, thread_func, args + pi))) {
			return error;
		}
	}

	args[0].tid = pthread_self();
	thread_func(args + 0);

	for (pi = 1; pi < p; pi++) {
		pthread_join(args[pi].tid, nullptr);
	}
	
	printf("[+] RESULT %3d:", pp);
	for (pi = 0; pi < n; pi++) {
		printf(" %8.2e", array[pi]);
	}
	printf("\n");

	delete[] args;
	delete[] array;
	return 0;
}
