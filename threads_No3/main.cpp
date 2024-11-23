#include <iostream>
#include <pthread.h>
#include <cfloat>

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
	int k = 0; // номер потока
	pthread_t tid = 1;
	const char * name = nullptr;
	double res = 0; // количество элементов, больших среднего
	double sum = 0;
	double average = 0;
	double first;
	double second;
	double last;
	double second_last;
	double count = 0;
	double count_in_file = 0;
	double start_elements_use = 0;
	double finish_elements_use = 0;
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

void process_res(Args *r) {
	double sum = 0, a = 0, b = 0;
	int count_prev = 0, count = 0, a_used = 1, b_used = 1, c_used = 1;

	for (int i = 0; i < r[0].p; i++) {
		if (count_prev == 0) { // Если в предыдущих файлах было 0 элементов
			if (std::abs(r[i].count_in_file - 1) < DBL_EPSILON) {
				b = r[i].last;
				b_used = 0; // b не использован в сумме
			} else if (r[i].count_in_file >= 1.9) {
				a = r[i].second_last;
				b = r[i].last;
				a_used = 0;
				b_used = 0;
			}
			count_prev = (int) std::min(2., r[i].count_in_file);
		} else if (count_prev == 1) { // Если в предыдущих файлах был только 1 элемент
			if (std::abs(r[i].count_in_file - 1) < DBL_EPSILON) { // Если в текущем файле 1 элемент
				a = b;
				b = r[i].last;
				a_used = b_used;
				b_used = 0;
				count_prev = 2;
			} else if (std::abs(r[i].count_in_file - 2) < DBL_EPSILON) { // Если в текущем 2 элемента
				if (std::abs(b - r[i].first - r[i].second) < DBL_EPSILON) {
					if (b_used == 0) {
						sum += b;
						count++;
					}
					if (std::abs(r[i].start_elements_use - 1) < DBL_EPSILON) {
						sum += r[i].first;
						count++;
					} else if (std::abs(r[i].start_elements_use - 1) < DBL_EPSILON) {
						sum += r[i].first + r[i].second;
						count += 2;
					}
					a_used = 1;
					b_used = 1;
				} else {
					a_used = 0;
					b_used = 0;
				}
				a = r[i].second_last;
				b = r[i].last;
				count_prev = 2; 
			} else if (std::abs(r[i].count_in_file - 3) < DBL_EPSILON) { // Если в текущем 3 элемента
				if (std::abs(b - r[i].first - r[i].second) < DBL_EPSILON) {
					if (b_used == 0) {
						sum += b;
						count++;
					}
					if (std::abs(r[i].start_elements_use - 1) < DBL_EPSILON) {
						sum += r[i].first;
						count++;
					} else if (std::abs(r[i].start_elements_use - 1) < DBL_EPSILON) {
						sum += r[i].first + r[i].second;
						count += 2;
					}
					a_used = 1;
				} else {
					a_used = 0;
				}
				a = r[i].second_last;
				b = r[i].last;
				a_used = std::max((r[i].finish_elements_use < DBL_EPSILON ? 0 : 1), a_used);
				b_used = r[i].finish_elements_use < 1.1 ? 0 : 1;
				count_prev = 2;
			} else if (r[i].count_in_file >= 3.9) { // Если в текущем больше 3х
				if (std::abs(b - r[i].first - r[i].second) < DBL_EPSILON) {
					if (b_used == 0) {
						sum += b;
						count++;
					}
					if (std::abs(r[i].start_elements_use - 1) < DBL_EPSILON) {
						sum += r[i].first;
						count++;
					} else if (std::abs(r[i].start_elements_use - 1) < DBL_EPSILON) {
						sum += r[i].first + r[i].second;
						count += 2;
					}
				}
				a = r[i].second_last;
				b = r[i].last;
				a_used = r[i].finish_elements_use < DBL_EPSILON ? 0 : 1;
				b_used = r[i].finish_elements_use < 1.1 ? 0 : 1;
				count_prev = 2;
			}
		} else { // если раньше было хотя бы 2 элемента
			if (std::abs(r[i].count_in_file - 1) < DBL_EPSILON) { // Если в текущем файле 1 элемент
				if (std::abs(a - b - r[i].first) < DBL_EPSILON) {
					if (a_used == 0) {
						sum += a;
						count++;
					}
					if (b_used == 0) {
						sum += b;
						count++;
					}
					sum += r[i].first;
					count++;
					a_used = 1;
					b_used = 1;
				} else {
					a_used = b_used;
					b_used = 0;
				}
				a = b;
				b = r[i].last;
			} else if (r[i].count_in_file > 1.9) { // если количество в файле больше 1
				if (std::abs(a - b - r[i].first) < DBL_EPSILON) {
					if (a_used == 0) {
						sum += a;
						count++;
					}
					if (b_used == 0) {
						sum += b;
						count++;
						b_used = 1;
					}
					if (r[i].start_elements_use < 1.1) {
						sum += r[i].first;
						count++;
						c_used = 1;
					}
				} else {
					c_used = 0;
				}
				if (std::abs(b - r[i].first - r[i].second) < DBL_EPSILON) {
					printf("%lf = %lf + %lf\n", b, r[i].first, r[i].second);
					if (b_used == 0) {
						printf("Yep!\n");
						sum += b;
						count++;
					}
					if (c_used == 0 && std::abs(r[i].start_elements_use - 1) < DBL_EPSILON) {
						sum += r[i].first;
						count++;
					} else if (std::abs(r[i].start_elements_use) < DBL_EPSILON) {
						if (c_used == 0){ // использован ли first
							sum += r[i].first + r[i].second;
							count += 2;
						} else {
							sum += r[i].second;
							count++;
						}
					}
				}

				a = r[i].second_last;
				b = r[i].last;

				// вычисление того, использованы ли a и b в сумме
				if (std::abs(r[i].count_in_file - 2) < DBL_EPSILON) {
					if (std::abs(b - r[i].first - r[i].second) < DBL_EPSILON) {
						a_used = 1;
						b_used = 1;
					} else if (std::abs(a - b - r[i].first) < DBL_EPSILON) {
						a_used = 1;
						b_used = 0;
					} else {
						a_used = 0;
						b_used = 0;
					}
				} else if (std::abs(r[i].count_in_file - 3) < DBL_EPSILON) {
					if (r[i].finish_elements_use > 0.9){
						a_used = 1;
						b_used = 1;
					} else {
						a_used = (std::abs(b - r[i].first - r[i].second) < DBL_EPSILON ? 1 : 0);
						b_used = 0;
					}
				} else { // количество элементов больше 3
					a_used = r[i].finish_elements_use < 0.1 ? 0 : 1;
					b_used = r[i].finish_elements_use < 1.1 ? 0 : 1;
				}
			}
		}
		printf("a is %lf, used: %d\n", a, a_used);
		printf("b is %lf, used: %d\n", b, b_used);
		printf("sum is %.1lf, count is %d\n", sum, count);
		printf("finish_elem_use is %lf\n", r[i].finish_elements_use);
	} // end for
	r[0].sum += sum;
	r[0].count += count;
} // end process_res


io_status sum_counter(
	FILE * f,
	double* count,
	double* sum,
	double * first,
	double * second,
	double * second_last,
	double * last,
	double * count_in_file,
	double * start_elements_use,
	double * finish_elements_use
)
{
    double a, b, c;
    *count = 0;
    *sum = 0;
	*count_in_file = 0;
	*start_elements_use = 0;
	*finish_elements_use = 0;
    
	if (fscanf(f, "%lf", &b) == 1) {
		first[0] = b;
		last[0] = b;
		count_in_file[0] = 1;
	}

	if (fscanf(f, "%lf", &c) == 1) {
		second_last[0] = b;
		last[0] = c;
		second[0] = c;
		count_in_file[0] = 2;
    }
    
	for (;;) {
		a = b;
		b = c;
        if (fscanf(f, "%lf", &c) != 1) {
            break;
        }
		second_last[0] = b;
		last[0] = c;
		count_in_file[0]++;
		if (std::abs(a - b - c) <= DBL_EPSILON) {
			finish_elements_use[0] = 2;
			sum[0] += a;
			count[0]++;
			if (std::abs(count_in_file[0] - 3) < DBL_EPSILON){
				start_elements_use[0] = 2;
			}
			if (std::abs(count_in_file[0] - 4) < DBL_EPSILON && start_elements_use[0] < 1.1){
				start_elements_use[0] = 1;
			}
		} else {
			if (std::abs(finish_elements_use[0] - 2) < DBL_EPSILON) {
				sum[0] += a + b;
				count[0] += 2;
			}
			finish_elements_use[0] = std::max(finish_elements_use[0] - 1, 0.);
		}
    }
	if (std::abs(finish_elements_use[0] - 2) < DBL_EPSILON) {
		sum[0] += a + b;
		count[0] += 2;
	}
    if (feof(f)) {
        return io_status::success;
    } else {
        return io_status::error_read;
    }
}

io_status res_counter(FILE * f, double* count, double average)
{
    double num;
    *count = 0;
    
    for (;;) {
        if (fscanf(f, "%lf", &num) != 1) {
            break;
        }
		if (num > average) {
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
	FILE * fp = fopen(a->name, "r");
	
	if (fp == nullptr) {
		a->error_type = io_status::error_open;
		a->error_flag = 1;
	}
	else {
		a->error_flag = 0;
	}
	synchronize(a->p, & a->error_flag, 1);
	if (a -> error_flag > 0) {
		if (fp != nullptr) fclose(fp);
		return nullptr;
	}

	a->error_type = sum_counter(
		fp, 
		&a->count,
		&a->sum, 
		&a->first, 
		&a->second, 
		&a->second_last, 
		&a->last, 
		&a->count_in_file,
		&a->start_elements_use,
		&a->finish_elements_use
	);
	if (a->error_type != io_status::success) {
		a->error_flag = 1;
	} else {
		a->error_flag = 0;
	}
	synchronize(a->p, &a->error_flag, 1);
	if (a -> error_flag > 0) {
		if (fp != nullptr) fclose(fp);
		return nullptr;
	}
	synchronize(a->p);

	// printf("1 of %s is %lf\n", a->name, a->first);
	// printf("2 of %s is %lf\n", a->name, a->second);
	// printf("-1 of %s is %lf\n", a->name, a->last);
	// printf("-2 of %s is %lf\n", a->name, a->second_last);
	// printf("count nums in %s is %lf\n", a->name, a->count_in_file);
	printf("sum of %s is %lf\n", a->name, a->sum);
	printf("count of %s is %lf\n", a->name, a->count);
	synchronize(a->p);

	if (a->k == 0) {
		process_res(a);
	}

	synchronize(a->p, & a->sum, 1);
	synchronize(a->p, & a->count, 1);


	if (a->count > 0.5) {
		a->average = a->sum / a->count;
	} else {
		printf("Yep!\n");
		a->res = 0;
		fclose(fp);
		a -> error_type = io_status::success;
		synchronize(a->p, & a->res, 1);
		return nullptr;
	}

	if (a->k == 0) {
		printf("sum is %lf\n", a->sum);
		printf("count is %lf\n", a->count);
	}
	synchronize(a->p);

	rewind(fp);
	a->error_type = res_counter(fp, &a->res, a->average);
	if (a->error_type != io_status::success) {
		a->error_flag = 1;
	} else {
		a->error_flag = 0;
	}
	synchronize(a->p, &a->error_flag, 1);
	if (a -> error_flag > 0) {
		if (fp != nullptr) fclose(fp);
		return nullptr;
	}

	fclose(fp);
	a -> error_type = io_status::success;
	synchronize(a->p, & a->res, 1);

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


int main(int argc, char *argv[])
{
	Args *a;
	if (argc == 1) {
		printf("[-] Usage: %s files!\n", argv[0]);
		return 1;
	}

	int p = argc - 1, error;
    a = new Args[p];
	if (a == nullptr) {
        printf("[-] Not enough memory!\n");
        return -1;
    }

	int k;
	for (k = 0; k < p; k++) {
		a[k].k = k;
		a[k].p = p;
		a[k].name = argv[k + 1];
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
	// printf("Yep!\n");
	if (!process_args(a)) {
		printf("[+] Result = %lld\n", (long long) a[0].res);
	}
	delete[] a;
	return 0;
}
