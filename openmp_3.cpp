#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

double f(double x) {
    return 4 / (1 + x*x);
}

double integral_shared(size_t N, double start, double end, double dx) {
    double total = 0;
    double tmp_res = 0;

    #pragma omp parallel shared(total) firstprivate(tmp_res)
    {
        #pragma omp for
        for (size_t i = 0; i < N; ++i) {
            tmp_res += (f(start + i * dx) + f(start + (i + 1) * dx)) / 2;
        }
        #pragma omp atomic
        total += tmp_res;
    }
    total *= dx;
    return total;
}

double integral_reduction(size_t N, double start, double end, double dx) {
    double total = 0;

    #pragma omp parallel for reduction(+:total)
    for (size_t i = 1; i < N; ++i) {
        total += (f(start + i * dx) + f(start + (i + 1) * dx)) / 2;
    }
    total *= dx;
    return total;
}


int main(int argc, char** argv) {

    size_t N = 100'000'000;

    double start = 0;
    double end = 1;
    double dx = (end - start) / N;

    int p = 8;

    if (argc >= 2) {
        p = atoll(argv[1]);
    }

    omp_set_num_threads(p);
    printf("Using %d out of %d threads\n", p, omp_get_num_procs());

    auto fun = argc == 3 ? integral_reduction : integral_shared;

    double t_start = omp_get_wtime();
    double result = fun(N, start, end, dx);
    printf("Result: %lf (%lf s)\n", result, omp_get_wtime() - t_start);

    return 0;
}