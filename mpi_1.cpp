
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<iostream>
using namespace std;

double func(double x){
    return 4/(1 + x*x);
}

double integral(const double x_0, const int p, const double dx) {
    double value = 0;
    double x = x_0;
    for (int i=0; i<p; ++i) {
        value += (func(x) + func(x + dx)) / 2 * dx;
        x += dx;
    }
    return value;
}

int main(int argc, char* argv[]) {
    int N = 100, a = 0, b = 1;
    double n = 100
    double dx = (b - a) / n;
    double integral_full = 0;
    double integral_i = 0;
    int myrank;
    int p;
    double begin, end, total;
    double I_0 = integral(a,N,dx)
    cout << "I_0: " << I_0 << endl;
    MPI_Status Status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    int size_p = N/p;
    int size_remain = N-size_p*p;
    
    MPI_Barrier(MPI_COMM_WORLD);
    begin = MPI_Wtime()

    integral_i = integral(a + myrank*size_p*dx, size_p, dx); 
    cout << "I_"<< myrank << ": " << integral_i << endl; 

    if (myrank != 0){
        MPI_Send(&integral_i, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    if (myrank == 0) {
        integral_full += integral_i;
        for (int i=1; i<p; ++i) {
            MPI_Recv(&integral_i, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
            integral_full += integral_i;
        }
        integral_full += integral(a + p * size_p * dx, size_remain, dx);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0){
        cout <<"I: " <<integrall_full << endl;
    }
    end = MPI_Wtime();
    total = end - begin;
    cout <<"Time " << p <<": "<< total << endl;

    MPI_Finalize();

return 0;

}