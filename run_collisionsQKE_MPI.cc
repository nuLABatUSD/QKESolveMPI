#include "collisionsQKE.hh"
#include "arrays.hh"
#include "density.hh"
#include "collisionsQKE_MPI.hh"

#include "mpi.h"


#include <iostream>
#include <chrono>

using std::cout;
using std::endl;

using namespace std;

using namespace std::chrono;

int main(int argc, char* argv[])
{
    int myid, numprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        
    linspace_and_gl* eps = new linspace_and_gl(0., 20., 201, 5);
    int N_bins = eps->get_length();

    collisions* C_MPI = new collisions(myid, numprocs, eps);
    
    double* R_values = new double[N_bins * 4];

    auto start = high_resolution_clock::now();
    
    C_MPI->compute_R(32., 32., R_values); 
        
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
//    cout << myid << ", " << C_MPI->get_load_value() << endl;
    
    if(myid == 0){
        for(int i = 0; i < eps->get_length(); i++)
            cout << eps->get_value(i) << ", " << R_values[i] << endl;
        cout << "Time elapsed: " << duration.count() / 1000. << " seconds" << endl;
    }

    delete[] R_values;
    delete C_MPI;
        
    delete eps;
    
    MPI_Finalize();
    return 0;
}