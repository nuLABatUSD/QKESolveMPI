#include "QKEMPI.hh"
#include "arrays.hh"
#include "density.hh"
#include "mpi.h"

#include <iostream>

using std::cout;
using std::endl;

int main(int argc, char* argv[]){
    int myid, numprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    linspace_and_gl* eps = new linspace_and_gl(0., 20., 201, 5);
    density* ics = new density(eps, 32., 0.9, 2*0.9, 0.9, 2*0.9, 0.9);
    
    QKEMPI* sim = new QKEMPI(myid, numprocs, 0., 0., 0., 1.e10, eps, ics);
    
    sim->run(10, 1, 5.e20, "thermal.csv", true);
    
    
    delete eps;
    delete ics;
    delete sim;

    MPI_Finalize();
    return 0;
}