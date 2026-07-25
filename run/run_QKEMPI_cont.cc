#include "../restart.hh"

#include "../code/include.hh"
#include "mpi.h"


int main(int argc, char* argv[]){
    int myid, numprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    linspace_and_gl* eps = new linspace_and_gl(0., EPS_MAX_LINSPACE, EPS_LINSPACE_POINTS, 5);
    density* ics = new density(eps->length(), eps, density_restart);
        
    QKEMPI* sim = new QKEMPI(myid, numprocs, PARAM_SIN_2THETA, PARAM_DELTA_M_SQUARED, x_restart, dx_restart, eps, ics);

    sim->run(atoi(argv[2]), atoi(argv[3]), PARAM_T_FINAL, argv[1], true);
    
    delete eps;
    delete ics;
    delete sim;

    MPI_Finalize();
    return 0;
}
