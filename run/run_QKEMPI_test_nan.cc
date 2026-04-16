#include "../density_test.hh"

#include "../code/include.hh"
#include "mpi.h"


void QKEMPI::set_density_object_for_test(){
    for(int i = 0; i < y_values->get_length(); i++)
        y_values->set_value(i, test_density[i]);
}

int main(int argc, char* argv[]){
    int myid, numprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    linspace_and_gl* eps = new linspace_and_gl(0., EPS_MAX_LINSPACE, EPS_LINSPACE_POINTS, 5);
    density* ics = new density(eps, IC_TCM, IC_NU_E, IC_NU_MU, IC_NUBAR_E, IC_NUBAR_MU, IC_MAX_DISTFUN);
    ics->set_T_Tcm(IC_TEMP, IC_TCM);
        
    QKEMPI* sim = new QKEMPI(myid, numprocs, PARAM_SIN_2THETA, PARAM_DELTA_M_SQUARED, 0., PARAM_DT_INIT, eps, ics);

    ofstream file;
    file.open(argv[1]);
    
    sim->set_density_object_for_test();
    sim->print_test_steps(file, PARAM_DN);

    file.close();
    
    delete eps;
    delete ics;
    delete sim;

    MPI_Finalize();
    return 0;
}
