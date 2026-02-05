#include "../code/include.hh"
#include "mpi.h"


int main(int argc, char* argv[])
{
    int myid, numprocs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        
    linspace_and_gl* eps = new linspace_and_gl(0., 20., 201, 5);
    int N_bins = eps->get_length();
    
    density* ics = new density(eps, IC_TCM, IC_NU_E, IC_NU_MU, IC_NUBAR_E, IC_NUBAR_MU, IC_MAX_DISTFUN);
    density* der = new density(ics);

    collisions* C_MPI = new collisions(myid, numprocs, eps, true, true, false);
        
    auto start = high_resolution_clock::now();

    C_MPI->C(ics, der);
        
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    cout << myid << ", " << C_MPI->get_load_value() << endl;
    
    if(myid == 0){
        cout << "Time elapsed: " << duration.count() / 1000. << " seconds" << endl;

/*    ofstream file;
    file.open("der_nue_test.csv");

    der->print_csv(file);
    file << endl;
    
    file.close();*/
    }
    
    delete C_MPI;
        
    delete ics;
    delete der;

    delete eps;
    
    MPI_Finalize();
    return 0;
}