#include "collisionsQKE.hh"
#include "arrays.hh"
#include "density.hh"
#include "collisionsQKE_MPI.hh"
#include "constants.hh"
#include "matrices.hh"
#include "base_arrays.hh"
#include "thermodynamics.hh"

#include "mpi.h"


#include <iostream>
#include <chrono>
#include <fstream>

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
    
    const std::string& R_output_file = std::string(argv[1]);
    
    std::ofstream R_output;
    if(myid==0){
        R_output.open(R_output_file);
    }
    
    double emax[4] = {10,15,20,25};
    double nbins[3] = {100,200,300};
    
    
    for(int emax_index=0; emax_index<4; emax_index++){
        for(int nbins_index=0; nbins_index<3; nbins_index++){
            
            linspace_and_gl* eps = new linspace_and_gl(0., emax[emax_index], nbins[nbins_index], 5);
            int N_bins = eps->get_length();

            collisions* C_MPI = new collisions(myid, numprocs, eps);

            double* R_values = new double[N_bins * 4];

            auto start = high_resolution_clock::now();

            C_MPI->compute_R(32., 32., R_values); 

            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(stop - start);
            
            if(myid == 0){
                R_output << emax[emax_index] << ", " << nbins[nbins_index] << ", ";
                for(int i = 0; i < eps->get_length()-1; i++){
                    R_output << eps->get_value(i) << ", ";
                    
                }
                R_output << eps->get_value(eps->get_length()-1) << endl;
                for(int i=0; i<eps->get_length(); i++){
                    R_output << R_values[i] << ", ";
                }
                R_output << R_values[eps->get_length()-1] << endl;
                cout << "Emax=" << emax[emax_index] << ", nbins=" << nbins[nbins_index] << endl << "Time elapsed: " << duration.count() / 1000. << " seconds" << endl;
            }
            
            delete[] R_values;
            delete C_MPI;

            delete eps;
        }
    }
    
    
    if(myid==0){
        R_output.close();
    }
    
    MPI_Finalize();
    return 0;
}