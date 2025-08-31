#include "collisionsQKE_MPI.hh"
#include "collisionsQKE.hh"
#include "arrays.hh"
#include "density.hh"
#include "mpi.h"

#include <iostream>

using std::cout;
using std::endl;

collisions::collisions(int rank, int num_ranks, linspace_and_gl* e){
    myid = rank;
    numprocs = num_ranks;
    
    eps = new linspace_and_gl(e);
    N_bins = eps->get_length();
    
    num_integrators = 0;
    load_value = 0;
    
    max_worker_bins = (2 * N_bins) / (numprocs - 1);
    if(2*N_bins % (numprocs-1) != 0)
        max_worker_bins++;
    
    worker_values = new int*[numprocs];
    worker_result_indexes = new int*[numprocs];
    
    int k;
    for(int j = 1; j < numprocs; j++){
        worker_values[j] = new int[2];
        worker_values[j][0] = 0;
        for(int i = j-1; i < N_bins * 2; i += numprocs-1)
            worker_values[j][0]++;
        worker_result_indexes[j] = new int[worker_values[j][0]];
        
        k=0;
        for(int i = j-1; i < N_bins * 2; i += numprocs-1){
            worker_result_indexes[j][k] = i;
            k++;
        }
    }    
    
    if(myid != 0){
/*        for(int i = myid-1; i < N_bins*2; i += numprocs-1)
            num_integrators++;
        integrators = new collision_integral*[num_integrators];
        int j = 0;
        for(int i = myid-1; i < N_bins*2; i += numprocs-1){
            integrators[j] = new nu_nu_collision(i%N_bins, eps, i < N_bins);
            load_value += integrators[j]->estimate_load();
            j++;
        }
*/
        num_integrators = worker_values[myid][0];
        integrators = new collision_integral*[num_integrators];
        for(int j = 0; j < worker_values[myid][0]; j++){
            integrators[j] = new nu_nu_collision(worker_result_indexes[myid][j] % N_bins, eps, worker_result_indexes[myid][j] < N_bins);
            load_value += integrators[j]->estimate_load();
        }
    }
    else{
    }
}

collisions::~collisions(){
    if(myid == 0){
        for(int i = 1; i < numprocs; i++){
            cout << i << ", " ;
            for(int j = 0; j < worker_values[i][0]; j++)
                cout << worker_result_indexes[i][j] << ", ";
            cout << endl;
        }
            
    }
    delete eps;

    for(int i = 1; i < numprocs; i++){
        delete[] worker_result_indexes[i];
        delete[] worker_values[i];
    }
    
    delete[] worker_result_indexes;
    delete[] worker_values;
    
    if(myid != 0){    
        for(int j = 0; j < num_integrators; j++)
            delete integrators[j];
        delete[] integrators;
    }
    
}

int collisions::get_load_value()
{   return load_value;}

/**************************

results is an array of length 4 * N_bins
0 - N_bins-1 : R for C0, neutrinos (true)
N_bins - 2*N_bins -1 : R for Cz, neutrinos (true)
2*N_bins - 3*N_bins-1 : R for C0, anti-neutrinos (false)
3*N_bins - 4*N_bins-1 : R for Cz, anti-neutrinos (false)

**************************/

void collisions::compute_R(double Tcm, double T, double* result){
    
    double* out_vals = new double[4 * N_bins]();
    double my_ans = 0.;
    int sender, tag;
    MPI_Status status;
    
    double dummy_int[2];
    
    if(myid == 0){
        for(int i = 0; i < N_bins * 2; i++){
            MPI_Recv(dummy_int, 2, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            sender = status.MPI_SOURCE;
            tag = status.MPI_TAG;
            
            out_vals[tag] = dummy_int[0];
            out_vals[tag + N_bins] = dummy_int[1];
        }
    }
    else{
        for(int j = 0; j < num_integrators; j++){
            tag = integrators[j]->get_bin();
            if(!integrators[j]->is_neutrino())
                tag += 2 * N_bins;
            integrators[j]->compute_R(Tcm, T, dummy_int);
            MPI_Send(dummy_int, 2, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
        }
    }
    
    MPI_Bcast(out_vals, 4 * N_bins, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    for(int i = 0; i < 4 * N_bins; i++)
        result[i] = out_vals[i];
    
    delete[] out_vals;
}

void collisions::C(density* dens, density* output, bool net){
       double* out_vals = new double[8 * N_bins]();
       double my_ans = 0.;
       int sender, tag;
       MPI_Status status;
       
       double* dummy_int = new double[4*max_worker_bins]();
       double dummy[4];
       
       if(myid == 0){
           for(int i = 1; i < numprocs; i++){
               MPI_Recv(dummy_int, 4*max_worker_bins, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
               sender = status.MPI_SOURCE;
               tag = status.MPI_TAG;
               
               for(int j = 0; j < worker_values[tag][0]; j++){
                   for(int k = 0; k < 4; k++)
                       out_vals[4 * worker_result_indexes[tag][j] + k] += dummy_int[4*j + k];
               }
           }
       }
       else{
           for(int j = 0; j < num_integrators; j++){
               tag = myid;
               integrators[j]->whole_integral(dens, dummy, net);
               for(int k = 0; k < 4; k++)
                   dummy_int[4*j + k] = dummy[k];
           }
           MPI_Send(dummy_int, 4*max_worker_bins, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
       }
/*       double dummy_int[4];
       
       if(myid == 0){
           for(int i = 0; i < N_bins * 2; i++){
               MPI_Recv(dummy_int, 4, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
               sender = status.MPI_SOURCE;
               tag = status.MPI_TAG;
               
               for(int j = 0; j < 4; j++)
                   out_vals[4 * tag + j] = dummy_int[j];
           }
       }
       else{
           for(int j = 0; j < num_integrators; j++){
               tag = integrators[j]->get_bin();
               if(!integrators[j]->is_neutrino())
                   tag += N_bins;
               integrators[j]->whole_integral(dens, dummy_int, net);
               MPI_Send(dummy_int, 4, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
           }
       }
*/       
       MPI_Bcast(out_vals, 8 * N_bins, MPI_DOUBLE, 0, MPI_COMM_WORLD);
       
       for(int i = 0; i < 8 * N_bins; i++)
           output->set_value(i, out_vals[i]);
       
       delete[] out_vals;
       delete[] dummy_int;
 
}