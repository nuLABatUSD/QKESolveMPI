#include "collisionsQKE_MPI.hh"
#include "collisionsQKE.hh"
#include "arrays.hh"
#include "density.hh"
#include "mpi.h"

#include <iostream>

using std::cout;
using std::endl;

collisions::collisions(int rank, int num_ranks, linspace_and_gl* e, bool nu_nu, bool nu_e, bool nu_e_ann){
    myid = rank;
    numprocs = num_ranks;
    
    eps = new linspace_and_gl(e);
    N_bins = eps->get_length();
    
    num_integrators = 0;
    load_value = 0;
    
    total_integrators = 2 * (nu_nu + nu_e + nu_e_ann) * N_bins;
    
/*    max_worker_bins = total_integrators / (numprocs - 1);
    if(total_integrators % (numprocs-1) != 0)
        max_worker_bins++;
*/    
    worker_values = new int*[numprocs];
    worker_result_indexes = new int*[numprocs];
    
    int** scat = new int*[numprocs];
    
    int k;
    for(int j = 1; j < numprocs; j++){
        worker_values[j] = new int[2];
        worker_values[j][0] = 0;
        for(int i = j-1; i < total_integrators; i += numprocs-1)
            worker_values[j][0]++;
        worker_result_indexes[j] = new int[worker_values[j][0]];
        scat[j] = new int[worker_values[j][0]];
        
        k=0;
        for(int i = j-1; i < total_integrators; i += numprocs-1){
            worker_result_indexes[j][k] = i % (2 * N_bins);
            
            if (nu_nu){
                if (i < 2 * N_bins)
                    scat[j][k] = 0;
                else if (i < 4 * N_bins){
                    if(nu_e)
                        scat[j][k] = 1;
                    else
                        scat[j][k] = 2;
                }
                else
                    scat[j][k] = 2;
            }
            else{
                if (i < 2 * N_bins){
                    if (nu_e)
                        scat[j][k] = 1;
                    else
                        scat[j][k] = 2;
                }
                else
                    scat[j][k] = 2;
            }
            
            
            k++;
        }
    }    
    max_worker_bins = worker_values[1][0];
        
    if(myid != 0){
        num_integrators = worker_values[myid][0];
        integrators = new collision_integral*[num_integrators];
        for(int j = 0; j < num_integrators; j++){
            switch(scat[myid][j]){
                case(0):
                    integrators[j] = new nu_nu_collision(worker_result_indexes[myid][j] % N_bins, eps, worker_result_indexes[myid][j] < N_bins);                
                    break;
                case(1):
                    integrators[j] = new nu_e_collision(worker_result_indexes[myid][j] % N_bins, eps, worker_result_indexes[myid][j] < N_bins, 32.);                
                    break;
                case(2):
                    break;
            }
            load_value += integrators[j]->estimate_load();
        }
    }
    else{
    }
}

collisions::~collisions(){
/*    if(myid == 0){
        for(int i = 1; i < numprocs; i++){
            cout << i << ", " ;
            for(int j = 0; j < worker_values[i][0]; j++)
                cout << worker_result_indexes[i][j] << ", ";
            cout << endl;
        }
            
    }*/
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

void collisions::print_coll(){
    cout << "myid = " << myid << ", num_integrators = " << num_integrators << endl;
    for(int i = 0; i < num_integrators; i++)
        cout << "\t" << worker_result_indexes[myid][i] << ", " << integrators[i]->get_eps_value() << endl;
}

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
        for(int i = 0; i < total_integrators; i++){
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
              
       double* dummy_int = new double[4];
       
       if(myid == 0){
           for(int i = 0; i < total_integrators; i++){
               MPI_Recv(dummy_int, 4, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
               sender = status.MPI_SOURCE;
               tag = status.MPI_TAG;
               
               for(int j = 0; j < 4; j++)
                   out_vals[4 * tag + j] += dummy_int[j];
           }
       }
       else{
           for(int j = 0; j < num_integrators; j++){
                tag = worker_result_indexes[myid][j];
               integrators[j]->whole_integral(dens, dummy_int, net);
               MPI_Send(dummy_int, 4, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
           }
       }

// Code to communicate once per worker, as opposed to after every calculation (see above)... slightly slower on Mac Mini
/*        double* worker_res = new double[4*max_worker_bins];
        if(myid==0){
            for(int i = 0; i < numprocs-1; i++){
                MPI_Recv(worker_res, 4*max_worker_bins, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                sender = status.MPI_SOURCE;
                tag = status.MPI_TAG;
                
                for(int j = 0; j < worker_values[tag][0]; j++){
                    for(int k = 0; k < 4; k++)
                        out_vals[4 * worker_result_indexes[tag][j]+k] += worker_res[4*j+k];
                }
            }
        
        }
        else{
            for(int j = 0; j < num_integrators; j++){
                integrators[j]->whole_integral(dens, dummy_int, net);
                for(int k = 0; k < 4; k++)
                    worker_res[4*j+k] = dummy_int[k];
            }
            MPI_Send(worker_res, 4*max_worker_bins, MPI_DOUBLE, 0, myid, MPI_COMM_WORLD);
        }

        delete[] worker_res;*/
       
       MPI_Bcast(out_vals, 8 * N_bins, MPI_DOUBLE, 0, MPI_COMM_WORLD);
       
       
       double P0val;
       three_vector* Pval = new three_vector();
       
       for(int i = 0; i < N_bins; i++){
           output->set_value(i, true, 0, out_vals[4 * i]);
           output->set_value(i, false, 0, out_vals[4 * (i+N_bins)]);
           
           P0val = dens->p0(i, true);
           dens->p_vector(i, true, Pval);
           for(int j = 0; j < 3; j++)
               output->set_value(i, true, j+1, (-out_vals[4*i] * Pval->get_value(j) + out_vals[4*i+j+1]) / P0val);
               
           P0val = dens->p0(i, false);
           dens->p_vector(i, false, Pval);
           for(int j = 0; j < 3; j++)
               output->set_value(i, false, j+1, (-out_vals[4*(i+N_bins)] * Pval->get_value(j) + out_vals[4*(i+N_bins)+j+1]) / P0val);
       }
       
       delete Pval;
       
       delete[] out_vals;
       delete[] dummy_int;
 
}