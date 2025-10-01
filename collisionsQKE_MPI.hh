#ifndef _COLLISIONS_MPI_HH_
#define _COLLISIONS_MPI_HH_

#include "collisionsQKE.hh"
#include "arrays.hh"
#include "density.hh"

class collisions{
    protected:
        int myid, numprocs, num_integrators, N_bins;
        
        linspace_and_gl* eps;
        collision_integral** integrators;
        
        int load_value;
        
        int max_worker_bins;
        int** worker_values;
        int** worker_result_indexes;
        
    public:
        collisions(int, int, linspace_and_gl*);
        ~collisions();
        
        int get_load_value();
        
        void compute_R(double, double, double*);
        void C(density*, density*, bool=true);
        double number_dens_sum_rule(density*);
        double energy_dens_sum_rule(density*);
};


#endif