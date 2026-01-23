#ifndef _COLLISIONS_MPI_HH_
#define _COLLISIONS_MPI_HH_

#include "include.hh"

class collisions{
    protected:
        int myid, numprocs, num_integrators, N_bins, total_integrators;
        
        linspace_and_gl* eps;
        collision_integral** integrators;
        
        int load_value;
        
        int max_worker_bins;
        int** worker_values;
        int** worker_result_indexes;
        
    public:
        collisions(int, int, linspace_and_gl*, bool, bool, bool);
        ~collisions();
        
        int get_load_value();
        
        void compute_R(double, double, double*);
        void C(density*, density*, bool=true);
        
        void set_min_rate(density*);
        
        void print_coll();
};


#endif