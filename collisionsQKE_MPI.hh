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
        
    public:
        collisions(int, int, linspace_and_gl*);
        ~collisions();
        
        void compute_R(double, double, double*);
        void C(density*, density*, bool);
};


#endif