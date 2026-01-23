#ifndef _QKEMPI_HH_
#define _QKEMPI_HH_

#include "include.hh"

class QKEMPI{
    protected:
        double x_value;
        double dx_value;
        density* y_values;
        linspace_and_gl* epsilon;
        
        double delta_m_squared;
        double cos_2theta;
        double sin_2theta;
        
        QKE* just_h;
        
        bool nu_nu_coll;
        bool nu_e_coll;
        bool nu_e_ann;
        collisions* coll_integrator;

        density* H_cross_P;
        density* k1;
        density* k2;
        density* k3;
        density* k4;
        density* k5;
        density* k6;

        density* z2; 
        density* z3;
        density* z4;
        density* z5;
        density* z6;
        
        int myid;
        int numprocs;
        
        double tol;
        double TINY;
        double Safety;
        int total_ODE_steps;
        int total_ODE_rejected_steps;
    
    public:
        // (int rank, int numranks, double sin2theta, double dm2, double x0, double dx0, linspace_and_gl* e, density* ic)
        QKEMPI(int, int, double, double, double, double, linspace_and_gl*, density*);
        ~QKEMPI();
        
        void print_state();
        void print_csv(ostream&);

        void f(double, density*, density*);
        void f_evaluate(density*);
        double first_derivative(double, density*, density*, double, double*);

        void RKCash_Karp(double, density*, double, double*, density*, density*);
        bool step_accept(density*, density*, density*, double, double*, bool=false, bool=false);
        
        bool RKCK_step(double, density*, double, double*, density*, double*);
        bool RKCK_step_advance();
        bool ODEOneRun(int, int, double, const string&, bool=false, bool=true);
        
        bool run(int, int, double, const string&, bool=true);
};


#endif