

#ifndef _COLLISIONS_QKE_HH_
#define _COLLISIONS_QKE_HH_

#include "arrays.hh"
#include "density.hh"
#include "matrices.hh"

#include <iostream>

using std::cout;
using std::endl;

const double _COMPUTE_R_ERROR_ = -999.;
const double tolerance_min_rate = 100.;

#define EPS_CUT_1 0
#define EPS_CUT_2 5
#define EPS_CUT_3 1
#define EPS_TRANS_2 2
#define EPS_LIM_1 3
#define EPS_LIM_2 4

#define NUE_Q2 0
#define NUE_E2 1
#define NUE_Q3 2
#define NUE_E3 3

class collision_integral{
    protected:
        int bin, p1, N_bins;
        double eps_value;
        bool neutrino;
        
        linspace_and_gl* eps;
        
        dep_vars* outer_vals;
        dummy_vars* outer_dummy_vars;
        
        dep_vars** inner_vals;
        sub_dummy_vars** inner_dummy_vars;
        
        int num_F;
        double*** F_values;
        
        double min_rate;
        
        matrix** mat_p;
        matrix** mat_minus_p;
        matrix** F_dummy;
        
        virtual void populate_F(density*, bool) = 0; //density* dens, bool net
        virtual double interior_integral(int, int) = 0; //int p2, int which_term
        
        void get_inner_matrix(density*, double, sub_dummy_vars*, int, bool, matrix*, bool);
        void get_inner_matrix(density*, double, int, int, bool, matrix*, bool);
        

    public:
        collision_integral(int, linspace_and_gl*, bool); //(int bin, dummy_vars* eps);
        virtual ~collision_integral();
        
        virtual int estimate_load() = 0;
        
        int get_bin();
        int get_num_bins();
        double get_eps_value();
        
        bool is_neutrino();
        
        virtual void whole_integral(density*, double*, bool) = 0; //density* dens, double* results, bool net
        void C(density*, double*);

        
        virtual void compute_R(double, double, double*) = 0;
        void set_min_rate(density*);
        double get_min_rate();
        
        void show_idv(int, int); //(info, outer), info = 0 (values); 1 (need interp); 2 (interp index)
};

class electron_collision_integral : public collision_integral{
    protected:
        double Tcm, scaled_me, me_squared;
        
        dep_vars* outer_vals_2;
        dummy_vars* outer_dummy_vars_2;
        
        dep_vars** inner_vals_2;
        sub_dummy_vars** inner_dummy_vars_2;
        
        matrix* G_L;
        matrix* G_R;

        virtual double interior_integral_2(int, int) = 0;
        using collision_integral::get_inner_matrix;
        void get_inner_matrix(density*, double, int, int, bool, matrix*, bool, bool);

    public:
        electron_collision_integral(int, linspace_and_gl*, bool, double);
        virtual ~electron_collision_integral();
        
        double mom_to_eps(double);
        double eps_to_mom(double);
        
        double get_Tcm();
        void set_Tcm(double);
        
        
        using collision_integral::show_idv;
        void show_idv(int, int, int); //(R1/2, info, outer)
};

class nu_nu_collision : public collision_integral{
    protected:
        sub_dummy_vars** p4_values;
            
        const int load_factor = 4 + 4; // Fvvsc plus Fvvbarsc
        void populate_F(density*, bool); 
        double interior_integral(int, int); 
        double interior_integral_2(int, int);

        double J(double, double, double);
        double K(double, double, double);
        
        void get_p4_matrix(density*, int, int, bool, matrix*, bool);
        
        
        void Fvvsc_components_term_1(density*, int, int, double*, three_vector*);
        void Fvvsc_components_term_2(density*, int,int, double*, three_vector*);
        void Fvvsc_components(density*, int, int, double*, three_vector*, bool);
        void Fvvsc_for_p1(density*, bool);

        void Fvvbarsc_components_term_1(density*, int, int, double*, three_vector*);
        void Fvvbarsc_components_term_2(density*, int,int, double*, three_vector*);
        void Fvvbarsc_components(density*, int, int, double*, three_vector*, bool);
        void Fvvbarsc_for_p1(density*, bool);
    public:
        nu_nu_collision(int, linspace_and_gl*, bool);
        ~nu_nu_collision();
        
        int estimate_load();
        
        void whole_integral(density*, double*, bool); 
        
        void compute_R(double, double, double*);
        
};


class nu_e_collision : public electron_collision_integral{
    protected:
        double**** R_elec_values;
        double*** epslim_values;
        
        void populate_F(density*, bool); 
        double interior_integral(int, int); 
        double interior_integral_2(int, int);

        double M_11(int, double*);
        double M_12(int, double*);
        double M_21(int, double*);
        double M_22(int, double*);              
        
        void epslim_R1(double, double*);
        void epslim_R2(double, double*);
        
        void F_LL_F_RR(density*, bool, double, double, int, int, bool, double*, three_vector*);
        void F_LR_F_RL(density*, bool, double, double, int, int, bool, double*, three_vector*);
        
        void F_R1_for_p1(density*, bool);
        void F_R2_for_p1(density*, bool);
    public:
        nu_e_collision(int, linspace_and_gl*, bool, double);
        ~nu_e_collision();

        int estimate_load();
        
        void whole_integral(density*, double*, bool); 
        
        void compute_R(double, double, double*);
        
};

double J1(double, double, double);
double J2(double, double);
double J3(double, double, double);

double K1(double, double);
double K2(double, double, double);
double K3(double, double, double);

#endif