#ifndef _DENSITY_HH_
#define _DENSITY_HH_

#include "arrays.hh"
#include "constants.hh"

#include <iostream>
#include <cmath>

using std::ostream;

class density;

class three_vector_for_QKE : public three_vector
{
    protected:
        const double v_th_const = 8*sqrt(2)*_GF_/(3*pow(_Z_boson_,2));
        const double v_dens_const = sqrt(2)*_GF_ / (2 * pow(_PI_,2));
    public:
    void v_vacuum(double, double, double);
    void v_density(dummy_vars*, density*);
    void v_thermal(dummy_vars*, density*);
};

class density : public dep_vars
{
    protected:
    int N_bins;
    dummy_vars* E;
    
    public:
    
    density(int, dummy_vars*);
    density(int, dummy_vars*, double*);
    density(dummy_vars*, double, double);
    density(dummy_vars*, int, int); // old
    density(dummy_vars*, double, double, double, double, double, double); //new and improved
    density(density*);
    ~density();
    
    double get_E_value(int);
    dummy_vars* get_E();
    double get_T();
    double get_Tcm();
    void set_Tcm(double);
    int num_bins();
    
    void set_T(double);
    void set_T_Tcm(double, double);
    using array::set_value;
    void set_value(int, bool, int, double);
    
    double p0(int, bool);
    void p_vector(int, bool, three_vector*);
    void p0_p(int, bool, three_vector*);

    void number_density(double*);
    void energy_density(double*);
    double von_neumann_entropy();
    double thermodynamic_entropy(bool);
    
    double interpolated_matrix(bool, int, double, three_vector*);
};


#endif