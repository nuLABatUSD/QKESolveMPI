#include "include.hh"

#ifndef USE_V_THERMAL
#define USE_V_THERMAL 0
#endif

QKE::QKE(linspace_and_gl* e, double sin_2theta, double delta_m_squared, density* dens) : ODESolve()
{
    epsilon = new linspace_and_gl(e);
    sin_2theta = sin_2theta;
    cos_2theta = sqrt(1 - pow(sin_2theta, 2));
    delta_m_squared = delta_m_squared;

    y_values = new density(dens);

    dummy_v_vac = new three_vector_for_QKE;
    dummy_v_vac->v_vacuum(delta_m_squared, cos_2theta, sin_2theta);

    x_value = 0.;
    dx_value = dummy_v_vac->magnitude();
    dx_value /= dens->get_T() * 0.1;
    dx_value = 0.001 / dx_value;
}

QKE::~QKE()
{
    delete dummy_v_vac;
    delete epsilon;
}


void QKE::f(double t, density* d1, density* d2)
{
    d2->zeros();
    
    three_vector_for_QKE* dummy_v_dens = new three_vector_for_QKE;
    three_vector_for_QKE* dummy_v_therm = new three_vector_for_QKE;

    three_vector* V_nu = new three_vector;
    three_vector* V_nubar = new three_vector;
    three_vector* p = new three_vector;
    three_vector* vcrossp = new three_vector;
    
    dummy_v_dens->v_density(epsilon, d1);
    if(USE_V_THERMAL == 1)
        dummy_v_therm->v_thermal(epsilon, d1);

    double Tcm = d1->get_Tcm();
    
    double en = 0.;
    double* integral_vals = new double[4];
    for (int i=1; i< epsilon->get_len(); i++){
        en = epsilon->get_value(i) * Tcm;
        
        V_nu->copy(dummy_v_dens);
        V_nu->add_to(1./en, dummy_v_vac);
        V_nu->add_to(en, dummy_v_therm);

        d1->p_vector(i, true, p);

        vcrossp->set_cross_product(V_nu, p);
        d2->set_value(4*i+1, vcrossp->get_value(0));
        d2->set_value(4*i+2, vcrossp->get_value(1));
        d2->set_value(4*i+3, vcrossp->get_value(2));
        
        V_nubar->copy(dummy_v_dens);
        V_nubar->add_to(-1./en, dummy_v_vac);
        V_nubar->add_to(-en, dummy_v_therm);

        d1->p_vector(i, false, p);
        vcrossp->set_cross_product(V_nubar, p);
        d2->set_value(4*(epsilon->get_len())+4*i+1, vcrossp->get_value(0));
        d2->set_value(4*(epsilon->get_len())+4*i+2, vcrossp->get_value(1));
        d2->set_value(4*(epsilon->get_len())+4*i+3, vcrossp->get_value(2));

    }
    
        
    delete dummy_v_dens;
    delete dummy_v_therm;
    delete V_nu;
    delete V_nubar;
    delete vcrossp;
    delete p;
    delete[] integral_vals;
    
   
   
}