#include "include.hh"

double energy_f(double u, double m, double T){
    return 2 / (pow(_PI_,2)) * pow(u,2) * pow(T,3) * sqrt(pow(u,2)*pow(T,2) + pow(m,2)) / (exp(sqrt(pow(u,2)*pow(T,2) + pow(m,2))/T) + 1);
}

double pressure_f(double u, double m, double T){
    return 2 / (3 * pow(_PI_,2)) * pow(u,4) * pow(T,5) / sqrt(pow(u,2)*pow(T,2) + pow(m,2)) / (exp(sqrt(pow(u,2)*pow(T,2) + pow(m,2))/T) + 1);
}
                                                               
void energy_and_pressure(double m, double T, double* rho, double* P){
    dep_vars* rho_vals = new dep_vars(50);
    dep_vars* P_vals = new dep_vars(50);
    gl_dummy_vars* x = new gl_dummy_vars(50);

    for(int i=0; i<50; i++){
        rho_vals->set_value(i, energy_f(x->get_value(i),m,T));
        P_vals->set_value(i, pressure_f(x->get_value(i),m,T));
    }
    
    *rho = x->integrate(rho_vals);
    *P = x->integrate(P_vals);
    
    delete rho_vals;
    delete P_vals;
    delete x;

}