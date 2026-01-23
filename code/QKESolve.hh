#ifndef __QKESOLVE_HH__
#define __QKESOLVE_HH__

#include "ODESolve.hh"
#include "arrays.hh"
#include "density.hh"

class QKE : public ODESolve<density>
{

    protected:
        double delta_m_squared;
        double cos_2theta;
        double sin_2theta;

        linspace_and_gl* epsilon;
        three_vector_for_QKE* dummy_v_vac;


    public:
        QKE(linspace_and_gl* epsilon, double cos_2theta, double delta_m_squared, density* dens);
        ~QKE();
        void f(double, density*, density*);
};

#endif