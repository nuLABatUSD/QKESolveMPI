#ifndef _MATRICES_HH_
#define _MATRICES_HH_

#include <complex>
#include "arrays.hh"
#include "density.hh"

using std::complex;

class matrix
{
    protected:
    complex<double> A0;
    complex_three_vector* A;
    
    public:
    matrix();
    matrix(complex<double>, complex_three_vector*);
    matrix(bool);
    
    void convert_p_to_matrix(double, three_vector*);
    void convert_p_to_matrix(density*, bool, int);
    void convert_p_to_identity_minus_matrix(double, three_vector*);
    void convert_p_to_identity_minus_matrix(density*, bool, int);
    void convert_this_to_identity_minus_this();
    
    void matrix_add(matrix*, matrix*);
    void multiply_by(complex<double>);
    void matrix_multiply(matrix*, matrix*);
    complex_three_vector* get_A();
    complex<double> get_A0();
    void print_all();
    
    ~matrix();
    
};


#endif