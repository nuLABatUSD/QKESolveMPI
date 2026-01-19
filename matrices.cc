#include "matrices.hh"
#include "arrays.hh"
#include "density.hh"
#include <complex>
#include <iostream>


using std::cout;
using std::endl;
using std::complex;
using std::abs;

double interpolate(double, double, double, double, double);
double extrapolate_exponential(double, double, double, double, double);
double extrapolate_linear(double, double, double, double, double);

double interpolate(double, int, double*, double*);


matrix::matrix(){
    A0 = complex<double> (0,0);
    A = new complex_three_vector();
    
}

matrix::matrix(complex<double> c, complex_three_vector* C){
    A0 = c;
    A = new complex_three_vector();
    for(int i=0; i<3; i++){
       A->set_value(i,C->get_value(i));
    }
    
}

matrix::matrix(bool id){
    if (id==true){
        A0 = complex<double> (1,0);
        A = new complex_three_vector();
    }
    
}

matrix::matrix(matrix* c){
    A0 = c->get_A0();
    A = new complex_three_vector(c->get_A());
}

complex<double> matrix::get_A0(){
    return A0;
}

complex_three_vector* matrix::get_A(){
    return A;
}

void matrix::convert_p_to_matrix(double p0, three_vector* p0p){
    A0 = complex<double> (0.5 * p0,0);
    A->make_complex(p0p);
    A->multiply_by(complex<double>(0.5,0));
}

void matrix::convert_p_to_matrix(density* dens, bool neutrino, int i){
    A0 = complex<double> (0.5 * dens->p0(i, neutrino),0);
    three_vector* p0p = new three_vector();
    dens->p0_p(i, neutrino, p0p);
    A->make_complex(p0p);
    A->multiply_by(complex<double>(0.5,0));
    delete p0p;
}

void matrix::convert_p_to_identity_minus_matrix(double p0, three_vector* p0p){
    A0 = complex<double> (1 - 0.5 * p0,0);
    A->make_complex(p0p);
    A->multiply_by(complex<double>(-0.5,0));
}

void matrix::convert_p_to_identity_minus_matrix(density* dens, bool neutrino, int i){
    A0 = complex<double> (1 - 0.5 * dens->p0(i, neutrino),0);
    three_vector* p0p = new three_vector();
    dens->p0_p(i, neutrino, p0p);
    A->make_complex(p0p);
    A->multiply_by(complex<double>(-0.5,0));
    delete p0p;
}

void matrix::convert_this_to_identity_minus_this()
{
    A0 = complex<double> (1.0) - A0;
    A->multiply_by(complex<double>(-1.,0.));
}

void matrix::print_all(){
    cout << "A_0: " << A0 << endl;
    cout << "A: (" << A->get_value(0) << ", " << A->get_value(1) << ", " << A->get_value(2) << ")" << endl;
}

//add two matrices, modifies whatever matrix it is called on
void matrix::matrix_add(matrix* C1, matrix* C2){
    complex_three_vector* C1_A = C1->get_A();
    complex_three_vector* C2_A = C2->get_A();
    A->add(C1_A, C2_A);
    A0 = C1->get_A0()+C2->get_A0();
    
    
}

//multiply whole matrix by complex number
void matrix::multiply_by(complex<double> c){
    A0 *= c;
    A->multiply_by(c);
}

//multiply two matrices, modifies whatever matrix it is called on
void matrix::matrix_multiply(matrix* C1, matrix* C2){
    complex_three_vector* C1_A = new complex_three_vector(C1->get_A());
    complex_three_vector* C2_A = new complex_three_vector(C2->get_A());
    complex<double> C1_A0 = C1->get_A0();
    complex<double> C2_A0 = C2->get_A0();
    
    complex_three_vector* cp = new complex_three_vector();
    cp->set_cross_product(C1_A, C2_A);
    
    
    A0 = C1_A0*C2_A0 + C1_A->dot_with(C2_A);
    C1_A->multiply_by(C2_A0);
    C2_A->multiply_by(C1_A0);
    A->add(C1_A, C2_A);
    A->add(A,cp);
    
    
    delete cp;
    delete C1_A;
    delete C2_A;
    
}

matrix::~matrix(){
   delete A; 
}