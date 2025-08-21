#ifndef _ARRAYS_HH_
#define _ARRAYS_HH_
#include <complex>
#include <iostream>

#include "base_arrays.hh"

class dummy_vars;
class dep_vars;
class complex_three_vector;

using std::complex;
using std::ostream;


class gl_dummy_vars : public dummy_vars
{
    public:
    gl_dummy_vars(int, double=0);
    gl_dummy_vars(gl_dummy_vars*);
};

class gel_dummy_vars : public dummy_vars
{
    public:
    gel_dummy_vars(int, double, double);
    gel_dummy_vars(gel_dummy_vars*);
};

class linspace_and_gl : public dummy_vars
{
    protected:
    int num_lin;
    int num_gl;

    public:
    linspace_and_gl(double, double, int, int);
    linspace_and_gl(linspace_and_gl*);
    int get_num_lin();
    int get_num_gl();
};

class linspace_for_trap : public linspace_and_gl
{
    public:
        linspace_for_trap(double, double, int);
        linspace_for_trap(linspace_for_trap*);
};


class three_vector : public dep_vars
{
    public:
    three_vector(int Nv=3);
    three_vector(double, double, double);
    three_vector(double*);
    three_vector(three_vector*);

    void add(three_vector*, three_vector*);
    double dot_with(three_vector*);
    double magnitude_squared();
    double magnitude();
    void set_cross_product(three_vector*, three_vector*);
    
    void make_real(complex_three_vector*);

};

class complex_three_vector{
    protected:
    complex<double>* values;
        
    public:
    complex_three_vector(int Nv=3);
    complex_three_vector(complex<double>, complex<double>, complex<double>);
    complex_three_vector(complex<double>);
    complex_three_vector(complex_three_vector*);
    
    void print_all();
    complex<double> get_value(int);
    void set_value(int, complex<double>);
    void make_complex(three_vector*);

    void multiply_by(complex<double>);
    void add(complex_three_vector*, complex_three_vector*);
    complex<double> dot_with(complex_three_vector*);
    complex<double> magnitude_squared();
    complex<double> magnitude();
    void set_cross_product(complex_three_vector*, complex_three_vector*);
    
    ~complex_three_vector();
    
};


#endif