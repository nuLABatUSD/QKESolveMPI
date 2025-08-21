#include "arrays.hh"
#include "gl_vals.hh"
#include "gel_vals.hh"

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using std::cout;
using std::endl;

gl_dummy_vars::gl_dummy_vars(int num_gl, double start) : dummy_vars(num_gl)
{
    const double* val;
    const double* w;
        
    switch(num_gl){
        case 2:
            val = xvals_2;
            w = wvals_2;
            break;
        case 5:
            val = xvals_5;
            w = wvals_5;
            break;
        case 10:
            val = xvals_10;
            w = wvals_10;
            break;
        case 50:
            val = xvals_50;
            w = wvals_50;
            break;
        default:
            cout << "Error: This Gauss Legendre number is not supported" << endl;
            return;    
    }
    
    for(int i = 0; i < num_gl; i++){
        values[i] = val[i] + start;
        weights[i] = w[i] * exp(val[i]);
    }
}

gl_dummy_vars::gl_dummy_vars(gl_dummy_vars* copy_me) : dummy_vars(copy_me)
{   ;}

gel_dummy_vars::gel_dummy_vars(int num_gel, double start, double end) : dummy_vars(num_gel)
{
    double half_width = (end - start) / 2.;
    double slope_shift = half_width;
    double shift = (end + start) / 2.;
    
    const double* val;
    const double* w;
        
    switch(num_gel){
        case 2:
            val = gel_vals_2;
            w = gel_weights_2;
            break;
        case 5:
            val = gel_vals_5;
            w = gel_weights_5;
            break;
        case 10:
            val = gel_vals_10;
            w = gel_weights_10;
            break;
        case 50:
            val = gel_vals_50;
            w = gel_weights_50;
            break;
        case 100:
            val = gel_vals;
            w = gel_weights;
            break;
        default:
            cout << "Error: This Gauss Legendre number is not supported" << endl;
            return;
    }
    
    for(int i = 0; i < num_gel; i++){
        values[i] = slope_shift * val[i] + shift;
        weights[i] = half_width * w[i];
    }
    
}

gel_dummy_vars::gel_dummy_vars(gel_dummy_vars* copy_me) : dummy_vars(copy_me)
{   ;}

linspace_and_gl::linspace_and_gl(double xmin, double xmax, int numlin, int numgl) : dummy_vars(numlin+numgl)
{
    num_lin = numlin;
    num_gl = numgl;
    
    double dx_val = (xmax - xmin) / (num_lin -1);
    for (int i = 0; i<num_lin; i++){
        values[i] = xmin + dx_val * i;
        weights[i] = dx_val;
    }
    
    weights[0] = dx_val / 2;
    weights[num_lin-1] = dx_val / 2;

    if (numgl > 0){
        gl_dummy_vars* gl = new gl_dummy_vars(numgl, xmax);
        
        for (int i = 0; i < numgl; i++){
            values[num_lin+i] = gl->get_value(i);
            weights[num_lin+i] = gl->get_weight(i);
        }
        delete gl;
    }
    
    max_linspace = values[num_lin-1];
}

linspace_and_gl::linspace_and_gl(linspace_and_gl* copy_me) : dummy_vars(copy_me)
{
    num_lin = copy_me->get_num_lin();
    num_gl = copy_me->get_num_gl();
}

int linspace_and_gl::get_num_lin()
{   return num_lin;}

int linspace_and_gl::get_num_gl()
{   return num_gl;}


linspace_for_trap::linspace_for_trap(double xmin, double xmax, int num) : linspace_and_gl(xmin, xmax, num, 0)
{ ; }

linspace_for_trap::linspace_for_trap(linspace_for_trap* copy_me) : linspace_and_gl(copy_me)
{ ;}

three_vector::three_vector(int Nv):dep_vars(3)
{;}

three_vector::three_vector(double x, double y, double z):dep_vars(3)
{
    values[0] = x;
    values[1] = y;
    values[2] = z;
}

three_vector::three_vector(double* copy_me):dep_vars(copy_me, 3)
{;}

three_vector::three_vector(three_vector* copy_me):dep_vars(copy_me)
{;}

void three_vector::add(three_vector* A, three_vector*B){
    values[0] = A->get_value(0) + B->get_value(0);
    values[1] = A->get_value(1) + B->get_value(1);
    values[2] = A->get_value(2) + B->get_value(2);
}

double three_vector::dot_with(three_vector* B)
{
    double dot = 0;
    for(int i = 0; i < 3; i++)
        dot += values[i] * B->get_value(i);
    return dot;
}

double three_vector::magnitude_squared()
{
    return dot_with(this);
}

double three_vector::magnitude()
{
    double sum = 0;
    for(int i =0; i < 3; i++)
        sum += pow(this->get_value(i),2);
    return sqrt(sum);
}

void three_vector::set_cross_product(three_vector* A, three_vector* B)
{
    values[0] = A->get_value(1) * B->get_value(2) - A->get_value(2) * B->get_value(1);
    values[1] = A->get_value(2) * B->get_value(0) - A->get_value(0) * B->get_value(2);
    values[2] = A->get_value(0) * B->get_value(1) - A->get_value(1) * B->get_value(0);
}

void three_vector::make_real(complex_three_vector* C)
{
    values[0] = real(C->get_value(0));
    values[1] = real(C->get_value(1));
    values[2] = real(C->get_value(2));
}





complex_three_vector::complex_three_vector(int Nv){
    values = new complex<double>[3]();
    
}

complex_three_vector::complex_three_vector(complex<double> x, complex<double> y, complex<double> z){
    values = new complex<double>[3]();
    
    values[0] = x;
    values[1] = y;
    values[2] = z;
    
}

complex_three_vector::complex_three_vector(complex<double> c){
    values = new complex<double>[3]();
    for (int i = 0; i < 3; i++)
        values[i] = c;
    
    
}

complex_three_vector::complex_three_vector(complex_three_vector* c){
    values = new complex<double>[3]();
    for (int i = 0; i < 3; i++)
        values[i] = c->get_value(i);
    
}

void complex_three_vector::print_all(){
    for (int i=0; i<3; i++){
        cout << values[i] << endl;
    }
    
}

complex<double> complex_three_vector::get_value(int i){
    return values[i];
}

void complex_three_vector::set_value(int i, complex<double> d){
    values[i] = d;    
}

void complex_three_vector::add(complex_three_vector* A, complex_three_vector* B){
    values[0] = A->get_value(0) + B->get_value(0);
    values[1] = A->get_value(1) + B->get_value(1);
    values[2] = A->get_value(2) + B->get_value(2);
}

void complex_three_vector::multiply_by(complex<double> a){
   for (int i=0; i<3; i++){
       values[i] *= a;
   }
}

complex<double> complex_three_vector::dot_with(complex_three_vector* B)
{
    complex<double> dot = 0;
    for(int i = 0; i < 3; i++)
        dot += values[i] * B->get_value(i);
    return dot;
}

complex<double> complex_three_vector::magnitude_squared()
{
    return dot_with(this);
}

complex<double> complex_three_vector::magnitude()
{
    complex<double> sum = 0;
    for(int i =0; i < 3; i++)
        sum += pow(this->get_value(i),2);
    return sqrt(sum);
}

void complex_three_vector::set_cross_product(complex_three_vector* A, complex_three_vector* B)
{
    values[0] = A->get_value(1) * B->get_value(2) - A->get_value(2) * B->get_value(1);
    values[1] = A->get_value(2) * B->get_value(0) - A->get_value(0) * B->get_value(2);
    values[2] = A->get_value(0) * B->get_value(1) - A->get_value(1) * B->get_value(0);
}

void complex_three_vector::make_complex(three_vector* A){
    values[0] = complex<double> (A->get_value(0),0);
    values[1] = complex<double> (A->get_value(1),0);
    values[2] = complex<double> (A->get_value(2),0);
}

complex_three_vector::~complex_three_vector()
{   delete[] values; }











