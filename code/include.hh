#include "../run_params.hh"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <chrono>

#include <cmath>
#include <complex>

#include "../../base_code/CashKarp_vals.hh"
#include "../../base_code/gel_vals.hh"
#include "../../base_code/gl_vals.hh"
#include "constants.hh"


#include "thermodynamics.hh"

#include "../../base_code/base_arrays.hh"
    class array;
    class dummy_vars;
    class gl_dummy_vars;
    class gel_dummy_vars;
    class dep_vars;
    class linspace_and_gl;
    class linspace_for_trap;
    class sub_dummy_vars;
    class three_vector;
    class complex_three_vector;
#include "density.hh"
    class three_vector_for_QKE;
    class density;
#include "arrays.hh"
#include "matrices.hh"
    class matrix;

#include "collisionsQKE.hh"
    class collision_integral;
    class electron_collision_integral;
    class nu_nu_collision;
    class nu_e_collision;
#include "collisionsQKE_MPI.hh"
    class collisions;

#include "../../base_code/ODESolve.hh"
#include "QKESolve.hh"
    class QKE;
#include "QKEMPI.hh"
    class QKEMPI;

using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;

using namespace std::chrono;

using std::atof;
using std::string;

using std::abs;
using std::max;
using std::complex;

