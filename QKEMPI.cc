#include "QKESolve.hh"
#include "arrays.hh"
#include "density.hh"
#include "collisionsQKE_MPI.hh"
#include "collisionsQKE.hh"
#include "CashKarp_vals.hh"
#include "QKEMPI.hh"
#include "run_params.hh"
#include "mpi.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <string>

using std::cout;
using std::endl;
using std::ofstream;
using std::abs;

using namespace std::chrono;

#ifndef ODE_SOLVER_TOLERANCE
#define ODE_SOLVER_TOLERANCE 1.e-8
#endif

#ifndef ODE_SOLVER_TINY
#define ODE_SOLVER_TINY (ODE_SOLVER_TOLERANCE * 1.e-16)
#endif

#ifndef ODE_SOLVER_SAFETY
#define ODE_SOLVER_SAFETY 0.9
#endif

QKEMPI::QKEMPI(int rank, int numranks, double sin2theta, double dm2, double x0, double dx0, linspace_and_gl* e, density* ic) {
    myid = rank;
    numprocs = numranks;

    epsilon =  new linspace_and_gl(e);
    
    sin_2theta = sin2theta;
    cos_2theta = sqrt(1 - pow(sin2theta, 2));
    delta_m_squared = dm2;

    x_value = x0;
    dx_value = dx0;
    y_values = new density(ic);
    
    just_h = new QKE(e, sin_2theta, delta_m_squared, y_values);
    
    nu_nu_coll = true;
    nu_e_coll = false;
    nu_e_ann = false;
    
    coll_integrator = new collisions(myid, numprocs, epsilon);
    
    tol = ODE_SOLVER_TOLERANCE;
    TINY = ODE_SOLVER_TINY;
    Safety = ODE_SOLVER_SAFETY;
    
    total_ODE_steps = 0;
    total_ODE_rejected_steps = 0;

    
    H_cross_P = new density(ic);
    H_cross_P->zeros();

    k1 = new density(ic);
    k2 = new density(ic);
    k3 = new density(ic);
    k4 = new density(ic);
    k5 = new density(ic);
    k6 = new density(ic);
    
    z2 = new density(ic);
    z3 = new density(ic);
    z4 = new density(ic);
    z5 = new density(ic);
    z6 = new density(ic);   
}

QKEMPI::~QKEMPI(){
    delete epsilon;
    delete y_values;
    
    delete just_h;
    delete coll_integrator;
    
    delete H_cross_P;
    delete k1;
    delete k2;
    delete k3;
    delete k4;
    delete k5;
    delete k6;
   
    delete z2;
    delete z3;
    delete z4;
    delete z5;
    delete z6;   
}

void QKEMPI::print_state()
{
    if(myid == 0){
        cout << "x = " << x_value << "; dx = " << dx_value << endl;
        y_values->print();
    }
}

void QKEMPI::print_csv(ostream& os)
{
    if(myid == 0){
        os.precision(std::numeric_limits<double>::max_digits10 - 1);
        os << x_value << ", " << dx_value << ", ";
        y_values->print_csv(os);
        os << endl;
    }
}

void QKEMPI::f(double t, density* dens, density* der){
    der->zeros();
    
    coll_integrator->C(dens, der);
    just_h->f(t, dens, H_cross_P);
    
    der->add_to(1., H_cross_P);
}

void QKEMPI::f_evaluate(density* der){
    f(x_value, y_values, der);
}

double QKEMPI::first_derivative(double t, density* d1, density* der, double dx, double* x_next){
    double new_dx = 0.;
    if(myid == 0){
        double next_dx = 0.;
        density* y_next = new density(d1);
        
        just_h->RKCK_step(t, d1, dx, x_next, y_next, &next_dx);
        new_dx = *x_next-t;
        
        delete y_next;
    }
    der->zeros();
    
    coll_integrator->C(d1, der);
    just_h->f(t, d1, H_cross_P);
    
    der->add_to(1., H_cross_P);

    MPI_Bcast(&new_dx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return new_dx;
}


void QKEMPI::RKCash_Karp(double x, density* y, double dx, double* x_stepped, density* y_5th, density* y_4th)
{
    // k1 = dx * f(x, y)
    dx = first_derivative(x, y, k1, dx, x_stepped);
    k1 -> multiply_by(dx);  //k1 = dx * f(x,y)
  
    // k2 = dx * f(x + a2*dx, y + b21*k1)
    z2 -> copy(y);           //z2 = y
    z2 -> add_to(b21, k1);      //z2 = y + b21*k1
    f(x + a2*dx, z2, k2);          //k2 = f(x+a2*dx, z2)
    k2 -> multiply_by(dx);     //dx*f(..)

    //k2->print(8,1);
    // k3 = dx * f(x + a3*dx, y + b31*k1 + b32*k2)
    z3 -> copy(y);           //z3 = y
    z3 -> add_to(b31, k1); //z3 = y + b31*k1
    z3 -> add_to(b32, k2);
    f(x + a3*dx, z3, k3);         // k3 = f(x + a3*dx, z3)
    k3 -> multiply_by(dx);  // k3 = dx*f(x + a3*dx, z3)
 
    // k4 = dx * f(x + a4*dx, y + b41*k1 + b42*k2 +b43*k3)
    z4 -> copy(y);           //z4 = y
    z4 -> add_to(b41, k1);  //z4 = y + b41*k1
    z4 -> add_to(b42, k2); //z4 = y + b41*k1 + b42*k2
    z4 -> add_to(b43, k3); //z4 = y + b41*k1 + b42*k2 + b43*k3
    f(x + a4*dx, z4, k4);         // k4 = f(x + a4*dx, z4)
    k4 -> multiply_by(dx);
        
    // k5 = dx * f(x + a5*dx, y + b51*k1 + b52*k2 + b53*k3 + b54*k4)
    z5 -> copy(y);           //z5 = y
    z5 -> add_to(b51, k1);      //z5 = y + b51*k1
    z5 -> add_to(b52, k2);      //z5 = y + b51*k1 + b52*k2
    z5 -> add_to(b53, k3);      //z5 = y + b51*k1 + b52*k2 + b53*k3
    z5 -> add_to(b54, k4);      //z5 = y + b51*k1 + b52*k2 + b53*k3 + b54*k4    
    f(x + a5*dx, z5, k5);         // k5 = f(x + a5*dx, z5)
    k5 -> multiply_by(dx);
    
    // k6 = dx * f(x + a6*dx, y + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5)
    z6 -> copy(y);           //z6 = y
    z6 -> add_to(b61, k1);      //z6 = y + b61*k1
    z6 -> add_to(b62, k2);      //z6 = y + b61*k1 + b62*k2
    z6 -> add_to(b63, k3);      //z6 = y + b61*k1 + b62*k2 + b63*k3
    z6 -> add_to(b64, k4);      //z6 = y + b61*k1 + b62*k2 + b63*k3 + b64*k4
    z6 -> add_to(b65, k5);      //z6 = y + b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5 
    f(x + a6*dx, z6, k6);         // k6 = f(x + a6*dx, z6)
    k6 -> multiply_by(dx);
     
    //y_5th = y + c1*k1 + c2*k2 + c3*k3 + c4*k4 + c5*k5 + c6*k6
    y_5th -> copy(y); //y_5th = y
    y_5th -> add_to(c1, k1);
    y_5th -> add_to(c2, k2);
    y_5th -> add_to(c3, k3);
    y_5th -> add_to(c4, k4);
    y_5th -> add_to(c5, k5);
    y_5th -> add_to(c6, k6);


    // y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4 + cstar5*k5 + cstar6*k6
    y_4th -> copy(y); //y_4th = y           
    y_4th -> add_to(cstar1, k1); //y_4th = y + cstar1*k1
    y_4th -> add_to(cstar2, k2); //y_4th = y + cstar1*k1 + cstar2*k2
    y_4th -> add_to(cstar3, k3); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3
    y_4th -> add_to(cstar4, k4); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4
    y_4th -> add_to(cstar5, k5); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4 + cstar5*k5
    y_4th -> add_to(cstar6, k6); //y_4th = y + cstar1*k1 + cstar2*k2 + cstar3*k3 + cstar4*k4 + cstar5*k5 + cstar6*k6

    // x_stepped = x + dx
    *x_stepped = x + dx;
    
    return;
}

bool QKEMPI::step_accept(density* y, density* y5, density* y4, double dx, double* dx_new, bool error_verbose, bool print_error_file)
{
    //everyone needs to have access to this
    //do i just have one processor do the work and then broadcast the results?
    
    bool accept;
    
    double result[2];
    
    if (myid == 0){
        int N = y->length();

        int problem = 0;

        double dsm = 0;
        double delta1 = 0;
        double delta0 = 0;

        for (int i = 0; i<N; i++)
        { 
            delta1 = abs(y5 -> get_value(i) - y4 -> get_value(i));
            delta0 = eps*(abs(y -> get_value(i)) + abs(y5 -> get_value(i) - y -> get_value(i))) + TINY;

            if (delta1/delta0 > dsm)
            { 
                dsm = delta1/delta0;
                problem = i;

             }
         }

        if (dsm == 0)
        {
            *dx_new = 5 * dx;
            accept = true;
        } 
        else if (dsm < 1){
            *dx_new = Safety * dx * pow(dsm, -0.2);
            *dx_new = std::min(5.0 * dx, *dx_new); 
            accept = true;
        }
        else{
            *dx_new = Safety * dx * pow(dsm, -0.25);
            if (error_verbose)
                cout << "dsm = " << dsm << ", dx = " << dx << endl << "problem index = " << problem << "; y5 = " << y5->get_value(problem) << "; y4 = " << y4->get_value(problem) << endl;
            
            if (print_error_file)
            {
                ofstream error_file("ODESolve_ERROR_STEP_ACCEPT.csv");
                print_csv(error_file);
                error_file.close();
                
                ofstream deriv_error_file("ODESolve_ERROR_fproblem.csv");
                
                deriv_error_file.precision(std::numeric_limits<double>::max_digits10 - 1);
                double x_temp = x_value;
                density* y_temp = new density(y_values);
                density* f_temp = new density(y_values);
                int N_values = 101;
                double dx_temp = dx_value / (N_values-1);
                for (int i = 0; i<N_values; i++)
                {
                    f(x_temp, y_temp, f_temp);
                    y_temp->add_to(dx_temp, f_temp);
                    deriv_error_file << x_temp << ", " << f_temp->get_value(problem) << ", " << y_temp->get_value(problem) << ", " << f_temp->get_value(problem-1) << ", " << f_temp->get_value(problem+1) << endl;
                    x_temp += dx_temp;
                }
                deriv_error_file.close();
                
                delete y_temp;
                delete f_temp;
                
            }
            total_ODE_rejected_steps++;

            accept = false;
        }
        
        if(accept)
            result[0] = 1.;
        else
            result[0] = -1.;
            
        result[1] = *dx_new;
    }
    MPI_Bcast(result, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if(myid > 0)
        *dx_new = result[1];
    return (result[0] > 0);
}


bool QKEMPI::RKCK_step(double x, density* y, double dx, double* x_next, density* y_next, double* dx_next)
{
    //everyone does everything except print out the error message
    double dx_try = dx;
    density* y5 = new density(y); //???
    density* y4 = new density(y);
    
    bool accept = false;
    for (int i = 0; i<10; i++){
        RKCash_Karp(x, y, dx_try, x_next, y5, y4);
        
        if (step_accept(y, y5, y4, dx_try, dx_next)){
            y_next -> copy(y5);
            accept = true;
            break;
        } 
        else{
           dx_try = *dx_next; 
        }

    }
    if(myid == 0){
        if (!accept)
        {
            cout << "ERROR:  10 iterations without acceptable step" << endl;
            cout << "x = " << x << "; dx = " << dx_try << endl;
        }
        
        if(y_next->isnan())
            cout << "ERROR: nans present in density object" << endl;
        
    }
    
    if(y_next->isnan()){
        return false;
    }

    delete y5;
    delete y4;
    
    return accept;
}

bool QKEMPI::RKCK_step_advance()
{   return RKCK_step(x_value, y_values, dx_value, &x_value, y_values, &dx_value);}

bool QKEMPI::ODEOneRun(int N_step, int dN, double x_final, const std::string& file_name, bool verbose, bool print_csv_file)
{
    //everyone does everything except write the results to the file
    
    bool no_error= true;
    bool done = false;
       
    if (x_final <= x_value){
        if(myid == 0)
            cout << "x_final = " << x_final << " is less than initial condition, x_value = " << x_value << endl;
        return true;
    }    
    
    ofstream file;
    auto start = high_resolution_clock::now();

    if(myid == 0){    
        if (print_csv_file){
            std::string data_filename = file_name + ".csv";
            file.open(data_filename);
        }

        if (verbose)
        {
            cout << "*******************" << endl;
            cout << "Running ODE Solver.  Initial Conditions:" << endl;
            print_state();
            if(print_csv_file)
                cout << "Output printed to " << file_name << endl;
            else
                cout << "No file output." << endl;
        }
        
        if (print_csv_file)
            print_csv(file);    
    }    
    
    for (int i = 0; i < N_step && no_error && !done; i++)
    {
        for(int j = 0; j < dN; j++)
        {
            if(x_value + dx_value > x_final)
                dx_value = x_final - x_value;
            
            if (!RKCK_step_advance())
            {
                no_error = false;
                break;
            }
            total_ODE_steps++;
            
            if (dx_value == 0 || x_value / dx_value > 1.e9){
                cout << "ISSUE: step size has crashed. Abort." << endl;
                no_error = false;
                break;
            }
            
            if (x_value == x_final)
            {
                if(myid == 0){
                    if(verbose)
                        cout << "Reached x_final" << endl;
                    if (print_csv_file)
                        print_csv(file);
                }
                done = true;
                break;
            }
        }
        if (myid == 0 && !done && print_csv_file)
            print_csv(file);
    }
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    if(myid == 0){
        if (verbose)
        {
            print_state();
            cout << endl << "Time elapsed: "
             << duration.count()/1000. << " seconds" << endl;
            cout << "steps rejected / total steps = " << total_ODE_rejected_steps << " / " << total_ODE_steps << " (" << (100 * total_ODE_rejected_steps) / (total_ODE_rejected_steps+total_ODE_steps) << "%)" << endl;
            cout << "First derivative rejected " << just_h->get_rejected_steps() << " steps without collisions" << endl;
        }
        
        if(print_csv_file){
            file.close();
            
            std::string eps_filename = file_name + "-eps.csv";
            file.open(eps_filename);
            epsilon->print_csv(file);            
        }
    }
    return done;
}

bool QKEMPI::run(int N_step, int dN, double x_final, const std::string& file_name, bool verbose)
{
    //we want everyone to run this
    return ODEOneRun(N_step, dN, x_final, file_name, verbose, true);
}


