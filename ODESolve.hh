#ifndef _ODESOLVE_HH_
#define _ODESOLVE_HH_
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include "CashKarp_vals.hh"
#include "arrays.hh"

using std::cout;
using std::endl;
using std::ofstream;
using std::abs;

using namespace std::chrono;

const double eps = 1e-8; 
const double TINY = 1e-22;
const double Safety = 0.9;


/*************************
/ template class ODESolve
/ Usage: create a class that inherits from ODESolve with a dep_vars inherited class (works only for a dep_vars inherited class, or dep_vars)
/ 
/ The constructor of the solver class must allocate the dep_vars object. The destructor should not. ODESolve destructor deletes y_values
/
/ protected variables (x, y, dx) is the current state of the system
/ set_ics sets these values
/ RKCK_step steps using Cash Karp, resulting in the next x, y, dx
/
/ f(double, dep_vars*, dep_vars*) must be defined for dy/dx = f(x, y)
/
*************************/

template <class dep>
class ODESolve
{
    protected:
        double x_value;
        dep* y_values;
        double dx_value;
        int total_ODE_steps;
        int total_ODE_rejected_steps;

        double eps = 1e-8; 
        double TINY = 1e-22;
        double Safety = 0.9;
        
    public:
        ODESolve();
        ~ODESolve();
        
        int get_rejected_steps();
        
        void set_tolerance(double);
        void set_tiny(double);
        void set_safety(double);
        
        void set_ics(double, dep*, double);
        void print_state();
        void print_csv(ostream&);
        
        virtual void f(double, dep*, dep*) = 0;
        void f_evaluate(dep*);
        
        void RKCash_Karp(double, dep*, double, double*, dep*, dep*);
        bool step_accept(dep*, dep*, dep*, double, double*, bool=false, bool=false);
        
        bool RKCK_step(double, dep*, double, double*, dep*, double*);
        bool RKCK_step_advance();
        bool ODEOneRun(int N_step, int dN, double x_final, const std::string& file_name, bool verbose = false, bool print_csv_file = true);
        
        bool run(int N_step, int dN, double x_final, const std::string& file_name, bool verbose = false);     
        
   
};


template <class dep>
ODESolve<dep>::ODESolve()
{
    x_value = 0.;
    dx_value = 1.;
    
    total_ODE_steps = 0;
    total_ODE_rejected_steps = 0;
}

template <class dep>
ODESolve<dep>::~ODESolve()
{   delete y_values; }

template <class dep>
void ODESolve<dep>::set_ics(double x0, dep* y0, double dx0)
{
    x_value = x0;
    y_values->copy(y0);
    dx_value = dx0;
}

template <class dep>
int ODESolve<dep>::get_rejected_steps()
{   return total_ODE_rejected_steps; }

template <class dep>
void ODESolve<dep>::set_tolerance(double tol)
{   eps = tol;}

template <class dep>
void ODESolve<dep>::set_tiny(double sm)
{   TINY = sm;}

template <class dep>
void ODESolve<dep>::set_safety(double s)
{   Safety = s;}

template <class dep>
void ODESolve<dep>::print_state()
{
    cout << "x = " << x_value << "; dx = " << dx_value << endl;
    y_values->print();
}

template <class dep>
void ODESolve<dep>::print_csv(ostream& os)
{
    os.precision(std::numeric_limits<double>::max_digits10 - 1);
    os << x_value << ", " << dx_value << ", ";
    y_values->print_csv(os);
    os << endl;
}

template <class dep>
void ODESolve<dep>::f_evaluate(dep* der)
{  f(x_value, y_values, der);  }

template <class dep>
void ODESolve<dep>::RKCash_Karp(double x, dep* y, double dx, double* x_stepped, dep* y_5th, dep* y_4th)
{
    int N = y->get_length();
    dep* k1 = new dep(y);
    dep* k2 = new dep(y);
    dep* k3 = new dep(y);
    dep* k4 = new dep(y);
    dep* k5 = new dep(y);
    dep* k6 = new dep(y);
    
    dep* z2 = new dep(y);
    dep* z3 = new dep(y);
    dep* z4 = new dep(y);
    dep* z5 = new dep(y);
    dep* z6 = new dep(y); 
    
    // k1 = dx * f(x, y)
    f(x, y, k1);
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
    return;
}

template <class dep>
bool ODESolve<dep>::step_accept(dep* y, dep* y5, dep* y4, double dx, double* dx_new, bool error_verbose, bool print_error_file)
{
    int N = y->get_length();

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
        
        return true;
    } 
    else if (dsm < 1){
        *dx_new = Safety * dx * pow(dsm, -0.2);
        *dx_new = std::min(5.0 * dx, *dx_new); 
        
        return true;
    }
    else{
        *dx_new = Safety * dx * pow(dsm, -0.25);
        *dx_new = std::min(0.5 * dx, *dx_new);
        
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
            dep* y_temp = new dep(y_values);
            dep* f_temp = new dep(y_values);
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
        return false;
    }
}

template <class dep>
bool ODESolve<dep>::RKCK_step(double x, dep* y, double dx, double* x_next, dep* y_next, double* dx_next)
{
    double dx_try = dx;
    double dx_future, x_future;
    int N = y->get_length();
    dep* y5 = new dep(y); 
    dep* y4 = new dep(y);
    bool accept = false;
    
    for (int i = 0; i<10; i++)        
    { 
        RKCash_Karp(x, y, dx_try, &x_future, y5, y4);
        if (step_accept(y, y5, y4, dx_try, &dx_future))
        {
            y_next -> copy(y5);
            *dx_next = dx_future;
            *x_next = x_future;
            accept = true;
            break;
        } 
        else {
            if (i < 10)
               dx_try = dx_future; 
        }
        
    }

    if (!accept)
    {
        cout << "ERROR:  10 iterations without acceptable step" << endl;
        cout << "x = " << x << "; dx = " << dx_try << endl;
        
        dx_try = dx;
        for (int i =0; i < 10; i++)
        {
            cout << "Step " << i << " ";
            RKCash_Karp(x, y, dx_try, &x_future, y5, y4);
            if (i < 9)  
                step_accept(y, y5, y4, dx_try, &dx_future, true, false);
            else
                step_accept(y, y5, y4, dx_try, &dx_future, true, true);
            dx_try = dx_future;
        }
    }
    
    if(y_next->isnan()){
        cout << "ERROR: nan present in dependent variables. Exit." << endl;
        return false;
    }
    
    delete y5;
    delete y4;
    
    return accept;
}

template <class dep>
bool ODESolve<dep>::RKCK_step_advance()
{   return RKCK_step(x_value, y_values, dx_value, &x_value, y_values, &dx_value);  }

template <class dep>
bool ODESolve<dep>::ODEOneRun(int N_step, int dN, double x_final, const std::string& file_name, bool verbose, bool print_csv_file)
{
    int N = y_values->get_length();
    
    // Initial values are (x_value, y_values, dx_value) as set by set_ics
    
    // Declare for RKCK_step
    
    bool no_error= true;
    bool done = false;
    
    if (x_final <= x_value)
    {
        cout << "x_final = " << x_final << " is less than initial condition, x_value = " << x_value << endl;
        return true;
    }
    
    ofstream file;
    if (print_csv_file)
        file.open(file_name);
    
    auto start = high_resolution_clock::now();
    
    if (verbose)
    {
        cout << "*******************" << endl;
        cout << "Running ODE Solver.  Initial Conditions:" << endl;
        print_state();
        cout << "Output printed to " << file_name << endl;
    }

    if (print_csv_file)
        print_csv(file);    
    
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
            
            if (x_value == x_final)
            {
                if(verbose)
                    cout << "Reached x_final" << endl;
                if (print_csv_file)
                    print_csv(file);
                done = true;
                break;
            }
        }
        if (!done && print_csv_file)
            print_csv(file);
    }
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    if (verbose)
    {
        print_state();
        cout << endl << "Time elapsed: "
         << duration.count()/1000. << " seconds" << endl;
        cout << "steps rejected / total steps = " << total_ODE_rejected_steps << " / " << total_ODE_steps << " (" << (100 * total_ODE_rejected_steps) / (total_ODE_rejected_steps+total_ODE_steps) << "%)" << endl;

    }

    if(print_csv_file)
        file.close();
    return done;
}

template <class dep>
bool ODESolve<dep>::run(int N_step, int dN, double x_final, const std::string& file_name, bool verbose)      
{
    return ODEOneRun(N_step, dN, x_final, file_name, verbose, true);
}

#endif