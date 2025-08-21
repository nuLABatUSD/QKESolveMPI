#ifndef _BASE_ARRAYS_HH_
#define _BASE_ARRAYS_HH_

#include <iostream>

using std::ostream;

/******************************

class array:
    the base class for all arrays
    is an abstract class (cannot be directly instantiated); serves as the base class to all array-type objects
    
    contains the basic functions:
    - get and set values
    - get the length (both as get_length() and get_len(), identical use)
    - print routines to cout
    - requires all derived classes to have a method void print_csv(ostream&); should not include a line break
    
class dummy_vars:
    base class for integration (and for bins)
    contains two arrays: abcissas and weights (for quadrature)
    function print_csv() prints two lines: first line for the abcissas (values) and second for the weights

class dep_vars:
    base class for dependent variables, primarily to ODE solving, and integration
    also a base class for arrays

******************************/
class array{
    protected:
        int N;
        double* values;
        
    public:
        array(int);
    
        double get_value(int);
        int get_length();
        int get_len();
        
        void set_value(int, double);
        
        void print_all();
        void print(int=3, int=1);
        virtual void print_csv(ostream&) = 0;
};

class dep_vars;

class dummy_vars : public array{
    protected:
        double* weights;
        
        double max_linspace;
        
    public:
        dummy_vars();
        dummy_vars(int);
        dummy_vars(dummy_vars*);
        ~dummy_vars();
        
        double get_weight(int);
        double get_max_linspace();
        
        void set_weight(int, double);
        
        void set_trap_weights();
        
        int bin_below(double);
        int index_below_for_interpolation(double);
        
        double integrate(dep_vars*);
        double partial_integrate_end(int, dep_vars*);
        double integrate_pow(dep_vars*, double);
        double partial_integrate_pow_end(int, dep_vars*, double);
        
        void print_csv(ostream&);
        
};

class dep_vars : public array{
        
    public:
        dep_vars(int);
        dep_vars(double*, int);
        dep_vars(dep_vars*);
        ~dep_vars();
                
        bool isnan();
                
        void zeros();
        void copy(dep_vars*);

        void multiply_by(double);
        void multiply_by(dep_vars*);
        void add_to(double, dep_vars*);
        
        void print_csv(ostream& os);
};



#endif
