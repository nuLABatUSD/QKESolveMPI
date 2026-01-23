#include "include.hh"

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

sub_dummy_vars::sub_dummy_vars(dummy_vars* copy_all) : dummy_vars(copy_all){
    orig_bins = new dummy_vars(copy_all);

    need_interpolation = new bool[N];
    interpolation_indices = new int[N];
    
    for(int i = 0; i < N; i++){
        need_interpolation[i] = false;
        interpolation_indices[i] = i;
    }
}

sub_dummy_vars::sub_dummy_vars(dummy_vars* dv, double A, double B, int N_GL) : dummy_vars(){
    bool err = false;
    
//    cout << "sub_dummy_vars, " << A << ", " << B << endl;
    
    if(A < 0)
        err = true;
    else{
        if(B != INNER_INTEGRAL_INFINITY){
            if(B < 0)
                err = true;
            else if(A > B)
                err = true;
        }
    }
    if(err){
        cout << "ERROR: sub_dummy_vars called with A = " << A << " and B = " << B << endl;
        N = -1;
        return;    
    }

    int count_min, count_max, bot_shift, top_shift, len;

    count_min = dv->index_below_for_interpolation(A) + 1;
    count_max = -1;
    bot_shift = 1;
    top_shift = 0;
    
    if (count_min < dv->get_length()){
//        cout << dv->get_value(count_min) << endl;
        if(abs(dv->get_value(count_min-1)-A) < SUBDV_INTERP_SMALL){
            count_min--;
            bot_shift = 0;
        }
        else if(dv->get_value(count_min) - A < SUBDV_INTERP_SMALL)
            bot_shift = 0;
    }
            
    if (B == INNER_INTEGRAL_INFINITY){
        if(count_min < dv->get_length() - N_GL)
            N = dv->get_length() - count_min + bot_shift;
        else
            N = 2;
        
    }
    else{
        count_max = dv->index_below_for_interpolation(B);
//            cout << B - dv->get_value(count_max) << endl;
        if(B - dv->get_value(count_max) > SUBDV_INTERP_SMALL){
            top_shift = 1;
        }
    
        N = count_max - count_min + 1 + bot_shift + top_shift;
    }
    

    values = new double[N]();
    weights = new double[N]();
    need_interpolation = new bool[N];
    interpolation_indices = new int[N];

    orig_bins = new dummy_vars(dv);
    
    if (B == INNER_INTEGRAL_INFINITY){
        if(count_min >= dv->get_length() - N_GL){
            gl_dummy_vars* gldv = new gl_dummy_vars(2, A);
            for(int i = 0; i < 2; i++){
                values[i] = gldv->get_value(i);
                weights[i] = gldv->get_weight(i);
                need_interpolation[i] = true;
                interpolation_indices[i] = dv->index_below_for_interpolation(values[i]);
            }
            delete gldv;
            return;
        }
        else{
            top_shift = 0;
            count_max = dv->get_length() - N_GL - 1;   
        } 
    }
    
//    cout << "sub_dummy_vars, " << count_min << ", " << bot_shift << ", " << top_shift << endl;
    
    if (count_min == dv->get_length()){
        values[0] = A;
        values[1] = B;
        set_trap_weights();
        
        need_interpolation[0] = true;
        interpolation_indices[0] = dv->get_length()-1;
        
        need_interpolation[1] = true;
        interpolation_indices[1] = dv->get_length()-1;    
        return;
    }
    
    if(bot_shift == 1){
        values[0] = A;
        need_interpolation[0] = true;
        interpolation_indices[0] = count_min-1;
    }
    
    for(int j = bot_shift; j <= count_max - count_min + bot_shift; j++){
        values[j] = dv->get_value(count_min + j - bot_shift);
        need_interpolation[j] = false;
        interpolation_indices[j] = count_min + j - bot_shift;
    }
    
    if(top_shift == 1){
        values[N-1] = B;
        need_interpolation[N-1] = true;
        interpolation_indices[N-1] = count_max;
    }
    set_trap_weights();
    
    
    if (B == INNER_INTEGRAL_INFINITY){
        for(int i = 1; i < 6; i++){
            values[N-i] = dv->get_value_from_end(i);
            weights[N-i] = dv->get_weight_from_end(i);
            need_interpolation[N-i] = false;
            interpolation_indices[N-i] = dv->get_length() - i;    
        }
    }
    
    if(abs(A-values[0]) > 1.e-12){
        cout << "sub_dummy_vars: " << A << ", " << B << ", " << values[0] << ", " << values[N-1] <<  endl;
        cout << "** " << dv->index_below_for_interpolation(A) << endl;
        }
}

sub_dummy_vars::sub_dummy_vars(dummy_vars* dv, int num) : dummy_vars() {
    orig_bins = new dummy_vars(dv);
    N = num;
    
    values = new double[N]();
    weights = new double[N]();
    need_interpolation = new bool[N];
    interpolation_indices = new int[N];
}

sub_dummy_vars::~sub_dummy_vars(){
    if(N!= -1){
        delete orig_bins;
        delete[] need_interpolation;
        delete[] interpolation_indices;
    }
}

bool sub_dummy_vars::get_need_interp(int i)
{   return need_interpolation[i];}

int sub_dummy_vars::get_interp_index(int i)
{   return interpolation_indices[i];}

void sub_dummy_vars::set_interp(){
    int bin_below;
    for(int i = 0; i < N; i++){
        bin_below = orig_bins->index_below_for_interpolation(values[i]);
        
        need_interpolation[i] = true;
        interpolation_indices[i] = bin_below;
        
        if(abs(values[i] - orig_bins->get_value(bin_below)) < SUBDV_INTERP_SMALL){
            need_interpolation[i] = false;
            interpolation_indices[i] = bin_below;
        }
        
        if(bin_below -1 >= 0)
            if(abs(values[i] - orig_bins->get_value(bin_below-1)) < SUBDV_INTERP_SMALL){
                need_interpolation[i] = false;
                interpolation_indices[i] = bin_below - 1;
            }
            
        if(bin_below + 1 < orig_bins->get_length())
            if(abs(values[i] - orig_bins->get_value(bin_below+1)) < SUBDV_INTERP_SMALL){
                need_interpolation[i] = false;
                interpolation_indices[i] = bin_below + 1;
            }
    }
}


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











