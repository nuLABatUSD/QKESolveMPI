#include "include.hh"

double interpolate_log_fifth(double, double*, double*);
double fifth_order_fit(double, double*, double*);
double interpolate_log_linear(double, double, double, double, double);
double linear(double, double, double, double, double);

double extrapolate_exponential(double, double, double, double, double);
double extrapolate_linear(double, double, double, double, double);

void three_vector_for_QKE::v_vacuum(double delta_m_squared, double cos_2theta, double sin_2theta ){
    values[0] = delta_m_squared / 2. * sin_2theta;
    values[1] = 0.;
    values[2] = - delta_m_squared / 2. * cos_2theta;
}


void three_vector_for_QKE::v_thermal(dummy_vars* q, density* d){

    cout << "V_T" << endl;

    dep_vars* d0 = new dep_vars(q->get_len()); 
    dep_vars* d1 = new dep_vars(q->get_len()); 
    dep_vars* d2 = new dep_vars(q->get_len());
    three_vector* dummy1 = new three_vector();
    three_vector* dummy2 = new three_vector();

    for(int i=0; i<q->get_len(); i++){
        d->p0_p(i, true, dummy1);
        d->p0_p(i, false, dummy2);

        d0->set_value(i, pow(q->get_value(i),3) * (dummy1->get_value(0) + dummy2->get_value(0)));
        d1->set_value(i, pow(q->get_value(i),3) * (dummy1->get_value(1) + dummy2->get_value(1)));
        d2->set_value(i, pow(q->get_value(i),3) * (dummy1->get_value(2) + dummy2->get_value(2))); 
    }
    values[0] = q->integrate(d0);
    values[1] = q->integrate(d1);
    values[2] = q->integrate(d2);

    delete d0;
    delete d1;
    delete d2;
    delete dummy1;
    delete dummy2;


    for(int i=0; i<3; i++){
        values[i] *= - v_th_const * pow(d->get_Tcm(), 4);
    }

    double energy_dens = 0;
    double pressure = 0;
    energy_and_pressure(_electron_mass_, d->get_T(), &energy_dens, &pressure);
    values[2] += -2 * sqrt(2) * pow(_W_boson_,-2) * _GF_ * (energy_dens+pressure);

}

void three_vector_for_QKE::v_density(dummy_vars* q, density* d){
    dep_vars* d0 = new dep_vars(q->get_len()); 
    dep_vars* d1 = new dep_vars(q->get_len()); 
    dep_vars* d2 = new dep_vars(q->get_len()); 
    three_vector* dummy1 = new three_vector();
    three_vector* dummy2 = new three_vector();

    for (int i=0; i<q->get_len(); i++){
        d->p0_p(i, true, dummy1);
        d->p0_p(i, false, dummy2);
        d0->set_value(i,pow(q->get_value(i),2) * (dummy1->get_value(0) - dummy2->get_value(0)));
        d1->set_value(i,pow(q->get_value(i),2) * (dummy1->get_value(1) - dummy2->get_value(1)));
        d2->set_value(i,pow(q->get_value(i),2) * (dummy1->get_value(2) - dummy2->get_value(2)));
    }
    values[0] = q->integrate(d0);
    values[1] = q->integrate(d1);
    values[2] = q->integrate(d2);

    delete d0;
    delete d1;
    delete d2;
    delete dummy1;
    delete dummy2;


    for (int i=0; i<3; i++){
        values[i] *=  v_dens_const * pow(d->get_Tcm(), 3);
    }
}

/*************************************

class density

*************************************/

density::density(int num, dummy_vars* eps):dep_vars(8*num+2)
{
    N_bins = num;
    E = new dummy_vars(eps);
}

density::density(int num, dummy_vars* eps, double* dvals):dep_vars(8*num+2){
    N_bins = num;
    E = new dummy_vars(eps);
    for(int i=0; i<N_bins*8+2; i++){
        values[i] = dvals[i];
    }
}

// creates a density object with FD spectra, with eta_nu and eta_mu as neutrino degeneracy parameters (and opposite for anti-neutrino)
density::density(dummy_vars* eps, double eta_nu, double eta_mu):dep_vars(8*eps->get_len()+2)
{
    N_bins = eps->get_len();
    E = new dummy_vars(eps);

    double fnu = 0;
    double fmu = 0;
    double fnubar = 0;
    double fmubar = 0;

    double eps_temp = 0.;

    for (int i=0; i<N_bins; i++){
        eps_temp = eps->get_value(i);
        fnu = 1 / (exp(eps_temp - eta_nu)+1);
        fmu = 1 / (exp(eps_temp - eta_mu)+1);
        values[4*i] = fnu + fmu;
        values[4*i+3] =  (fnu - fmu)/(fnu+fmu+1.e-240);

        fnubar = 1 / (exp(eps_temp + eta_nu)+1);
        fmubar = 1 / (exp(eps_temp + eta_mu)+1);
        values[4*N_bins + 4*i] = fnubar + fmubar;
        values[4*N_bins + 4*i+3] = (fnubar - fmubar)/(fnu+fmu+1.e-240);
    }

}

// Creates exponential distributions (replaced following constructor)
density::density(dummy_vars* eps, int A, int B):dep_vars(8*eps->get_len()+2)
{
    N_bins = eps->get_len();
    E = new dummy_vars(eps);

    double fnu = 0;
    double fmu = 0;
    double fnubar = 0;
    double fmubar = 0;

    double eps_temp = 0.;

    for (int i=0; i<N_bins; i++){
        eps_temp = eps->get_value(i);
        fnu = (double)(A) * exp(-1 * eps_temp) / 10.;
        fmu = (double)(B) * exp(-2 * eps_temp) / 10.;
        values[4*i] = fnu + fmu;
        values[4*i+3] =  (fnu - fmu)/(fnu+fmu+1.e-240);

        fnubar = fnu;
        fmubar = fmu;
        values[4*N_bins + 4*i] = fnubar + fmubar;
        values[4*N_bins + 4*i+3] = (fnubar - fmubar)/(fnu+fmu+1.e-240);
    }

}


/******************************

Creates density objects with exponential distributions of the form f = A k^3 e^{- k eps}. The normalization, A, is chosen so that the number density of all
four species are equal, A is set by the input max_f = A k^3 (for the largest k)

*******************************/
density::density(dummy_vars* eps, double T, double k_e, double k_m, double k_ebar, double k_mbar, double max_f) : dep_vars(8*eps->get_len()+2)
{
    N_bins = eps->get_len();
    E = new dummy_vars(eps);
    
    this->set_T(T);
    
    double A;
    
    A = max(k_e, k_m);
    A = max(A, k_ebar);
    A = max(A, k_mbar);
    
    A = max_f / pow(A, 3);
    
    
    double fnu = 0;
    double fmu = 0;
    double fnubar = 0;
    double fmubar = 0;

    double eps_temp = 0.;

    for (int i=0; i<N_bins; i++){
        eps_temp = eps->get_value(i);
        
        fnu = A * pow(k_e, 3) * exp(- k_e * eps_temp);
        fmu = A * pow(k_m, 3) * exp(- k_m * eps_temp);
        fnubar = A * pow(k_ebar, 3) * exp(- k_ebar * eps_temp);
        fmubar = A * pow(k_mbar, 3) * exp(- k_mbar * eps_temp);
        
        this->set_value(4*i, fnu+fmu);
        
        set_value(4*i+3, (fnu - fmu)/(fnu+fmu+1.e-240));
        set_value(4*N_bins + 4*i, fnubar + fmubar);
        set_value(4*N_bins + 4*i+3, (fnubar - fmubar)/(fnu+fmu+1.e-240));
    }
}



density::density(density* copy_me):dep_vars(copy_me)
{
    N_bins = copy_me->num_bins();
    E = new dummy_vars(copy_me->get_E());
}

density::~density(){    
    delete E; 
}

dummy_vars* density::get_E(){
    return E;
}

double density::get_E_value(int i){
    return E->get_value(i);
}

double density::get_T(){
    return values[N-2];
}

double density::get_Tcm(){
    return values[N-1];
}

int density::num_bins(){
    return N_bins;
}

void density::set_T(double T){ 
    values[N-2] = T;
    values[N-1] = T;
}

void density::set_Tcm(double T){
    values[N-1] = T;
}

void density::set_T_Tcm(double T, double Tcm){
    set_T(T);
    set_Tcm(Tcm);
}

void density::set_value(int b, bool nu, int comp, double val){
    int index = b;
    if (!nu)
        index += N_bins;
    set_value(4 * index + comp, val);
}

double density::p0(int t, bool neutrino){
    if(neutrino==true){
        return values[4*t];
    }

    else{
        if(4*t+N_bins*4>8*N_bins-1){
            std::cout << "Warning: p0 exceeded the end of the density array, attempting to use index " << 4*t+N_bins*4 << ", t=" << t << std::endl;
        }
    return values[4*t+N_bins*4];
    }
}

void density::p_vector(int t, bool neutrino, three_vector* p){
    if(neutrino==true){
        for(int i=0; i<3; i++){
            p->set_value(i, values[4*t+i+1]);
        }
    }
    else{
        for(int i=0; i<3; i++){
            if(N_bins*4+4*t+1+i>8*N_bins-1){
                std::cout << "Warning: p_vector exceeded the end of the density array" << std::endl;
            }
            p->set_value(i, values[N_bins*4+4*t+1+i]);
        }
    }
}

void density::p0_p(int t, bool neutrino, three_vector* p){
    if(neutrino==true){
        for(int i=0; i<3; i++){
            p->set_value(i,values[4*t+i+1]);
        }
        p->multiply_by(values[4*t]);
    }

    else{
        for(int i=0; i<3; i++){
            p->set_value(i, values[N_bins*4+4*t+i+1]);
        }
        p->multiply_by(values[4*t+N_bins*4]);
    }
}

void density::number_density(double* output)
    {
    dep_vars* nu_e = new dep_vars(N_bins);
    dep_vars* nu_mu = new dep_vars(N_bins);
    dep_vars* nubar_e = new dep_vars(N_bins);
    dep_vars* nubar_mu = new dep_vars(N_bins);

    double P0, P0bar, Pz, Pzbar, eps;
    for(int i = 0; i < N_bins; i++)
    {
        P0 = values[4*i];
        P0bar = values[4*i+N_bins*4];
        Pz = values[4*i+3];
        Pzbar = values[N_bins*4+4*i+3];
        eps = E->get_value(i);
        nu_e->set_value(i, 0.5 * P0 * (1 + Pz) * eps * eps);
        nu_mu->set_value(i, 0.5 * P0 * (1 - Pz) * eps * eps);
        nubar_e->set_value(i, 0.5 * P0bar * (1 + Pzbar) * eps * eps);
        nubar_mu->set_value(i, 0.5 * P0bar * (1 - Pzbar) * eps * eps);
    }

    double norm = pow(values[N-1], 3) / (2 * _PI_ * _PI_);
    output[0] = E->integrate(nu_e) * norm;
    output[1] = E->integrate(nu_mu) * norm;
    output[2] = E->integrate(nubar_e) * norm;
    output[3] = E->integrate(nubar_mu) * norm;

    delete nu_e;
    delete nu_mu;
    delete nubar_e;
    delete nubar_mu;
}

void density::energy_density(double* output)
    {
    dep_vars* nu_e = new dep_vars(N_bins);
    dep_vars* nu_mu = new dep_vars(N_bins);
    dep_vars* nubar_e = new dep_vars(N_bins);
    dep_vars* nubar_mu = new dep_vars(N_bins);

    double P0, P0bar, Pz, Pzbar, eps;
    for(int i = 0; i < N_bins; i++)
    {
        P0 = values[4*i];
        P0bar = values[4*i+N_bins*4];
        Pz = values[4*i+3];
        Pzbar = values[N_bins*4+4*i+3];
        eps = E->get_value(i);
        nu_e->set_value(i, 0.5 * P0 * (1 + Pz) * pow(eps,3));
        nu_mu->set_value(i, 0.5 * P0 * (1 - Pz) * pow(eps,3));
        nubar_e->set_value(i, 0.5 * P0bar * (1 + Pzbar) * pow(eps,3));
        nubar_mu->set_value(i, 0.5 * P0bar * (1 - Pzbar) * pow(eps,3));
    }

    double norm = pow(values[N-1], 4) / (2 * _PI_ * _PI_);
    output[0] = E->integrate(nu_e) * norm;
    output[1] = E->integrate(nu_mu) * norm;
    output[2] = E->integrate(nubar_e) * norm;
    output[3] = E->integrate(nubar_mu) * norm;

    delete nu_e;
    delete nu_mu;
    delete nubar_e;
    delete nubar_mu;
}



double density::von_neumann_entropy(){
    dep_vars* integrand = new dep_vars(N_bins);
    
    double eig1;
    double eig2;
    three_vector* p0p = new three_vector();
    for(int i=0; i<N_bins; i++){
        
        this->p0_p(i, true, p0p);
        eig1 = 0.5 * this->p0(i, true) + 0.5 * p0p->get_value(2);
        this->p0_p(i, false, p0p);
        eig2 = 0.5 * this->p0(i, false) - 0.5 * p0p->get_value(2);
        
        integrand->set_value(i, pow(E->get_value(i),2)*(eig1*log(eig1) + eig2 * log(eig2)));
        
    }
    
    double result = -pow(values[N-1],3)/(2 * pow(_PI_,2)) * E->integrate(integrand);
    delete p0p;
    delete integrand;
    
    return result;
}

double density::thermodynamic_entropy(bool neutrino){
    dep_vars* integrand = new dep_vars(N_bins);
    
    double eig1;
    double eig2;
    three_vector* p0p = new three_vector();
    for(int i=0; i<N_bins; i++){
        
        this->p0_p(i, neutrino, p0p);
        eig1 = 0.5 * this->p0(i, neutrino) + 0.5 * p0p->get_value(2);
        
        integrand->set_value(i, pow(E->get_value(i),2)*(eig1*log(eig1) + (1-eig1)*log(1-eig1)));
        
    }
    
    double result = -pow(values[N-1],3)/(2 * pow(_PI_,2)) * E->integrate(integrand);
    delete p0p;
    delete integrand;
    
    return result;
    
}

double density::interpolated_matrix(bool neutrino, int index, double p4_energy, three_vector* p0p){
    double* results = new double[4]();
    double p0;

    //if in linspace, do fifth order interpolation
    if(p4_energy <= E->get_max_linspace()){
        int ind = std::max(0, index-2);
        ind = std::min(E->get_len()-1-4, ind);
        double* eps_vals = new double[5]();
        double* matrix_vals= new double[5]();

        for(int i=0; i<5; i++){
            eps_vals[i] = E->get_value(ind+i);
        }

        //1/2(p0+p0pz)
        for(int j=0; j<5; j++){
            this->p0_p(ind+j, neutrino, p0p);
            matrix_vals[j] = 0.5 * (this->p0(ind+j, neutrino) + p0p->get_value(2));
        }
        results[0] = interpolate_log_fifth(p4_energy, eps_vals, matrix_vals);

        //1/2p0px
        for(int j=0; j<5; j++){
            this->p0_p(ind+j, neutrino, p0p);
            matrix_vals[j] = 0.5 * p0p->get_value(0);
        }
        results[1] = interpolate_log_fifth(p4_energy, eps_vals, matrix_vals);

        //1/2p0py
        for(int j=0; j<5; j++){
            this->p0_p(ind+j, neutrino, p0p);
            matrix_vals[j] = 0.5 * p0p->get_value(1);
        }
        results[2] = interpolate_log_fifth(p4_energy, eps_vals, matrix_vals);

        //1/2(p0-p0pz)
        for(int j=0; j<5; j++){
            this->p0_p(ind+j, neutrino, p0p);
            matrix_vals[j] = 0.5 * (this->p0(ind+j, neutrino) - p0p->get_value(2));
        }
        results[3] = interpolate_log_fifth(p4_energy, eps_vals, matrix_vals);


        delete[] eps_vals;
        delete[] matrix_vals;
    }


    //if not in linspace do linear interpolation
    else{
        if(index==E->get_len()-1){
            index = index-1;
        }
        double energy_one = E->get_value(index);
        double energy_two = E->get_value(index+1);
        three_vector* secondp0p = new three_vector();

        this->p0_p(index, neutrino, p0p);
        this->p0_p(index+1, neutrino, secondp0p);

        //1/2(p0+p0pz)
        results[0] = interpolate_log_linear(p4_energy, energy_one, energy_two, 0.5*(this->p0(index, neutrino) + p0p->get_value(2)), 0.5*(this->p0(index+1, neutrino) + secondp0p->get_value(2)));

        //1/2(p0px)
        results[1] = interpolate_log_linear(p4_energy, energy_one, energy_two, 0.5*p0p->get_value(0), 0.5*secondp0p->get_value(0));

        //1/2(p0py)
        results[2] = interpolate_log_linear(p4_energy, energy_one, energy_two, 0.5*p0p->get_value(1), 0.5*secondp0p->get_value(1));

        //1/2(p0-p0pz)
        results[3] = interpolate_log_linear(p4_energy, energy_one, energy_two, 0.5*(this->p0(index, neutrino) - p0p->get_value(2)), 0.5*(this->p0(index+1, neutrino) - secondp0p->get_value(2)));

        delete secondp0p;
    }


    p0p->set_value(0, 2*results[1]);
    p0p->set_value(1, 2*results[2]);
    p0p->set_value(2, results[0]-results[3]);

    
    //p0
    p0 = results[0] + results[3];
    delete[] results;

    
    return p0;

}

double interpolate_log_fifth(double x, double* x_vals, double* y_vals){
    double y_temp;


    for(int i=1; i<5; i++){
        if(y_vals[0] * y_vals[i] <= 0){
        return fifth_order_fit(x, x_vals, y_vals);
        }
    }

    double* y_log = new double[5]();
    for(int j=0; j<5; j++){
        y_log[j] = log(std::abs(y_vals[j]));
    }
    y_temp = fifth_order_fit(x, x_vals, y_log);
    
    delete[] y_log;

    if(y_vals[0] > 0){
        return exp(y_temp);
    }
    else{
        return -exp(y_temp);
    }

}
double fifth_order_fit(double x, double* x_vals, double* y_vals){
    double fit = 0;
    double Lj;

    for(int j=0; j<5; j++){
        Lj = 1.0;
        for(int i=0; i<5; i++){
            if(i != j){
                Lj *= (x - x_vals[i]) / (x_vals[j] - x_vals[i]);
            }
        }
        fit += y_vals[j] * Lj;
    }
    return fit;
}


double interpolate_log_linear(double x, double x_val1, double x_val2, double y_val1, double y_val2){

    if(y_val1 * y_val2 <= 0){
        return extrapolate_linear(x, x_val1, x_val2, y_val1, y_val2);
    }
    else{
        double y_temp = linear(x, x_val1, x_val2, log(std::abs(y_val1)), log(std::abs(y_val2)));
        if(y_val1 > 0){
            return exp(y_temp);
        }
        else{
            return -exp(y_temp);
        }
    }
}
double linear(double x, double x1, double x2, double y1, double y2){
    if(x2-x1==0){std::cout << "warning: attempting to divide by 0**" << x << std::endl;}
        double slope = (y2-y1)/(x2-x1);
        return slope * (x-x2) + y2;
    }

double extrapolate_exponential(double x, double x1, double x2, double y1, double y2){
    //note: this assumes x1<x2, so we expect y1>y2 because this is an exponential decay model
    if(y1==y2){
        return y1;
    }

    else{
        //model is Ce^(-ax)
        if(y1/y2 < 1){
            return extrapolate_linear(x, x1, x2, y1, y2);
        }
        else{
            if(x1-x2==0){std::cout << "warning: attempting to divide by 0" << x << std::endl;}
            double a = -log(y1/y2) / (x1-x2);
            double C = y1 * exp(a * x1);
            return C * exp(-a * x);
        }
    }
}


double extrapolate_linear(double x, double x1, double x2, double y1, double y2){
    if(x2-x1==0){std::cout << "warning: attempting to divide by 0**" << x << std::endl;}
    double slope = (y2-y1)/(x2-x1);
    double Delta = slope * (x - x2);
    double result = 0;

    if (Delta > 0){
        if (y2 < 1){
            return y2 + (1 - y2) * std::tanh(Delta/(1-y2));
        }
        else{
            return 1.0;
        }
    }
    else{
        if (y2 > -1){
            return y2 + (y2 + 1) * std::tanh(Delta/(y2+1));
        }
        else{
            return -1.0;
        }
    }
}
