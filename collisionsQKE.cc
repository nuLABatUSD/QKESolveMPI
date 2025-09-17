#include "arrays.hh"
#include "base_arrays.hh"
#include "density.hh"
#include "collisionsQKE.hh"
#include "constants.hh"
#include "matrices.hh"
#include "run_params.hh"

#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

using std::abs;
using std::max;

collision_integral::collision_integral(int b, linspace_and_gl* e, bool nu){
    bin = b;
    p1 = bin;
    N_bins = e->get_length();
    eps_value = e->get_value(bin);
    neutrino = nu;
    
    eps = new linspace_and_gl(e);
    
    min_rate = 0.;
    
/**********************************

    outer_vals and outer_dummy_vars; inner_vals and inner_dummy_vars; F_values matrix are allocated in the individual collision integrals
    BUT they are deleted in the base class, because they all need to be deleted. Only arrays particular to integrals should be deleted in the individual destructors
    
**********************************/
}

collision_integral::~collision_integral(){
    delete eps;
    
    int N_outer = outer_dummy_vars->get_length();
    
    for(int i = 0; i < N_outer; i++){
        delete inner_dummy_vars[i];
        delete inner_vals[i];
        
        for(int j = 0; j < 4 * num_F; j++)
            delete F_values[j][i];
    }
    
    for(int j = 0; j < 4 * num_F; j++)
        delete[] F_values[j];
        
    delete[] F_values;
    
    delete[] inner_dummy_vars;
    delete[] inner_vals;
    
    delete outer_dummy_vars;
    delete outer_vals;
}

int collision_integral::get_bin()
{   return bin;}

int collision_integral::get_num_bins()
{   return N_bins;}

double collision_integral::get_eps_value()
{   return eps_value;}

bool collision_integral::is_neutrino()
{   return neutrino;}

void collision_integral::set_min_rate(density* dens){
    double R_vals[2];
    compute_R(dens->get_Tcm(), dens->get_T(), R_vals);
    
    double C_vals[4];
    whole_integral(dens, C_vals, false);
    
    double c_max = max(C_vals[0], C_vals[1]);
    for(int i = 2; i < 4; i++)
        c_max = max(c_max, C_vals[i]);
        
    min_rate = max(R_vals[0], R_vals[1]) * c_max * TOLERANCE_MIN_RATE;
}

double collision_integral::get_min_rate()
{   return min_rate;}

void collision_integral::get_inner_matrix(density* dens, double nrg, sub_dummy_vars* sdv, int inner_index, bool nu, matrix* result, bool use_matrix){ 
    if(sdv->get_need_interp(inner_index)){
        three_vector* A = new three_vector();
        double A0 = dens->interpolated_matrix(nu, sdv->get_interp_index(inner_index), nrg, A);
        if(use_matrix)
            result->convert_p_to_matrix(A0, A);
        else
            result->convert_p_to_identity_minus_matrix(A0, A);
        delete A;
    }
    else{
        if(use_matrix)
            result->convert_p_to_matrix(dens, nu, inner_index);
        else
            result->convert_p_to_identity_minus_matrix(dens, nu, inner_index);
    }
}

void collision_integral::get_inner_matrix(density* dens, double nrg, int outer_index, int inner_index, bool nu, matrix* result, bool use_matrix){ 
    get_inner_matrix(dens, nrg, inner_dummy_vars[outer_index], inner_index, nu, result, use_matrix);
}


electron_collision_integral::electron_collision_integral(int b, linspace_and_gl* e, bool nu, double T_cm) : collision_integral(b, e, nu){
    Tcm = T_cm;

    scaled_me = _electron_mass_ / Tcm;
    me_squared = scaled_me * scaled_me;

    complex_three_vector* A = new complex_three_vector();
    A->set_value(2, complex<double> (0.5,0));
    G_L = new matrix(complex<double> (_sin_squared_theta_W_,0),A);
    delete A;
   
    G_R = new matrix(true);
    G_R->multiply_by(_sin_squared_theta_W_);
}

electron_collision_integral::~electron_collision_integral(){
    delete G_L;
    delete G_R;
    
    int N_outer = outer_dummy_vars_2->get_length();
    
    for(int i = 0; i < N_outer; i++){
        delete inner_dummy_vars_2[i];
        delete inner_vals_2[i];
    }
        
    delete[] inner_dummy_vars_2;
    delete[] inner_vals_2;
    
    delete outer_dummy_vars_2;
    delete outer_vals_2;
    
}

double electron_collision_integral::get_Tcm()
{   return Tcm;}

void electron_collision_integral::set_Tcm(double T_cm)
{   Tcm = T_cm;}

double electron_collision_integral::mom_to_eps(double q){
    return sqrt(q*q + me_squared);
}

double electron_collision_integral::eps_to_mom(double e){
    if(e < scaled_me){
        cout << "ERROR: TRYING TO USE ELECTRON ENERGY LESS THAN MASS" << endl;
        return -1.;
    }
    return sqrt(e*e - me_squared);
}

void electron_collision_integral::get_inner_matrix(density* dens, double nrg, int outer_index, int inner_index, bool nu, matrix* result, bool use_matrix, bool inner_dv1){
    if(inner_dv1)
        get_inner_matrix(dens, nrg, inner_dummy_vars[outer_index], inner_index, nu, results, use_matrix);
    else
        get_inner_matrix(dens, nrg, inner_dummy_vars_2[outer_index], inner_index, nu, results, use_matrix);    
}

nu_nu_collision::nu_nu_collision(int b, linspace_and_gl* e, bool nu) : collision_integral(b, e, nu){
    outer_dummy_vars = new dummy_vars(e);
    outer_vals = new dep_vars(N_bins);
    
    inner_dummy_vars = new sub_dummy_vars*[N_bins];
    inner_vals = new dep_vars*[N_bins];

    interpolation_indices = new int*[N_bins];
    
    p4_values = new sub_dummy_vars*[N_bins];    

    num_F = 2;
    F_values = new double**[4*num_F];
    for(int j = 0; j < 4*num_F; j++)
        F_values[j] = new double*[N_bins];

    int p4_len;
    for(int i = 0; i < N_bins; i++){
        inner_dummy_vars[i] = new sub_dummy_vars(e);
        inner_vals[i] = new dep_vars(N_bins);
        
        interpolation_indices[i] = new int[N_bins];
        
        p4_len = 0;
        for(int j = 0; j < N_bins; j++){
            if(eps_value + outer_dummy_vars->get_value(i) - inner_dummy_vars[i]->get_value(j) >= 0)
                p4_len++;
            else
                break;
        }
        
        
        
        for(int j = 0; j < 4*num_F; j++)
            F_values[j][i] = new double[N_bins]();
    }
    
    double p2_energy;
    double p3_energy;
    double p4_energy;
    
    for(int p2=0; p2<outer_dummy_vars->get_len(); p2++){
        p2_energy = outer_dummy_vars->get_value(p2);
        for(int p3=0; p3<inner_dummy_vars[p2]->get_len(); p3++){
            p3_energy = inner_dummy_vars[p2]->get_value(p3);
            p4_energy = eps_value + p2_energy - p3_energy;
            if(p4_energy>=0){
                interpolation_indices[p2][p3] = eps->index_below_for_interpolation(p4_energy);            
            }
        }
    }
}

nu_nu_collision::~nu_nu_collision(){
    for(int i = 0; i < N_bins; i++)
        delete interpolation_indices[i];
    delete[] interpolation_indices;
}

int nu_nu_collision::estimate_load(){
    int iter = 0;
    
    for(int i = 0; i < outer_dummy_vars->get_length(); i++)
        iter += inner_dummy_vars[i]->bin_below(eps_value + outer_dummy_vars->get_value(i));
        
    return iter * load_factor;

}

/**********************************

NOTE: populate_F is built upon p1, p2, p3 all being at the bins (no interpolation needed)

**********************************/

void nu_nu_collision::populate_F(density* dens, bool net){    
    Fvvsc_for_p1(dens, net);
    Fvvbarsc_for_p1(dens, net);
}

void nu_nu_collision::Fvvsc_for_p1(density* dens, bool net){
    double F0 = 0;
    three_vector* Fxyz = new three_vector();
    for(int p2=0; p2<outer_dummy_vars->get_len(); p2++){
        for(int p3=0; p3<inner_dummy_vars[p2]->get_len(); p3++){
            //this clause means F is filled in only if p3_energy is less than p1_energy+p2_energy
            //this won't affect results for p3 objects that go up to p1+p2, and this is useful in trap rule when we want to set integrand to 0 past p1+p2
            if(inner_dummy_vars[p2]->get_value(p3) <= eps_value + outer_dummy_vars->get_value(p2)){
                //if p4_energy>=0
                if(eps_value+outer_dummy_vars->get_value(p2)-inner_dummy_vars[p2]->get_value(p3)>=0){
                    Fvvsc_components(dens, p2, p3, &F0, Fxyz, net);

                    //factor of 1/4 corrects for use of Froustey's matrix form statistical factor in BURST integral
                    F_values[0][p2][p3] = F0;
                    F_values[1][p2][p3] = Fxyz->get_value(0);
                    F_values[2][p2][p3] = Fxyz->get_value(1);
                    F_values[3][p2][p3] = Fxyz->get_value(2);
                }
            }
        }
    }  
    delete Fxyz;
}

void nu_nu_collision::Fvvbarsc_for_p1(density* dens, bool net){
    double F0 = 0;
    three_vector* Fxyz = new three_vector();
    for(int p2=0; p2<outer_dummy_vars->get_len(); p2++){
        for(int p3=0; p3<inner_dummy_vars[p2]->get_len(); p3++){
            //only if p3_energy is less than p1_energy+p2_energy is this called--> will make integrand 0 past p1+p2 for trap 
            //won't affect results for other dummy var objects
            if(inner_dummy_vars[p2]->get_value(p3) <= eps_value + outer_dummy_vars->get_value(p2)){
                //p4_energy must be >=0
                if(eps_value + outer_dummy_vars->get_value(p2) - inner_dummy_vars[p2]->get_value(p3)>=0){
                    Fvvbarsc_components(dens, p2, p3, &F0, Fxyz, net);

                    F_values[4][p2][p3] = F0;
                    F_values[5][p2][p3] = Fxyz->get_value(0);
                    F_values[6][p2][p3] = Fxyz->get_value(1);
                    F_values[7][p2][p3] = Fxyz->get_value(2);

                }
            }
        }
    }  
    delete Fxyz;
}



void nu_nu_collision::Fvvsc_components(density* dens, int p2, int p3, double* F03, three_vector* F3, bool net){
    double F01;
    three_vector* F1 = new three_vector();
    double F02;
    three_vector* F2 = new three_vector();


    Fvvsc_components_term_1(dens, p2, p3, &F01, F1);
    Fvvsc_components_term_2(dens, p2, p3, &F02, F2);

    if(net==true){
        F2->multiply_by(-1);
        *F03 = F01 - F02;
    }
    else{
        *F03 = F01 + F02;
    }
    F3->add(F1, F2);

    delete F1;
    delete F2;
}

void nu_nu_collision::Fvvbarsc_components(density* dens, int p2, int p3, double* F03, three_vector* F3, bool net){
    double F01;
    three_vector* F1 = new three_vector();
    double F02;
    three_vector* F2 = new three_vector();

    Fvvbarsc_components_term_1(dens, p2, p3, &F01, F1);
    Fvvbarsc_components_term_2(dens, p2, p3, &F02, F2);

    if(net==true){
        F2->multiply_by(-1);
        *F03 = F01 - F02;
    }
    else{
        *F03 = F01 + F02;
    }

    F3->add(F1, F2);


    delete F1;
    delete F2;
}



void nu_nu_collision::Fvvsc_components_term_1(density* dens, int p2, int p3, double* F0, three_vector* F){
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();
    matrix* p_4 = new matrix();
    p_1->convert_p_to_identity_minus_matrix(dens, neutrino, p1);
    p_2->convert_p_to_identity_minus_matrix(dens, neutrino, p2);
    p_3->convert_p_to_matrix(dens, neutrino, p3);

    double p4_energy = eps_value + outer_dummy_vars->get_value(p2) - inner_dummy_vars[p2]->get_value(p3);
    if(p4_energy<0){
        p4_energy = 0;
    }
    
    three_vector* A = new three_vector();
    double A0 = dens->interpolated_matrix(neutrino, interpolation_indices[p2][p3], p4_energy, A);
    p_4->convert_p_to_matrix(A0,A);
    delete A;

    /*
       p_1 = 1-rho_1
       p_2 = 1-rho_2
       p_3 = rho_3
       p_4 = rho_4

       F_dummy1 = (1-rho_2)(rho_4)
       id = tr((1-rho_2)(rho_4)) = 2*A_0*1
       F_dummy2 = F_dummy1 + id
       F_dummy3 = (rho_3) * F_dummy2
       F_dummy4 = (1-rho_1) * F_dummy3
       F_dummy4 is the second term in 2.53 (F_sc) in Froustey
       
    */

    matrix* F_dummy1 = new matrix();
    F_dummy1->matrix_multiply(p_2, p_4);

    matrix* id = new matrix(true);
    id->multiply_by(F_dummy1->get_A0()*(complex<double> (2,0)));

    matrix* F_dummy2 = new matrix();
    F_dummy2->matrix_add(F_dummy1, id);

    matrix* F_dummy3 = new matrix();
    F_dummy3->matrix_multiply(p_3, F_dummy2);

    matrix* F_dummy4 = new matrix();
    F_dummy4->matrix_multiply(p_1, F_dummy3);

    complex<double> comp_F0 = F_dummy4->get_A0();
    complex_three_vector* comp_F = F_dummy4->get_A();

    comp_F->multiply_by(2);

    *F0 = 2*real(comp_F0);
    F->make_real(comp_F);

    delete F_dummy1;
    delete F_dummy2;
    delete id;
    delete F_dummy3;
    delete F_dummy4;
    delete p_1;
    delete p_2;
    delete p_3;
    delete p_4;

}


void nu_nu_collision::Fvvsc_components_term_2(density* dens, int p2, int p3, double* F0, three_vector* F){
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();

    matrix* p_4 = new matrix(true);
    p_1->convert_p_to_matrix(dens, neutrino, p1);
    p_2->convert_p_to_matrix(dens, neutrino, p2);
    p_3->convert_p_to_identity_minus_matrix(dens, neutrino, p3);
    
    double p4_energy = eps_value + outer_dummy_vars->get_value(p2) - inner_dummy_vars[p2]->get_value(p3);
    if(p4_energy<0){
        p4_energy = 0;
    }

    three_vector* A = new three_vector();
    double A0 = dens->interpolated_matrix(neutrino, interpolation_indices[p2][p3], p4_energy, A);
    p_4->convert_p_to_identity_minus_matrix(A0,A);
    delete A;

    /*
       p_1 = rho_1
       p_2 = rho_2
       p_3 = 1-rho_3
       p_4 = 1-rho_4

       F_dummy1 = (rho_2)(1-rho_4)
       id = tr((rho_2)(1-rho_4)) = 2*A_0*1
       F_dummy2 = F_dummy1 + id
       F_dummy3 = (1-rho_3) * F_dummy2
       F_dummy4 = (rho_1) * F_dummy3
       F_dummy4 is the fourth term in 2.53 (F_sc) in Froustey

    */

    matrix* F_dummy1 = new matrix();
    F_dummy1->matrix_multiply(p_2, p_4);

    matrix* id = new matrix(true);
    id->multiply_by(F_dummy1->get_A0()*(complex<double> (2,0)));

    matrix* F_dummy2 = new matrix();
    F_dummy2->matrix_add(F_dummy1, id);

    matrix* F_dummy3 = new matrix();
    F_dummy3->matrix_multiply(p_3,F_dummy2);

    matrix* F_dummy4 = new matrix();
    F_dummy4->matrix_multiply(p_1,F_dummy3); 


    complex<double> comp_F0 = F_dummy4->get_A0();
    complex_three_vector* comp_F = F_dummy4->get_A();

    comp_F->multiply_by(2);

    *F0 = 2*real(comp_F0);
    F->make_real(comp_F);

    delete F_dummy1;
    delete F_dummy2;
    delete id;
    delete F_dummy3;
    delete F_dummy4;
    delete p_1;
    delete p_2;
    delete p_3;
    delete p_4;
}

void nu_nu_collision::Fvvbarsc_components_term_1(density* dens, int p2, int p3, double* F0, three_vector* F){
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();

    matrix* p_4 = new matrix(true);
    p_1->convert_p_to_identity_minus_matrix(dens, neutrino, p1);
    p_2->convert_p_to_identity_minus_matrix(dens, not neutrino, p2);
    p_3->convert_p_to_matrix(dens, not neutrino, p3);
    
    double p4_energy = eps_value + outer_dummy_vars->get_value(p2) - inner_dummy_vars[p2]->get_value(p3);
    if(p4_energy<0){
        p4_energy = 0;
    }

    three_vector* A = new three_vector();
    double A0 = dens->interpolated_matrix(neutrino, interpolation_indices[p2][p3], p4_energy, A);
    p_4->convert_p_to_matrix(A0,A);
    delete A;

    /*
       F_dummy1 = (rho_3)(1-rho_2)
       id1 = 1*tr((rho_3)(1-rho_2))
       F_dummy2 = (rho_3)(1-rho_2)+1*tr((rho_3)(1-rho_2))
       F_dummy3 = (1-rho_1)(rho_4)
       F_dummy4 = (1-rho_1)(rho_4) * [(rho_3)(1-rho_2)+1*tr((rho_3)(1-rho_2))]
       F_dummy4 is the second term in 2.54 (F_sc) in Froustey

       F_dummy5 = (rho_3)(rho_4)
       id2 = 1*tr((rho_3)(rho_4))
       F_dummy6 = (rho_3)(rho_4)+1*tr((rho_3)(rho_4))
       F_dummy7 = (1-rho_1)(1-rho_2)
       F_dummy8 = (1-rho_1)(1-rho_2) * [(rho_3)(rho_4)+1*tr((rho_3)(rho_4))]
       F_dummy8 is the second term in 2.55 (F_ann) in Froustey

       F_dummy9 = F_dummy4+F_dummy8
    */

    matrix* F_dummy1 = new matrix();
    F_dummy1->matrix_multiply(p_3, p_2);

    matrix* id1 = new matrix(true);
    id1->multiply_by(F_dummy1->get_A0()*(complex<double> (2,0)));

    matrix* F_dummy2 = new matrix();
    F_dummy2->matrix_add(F_dummy1, id1);

    matrix* F_dummy3 = new matrix();
    F_dummy3->matrix_multiply(p_1, p_4);

    matrix* F_dummy4 = new matrix();
    F_dummy4->matrix_multiply(F_dummy3, F_dummy2);

    matrix* F_dummy5 = new matrix();
    F_dummy5->matrix_multiply(p_3, p_4);

    matrix* id2 = new matrix(true);
    id2->multiply_by(F_dummy5->get_A0()*(complex<double> (2,0)));

    matrix* F_dummy6 = new matrix();
    F_dummy6->matrix_add(F_dummy5, id2);

    matrix* F_dummy7 = new matrix();
    F_dummy7->matrix_multiply(p_1, p_2);

    matrix* F_dummy8 = new matrix();
    F_dummy8->matrix_multiply(F_dummy7, F_dummy6);

    matrix* F_dummy9 = new matrix();
    F_dummy9->matrix_add(F_dummy4, F_dummy8);

    complex<double> comp_F0 = F_dummy9->get_A0();
    complex_three_vector* comp_F = F_dummy9->get_A();

    comp_F->multiply_by(2);

    *F0 = 2*real(comp_F0);
    F->make_real(comp_F);
    

    delete F_dummy1;
    delete F_dummy2;
    delete F_dummy3;
    delete F_dummy4;
    delete F_dummy5;
    delete F_dummy6;
    delete F_dummy7;
    delete F_dummy8;
    delete F_dummy9;
    delete id1;
    delete id2;
    delete p_1;
    delete p_2;
    delete p_3;
    delete p_4;    
}

void nu_nu_collision::Fvvbarsc_components_term_2(density* dens, int p2, int p3, double* F0, three_vector* F){
    matrix* p_1 = new matrix();
    matrix* p_2 = new matrix();
    matrix* p_3 = new matrix();

    matrix* p_4 = new matrix(true);
    p_1->convert_p_to_matrix(dens, neutrino, p1);
    p_2->convert_p_to_matrix(dens, not neutrino, p2);
    p_3->convert_p_to_identity_minus_matrix(dens, not neutrino, p3);
    
    double p4_energy = eps_value + outer_dummy_vars->get_value(p2) - inner_dummy_vars[p2]->get_value(p3);
    if(p4_energy<0){
        p4_energy = 0;
    }

    three_vector* A = new three_vector();
    double A0 = dens->interpolated_matrix(neutrino, interpolation_indices[p2][p3], p4_energy, A);
    p_4->convert_p_to_identity_minus_matrix(A0,A);
    delete A;

    /*
       F_dummy1 = (1-rho_3)(rho_2)
       id1 = 1*tr((1-rho_3)(rho_2))
       F_dummy2 = (1-rho_3)(rho_2)+1*tr((1-rho_3)(rho_2))
       F_dummy3 = (rho_1)(1-rho_4)
       F_dummy4 = (rho_1)(1-rho_4) * [(1-rho_3)(rho_2)+1*tr((1-rho_3)(rho_2))]
       F_dummy4 is the fourth term of 2.54 (F_sc) in Froustey

       F_dummy5 = (1-rho_3)(1-rho_4)
       id2 = 1*tr((1-rho_3)(1-rho_4))
       F_dummy6 = (1-rho_3)(1-rho_4)+1*tr((1-rho_3)(1-rho_4))
       F_dummy7 = (rho_1)(rho_2)
       F_dummy8 = (rho_1)(rho_2) * [(1-rho_3)(1-rho_4)+1*tr((1-rho_3)(1-rho_4))]
       F_dummy 8 is the fourth term of 2.55 (F_ann) in Froustey

       F_dummy9 = F_dummy4+F_dummy8
    */

    matrix* F_dummy1 = new matrix();
    F_dummy1->matrix_multiply(p_3, p_2);

    matrix* id1 = new matrix(true);
    id1->multiply_by(F_dummy1->get_A0()*(complex<double> (2,0)));

    matrix* F_dummy2 = new matrix();
    F_dummy2->matrix_add(F_dummy1, id1);

    matrix* F_dummy3 = new matrix();
    F_dummy3->matrix_multiply(p_1, p_4);

    matrix* F_dummy4 = new matrix();
    F_dummy4->matrix_multiply(F_dummy3, F_dummy2);

    matrix* F_dummy5 = new matrix();
    F_dummy5->matrix_multiply(p_3, p_4);

    matrix* id2 = new matrix(true);
    id2->multiply_by(F_dummy5->get_A0()*(complex<double> (2,0)));

    matrix* F_dummy6 = new matrix();
    F_dummy6->matrix_add(F_dummy5, id2);

    matrix* F_dummy7 = new matrix();
    F_dummy7->matrix_multiply(p_1, p_2);

    matrix* F_dummy8 = new matrix();
    F_dummy8->matrix_multiply(F_dummy7, F_dummy6);

    matrix* F_dummy9 = new matrix();
    F_dummy9->matrix_add(F_dummy4, F_dummy8);

    complex<double> comp_F0 = F_dummy9->get_A0();
    complex_three_vector* comp_F = F_dummy9->get_A();

    comp_F->multiply_by(2);

    *F0 = 2*real(comp_F0);
    F->make_real(comp_F);

    delete F_dummy1;
    delete F_dummy2;
    delete F_dummy3;
    delete F_dummy4;
    delete F_dummy5;
    delete F_dummy6;
    delete F_dummy7;
    delete F_dummy8;
    delete F_dummy9;
    delete id1;
    delete id2;
    delete p_1;
    delete p_2;
    delete p_3;
    delete p_4;
}

double nu_nu_collision::interior_integral(int p2, int which_term){
    inner_vals[p2]->zeros();
    
    for(int p3 = 0; p3 < inner_vals[p2]->get_length(); p3++)
        inner_vals[p2]->set_value(p3, 0.25 * J(eps_value, outer_dummy_vars->get_value(p2), inner_dummy_vars[p2]->get_value(p3)) * F_values[which_term][p2][p3]
                                    + 0.125 * K(eps_value, outer_dummy_vars->get_value(p2), inner_dummy_vars[p2]->get_value(p3)) * F_values[which_term+4][p2][p3]);
    
    return inner_dummy_vars[p2]->integrate(inner_vals[p2]);
}

double nu_nu_collision::J(double p1, double p2, double p3){
    if(p2 < p1){
        if(p3 < p2)
            return J1(p1, p2, p3);
        else if (p3 < p1)
            return J2(p1, p2);
        else if (p3 < p1 + p2)
            return J3(p1, p2, p3);
        else
            return 0;
    }
    else{
        if(p3 < p1)
            return J1(p1, p2, p3);
        else if (p3 < p2)
            return J2(p2, p1);
        else if (p3 < p1 + p2)
            return J3(p1, p2, p3);
        else
            return 0;
    }
}

double nu_nu_collision::K(double p1, double p2, double p3){
    if(p2 < p1){
        if(p3 < p2)
            return K1(p1, p3);
        else if (p3 < p1)
            return K2(p1, p2, p3);
        else if (p3 < p1 + p2)
            return K3(p1, p2, p3);
        else
            return 0;
    }
    else{
        if(p3 < p1)
            return K1(p1, p3);
        else if (p3 < p2)
            return K1(p3, p1);
        else if (p3 < p1 + p2)
            return K3(p1, p2, p3);
        else
            return 0;
    }
}

double J1(double p1, double p2, double p3){
    return 16./15 * pow(p3,3) * (10 * pow(p1+p2,2) - 15 * (p1+p2) * p3 + 6 * pow(p3,2));  
}

double J2(double p1, double p2){
    return 16./15 * pow(p2,3) * (10 * pow(p1,2) + 5 * p1*p2 + pow(p2,2));  
}

double J3(double p1, double p2, double p3){
    return 16./15 * (pow(p1+p2,5) - 10 * pow(p1+p2,2) * pow(p3, 3) + 15 * (p1+p2) * pow(p3,4) - 6 * pow(p3,5));
}

double K1(double p1, double p3){
    return 16./15 * pow(p3,3) * (10 * pow(p1,2) - 5 * p1*p3 + pow(p3,2));
}

double K2(double p1, double p2, double p3){
    return 16./15 * pow(p2,3) * (10 * pow(p1-p3,2) + 15 * (p1-p3) * p2 + 6 * pow(p2,2));
}

double K3(double p1, double p2, double p3){
    return 16./15 * (pow(p1-p3,5) + 10 * pow(p1-p3,2) * pow(p2,3) + 15 * (p1-p3) * pow(p2,4) + 6 * pow(p2,5));
}


void nu_nu_collision::whole_integral(density* dens, double* results, bool net){
    if (eps_value == 0)
        for(int j = 0; j < 4; j++)
            results[j] = 0.;
    else{
        populate_F(dens, net);
        
        double coeff = pow(dens->get_Tcm(), 5) * pow(_GF_, 2) / (pow(2*_PI_, 3) * pow(eps_value,2 ));
        
        for(int j = 0; j < 4; j++){
            for(int p2 = 0; p2 < N_bins; p2++)
                outer_vals->set_value(p2, interior_integral(p2, j));
                
            results[j] = coeff * outer_dummy_vars->integrate(outer_vals);
            if (abs(results[j]) < min_rate)
                results[j] = 0.;
        }
    }
}

void nu_nu_collision::compute_R(double Tcm, double T, double* results){
    density* thermal = new density(eps, 0., 0.);
    thermal->set_T_Tcm(T, Tcm);
    
    double net[4];
    double frs[4];
    
    whole_integral(thermal, net, true);
    whole_integral(thermal, frs, false);
    
    for(int j = 0; j < 4; j+=3){
        if (frs[j] == 0)
            results[j%2] = _COMPUTE_R_ERROR_;
        else
            results[j%2] = abs(net[j] / frs[j]);
       }
    
    delete thermal;
}


nu_e_collision::nu_e_collision(int b, linspace_and_gl* e, bool nu, double T_cm) : electron_collision_integral(b, e, nu, T_cm){
    num_F = 4;
    F_values = new double**[4*num_F];

    int N_outer = 50;
    outer_dummy_vars = new gl_dummy_vars(N_outer, 0.);
    outer_vals = new dep_vars(N_outer);
    
    inner_dummy_vars = new dummy_vars*[N_outer];
    inner_vals = new dep_vars*[N_outer];
    
    double epslim[5];
    double eps2, eps_low_lim, p4_low, p4_high;
    
    int count_min, count_max, p4_len;
    
    int bot_shift, top_shift;
    
    interpolation_indices = new int**[2];
    
    interpolation_indices[0] = new int*[N_outer];
    
    for(int i = 0; i < N_outer; i++){
        epslim_R1(outer_dummy_vars->get_value(i), epslim);
        eps2 = mom_to_eps(outer_dummy_vars->get_value(i));
        
        eps_low_lim = scaled_me;
        if(eps_value < 0.5 * scaled_me && eps2 > epslim[EPS_CUT_1])
            eps_low_lim = epslim[EPS_LIM_2];
            
        p4_low = eps_value + eps2 - epslim[EPS_LIM_1];
        p4_high = eps_value + eps2 - eps_low_lim;
        
        inner_dummy_vars[i] = eps->new_inner_dummy_vars(p4_low, p4_high);
        
        

        inner_vals[i] = new dep_vars(inner_dummy_vars[i]->get_length());
        
        for(int j = 0; j < 8; j++){
            F_values[j] = new double*[N_outer];
            for(int k = 0; k < N_outer; k++)
                F_values[j][k] = new double[inner_vals[i]->get_length()]();
        }
    }
    
    outer_dummy_vars_2 = new gl_dummy_vars(N_outer, 0.);
    outer_vals_2 = new dep_vars(N_outer);
    
    inner_dummy_vars_2 = new dummy_vars*[N_outer];
    inner_vals_2 = new dep_vars*[N_outer];
    
    double epslim_2[6];
    double eps3;
    
    for(int i = 0; i < N_outer; i++){
        epslim_R2(outer_dummy_vars_2->get_value(i), epslim_2);
        eps3 = mom_to_eps(outer_dummy_vars_2->get_value(i));
        
        eps_low_lim = 0.;
        if(eps3 > epslim_2[EPS_CUT_2])
            eps_low_lim = epslim_2[EPS_LIM_2];

        p4_low = eps_value + eps_low_lim - eps3;
        
        if(eps_value < 0.5 * scaled_me && eps3 < epslim_2[EPS_CUT_1])
            p4_high = eps_value + epslim_2[EPS_LIM_1] - eps3;
        else
            p4_high = INNER_INTEGRAL_INFINITY;

        if(p4_high == INNER_INTEGRAL_INFINITY)
            inner_dummy_vars_2[i] = eps->new_inner_dummy_vars_infinite(p4_low, eps->get_num_gl());
        else
            inner_dummy_vars_2[i] = eps->new_inner_dummy_vars(p4_low, p4_high);
            
        inner_vals_2[i] = new dep_vars(inner_dummy_vars_2[i]->get_length());
        
        for(int j = 8; j < 16; j++){
            F_values[j] = new double*[N_outer];
            for(int k = 0; k < N_outer; k++)
                F_values[j][k] = new double[inner_vals_2[i]->get_length()]();
        }
    }
}

void nu_e_collision::epslim_R1(double q2, double* eps_lim){
    double eps_cut_1, eps_cut_3, eps_trans_2, eps_lim_1, eps_lim_2;
    eps_cut_1 = scaled_me + 2 * eps_value * eps_value / (scaled_me - 2 * eps_value);
    eps_cut_3 = sqrt(eps_value*eps_value + me_squared);

    double eps2 = mom_to_eps(q2);

    eps_trans_2 = 0.5*(2 * eps_value + eps2 - q2 + me_squared/(2*eps_value + eps2 - q2));
    eps_lim_1 = 0.5 * (2 * eps_value + eps2 + q2 + me_squared / (2 * eps_value + eps2 + q2));
    eps_lim_2 = eps_trans_2;

    eps_lim[EPS_CUT_1] = eps_cut_1;
    eps_lim[EPS_CUT_3] = eps_cut_3;
    eps_lim[EPS_TRANS_2] = eps_trans_2;
    eps_lim[EPS_LIM_1] = eps_lim_1;
    eps_lim[EPS_LIM_2] = eps_lim_2;
}

void nu_e_collision::epslim_R2(double q3, double* eps_lim){
    double eps_cut_1, eps_cut_2, eps_cut_3, eps_trans_2, eps_lim_1, eps_lim_2;

    eps_cut_1 = eps_value + me_squared / (4 * eps_value);
    eps_cut_2 = eps_value + scaled_me * (eps_value + scaled_me) / (2 * eps_value + scaled_me);
    eps_cut_3 = mom_to_eps(eps_value);
    
    double eps3 = mom_to_eps(q3);

    eps_trans_2 = 0.5 * (eps3 + q3 - 2 * eps_value + me_squared / (eps3 + q3 - 2 * eps_value));
    
    eps_lim_1 = 0.5 * (eps3 - q3 - 2 * eps_value + me_squared / (eps3 - q3 - 2 * eps_value));
    eps_lim_2 = eps_trans_2;
    
    eps_lim[EPS_CUT_1] = eps_cut_1;
    eps_lim[EPS_CUT_2] = eps_cut_2;
    eps_lim[EPS_CUT_3] = eps_cut_3;
    eps_lim[EPS_TRANS_2] = eps_trans_2;
    eps_lim[EPS_LIM_1] = eps_lim_1;
    eps_lim[EPS_LIM_2] = eps_lim_2;
}

int nu_nu_collision::estimate_load(){
    int iter = 0;
    
    for(int i = 0; i < outer_dummy_vars->get_length(); i++)
        iter += inner_dummy_vars[i]->get_length() * 11;
        
    for(int i = 0; i < outer_dummy_vars_2->get_length(); i++)
        iter += inner_dummy_vars_2[i]->get_length() * 11;
        
    return iter;
}

void nu_e_collision::populate_F(density* dens, bool net){
    F_LL_RR_for_p1(dens, net);
    F_LR_RL_for_p1(dens, net);
}

void nu_e_collision::F_LL_RR_for_p1(density* dens, bool net){
    double F0 = 0.;
    three_vector* Fxyz = new three_vector();
    
    for(int p2 = 0; p2 < outer_dummy_vars->get_length(); p2++)
}





