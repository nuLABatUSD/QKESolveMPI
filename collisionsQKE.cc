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
    
    if (eps_value == 0.)
        return;
    
    min_rate = 0.;
    
/**********************************

    outer_vals and outer_dummy_vars; inner_vals and inner_dummy_vars; F_values matrix are allocated in the individual collision integrals
    BUT they are deleted in the base class, because they all need to be deleted. Only arrays particular to integrals should be deleted in the individual destructors
    
**********************************/
}

collision_integral::~collision_integral(){
    delete eps;

    if (eps_value == 0.)
        return;
    
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
//    cout << "call base get_inner_matrix" << endl;
//    cout << "* " << sdv->get_need_interp(inner_index) << endl;
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
//        cout << "** " << sdv->get_interp_index(inner_index) << ", " << inner_index << ", " << sdv->get_length() << endl;
        if(use_matrix)
            result->convert_p_to_matrix(dens, nu, sdv->get_interp_index(inner_index));
        else
            result->convert_p_to_identity_minus_matrix(dens, nu, sdv->get_interp_index(inner_index));
    }
//    cout << "get_inner_matrix base done" << endl;
}

void collision_integral::get_inner_matrix(density* dens, double nrg, int outer_index, int inner_index, bool nu, matrix* result, bool use_matrix){ 
    get_inner_matrix(dens, nrg, inner_dummy_vars[outer_index], inner_index, nu, result, use_matrix);
}


electron_collision_integral::electron_collision_integral(int b, linspace_and_gl* e, bool nu, double T_cm) : collision_integral(b, e, nu){
    if (eps_value == 0.)
         return;
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
    if (eps_value == 0.)
         return;
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
    if (e > scaled_me)
        return sqrt(e*e - me_squared);
    else if (e - scaled_me > - 1.e-12 * scaled_me)
        return 0.;
    else{
        cout << "eps_to_mom:  ERROR: TRYING TO USE ELECTRON ENERGY LESS THAN MASS; eps_value = " << eps_value << endl;
        cout << "   e = " << e << ", m/Tcm = " << scaled_me << "; e - m/Tcm = " << e - scaled_me << endl;
        return -1.;
    }
}

void electron_collision_integral::get_inner_matrix(density* dens, double nrg, int outer_index, int inner_index, bool nu, matrix* result, bool use_matrix, bool inner_dv1){
//    cout << "get_inner_matrix" << endl;
    if(inner_dv1)
        get_inner_matrix(dens, nrg, inner_dummy_vars[outer_index], inner_index, nu, result, use_matrix);
    else
        get_inner_matrix(dens, nrg, inner_dummy_vars_2[outer_index], inner_index, nu, result, use_matrix);    
//    cout << "get_inner_matrix done" << endl;
}

void collision_integral::show_idv(int info, int outer){
    sub_dummy_vars* sdv = inner_dummy_vars[outer];
    for(int i = 0; i < sdv->get_length(); i++){
        if(info == 0)
            cout << sdv->get_value(i);
        else if(info == 1)
            cout << sdv->get_need_interp(i);
        else
            cout << sdv->get_interp_index(i);
        cout << ", ";
    }
    cout << endl;
}

void electron_collision_integral::show_idv(int R, int info, int outer){
    if(R == 1)
        collision_integral::show_idv(info, outer);
    else{
        sub_dummy_vars* sdv = inner_dummy_vars_2[outer];
        
        for(int i = 0; i < sdv->get_length(); i++){
            if(info == 0)
                cout << sdv->get_value(i);
            else if(info == 1)
                cout << sdv->get_need_interp(i);
            else
                cout << sdv->get_interp_index(i);
            cout << ", ";
        }
        cout << endl;
    }
}



nu_nu_collision::nu_nu_collision(int b, linspace_and_gl* e, bool nu) : collision_integral(b, e, nu){
    if (eps_value == 0.)
         return;
    outer_dummy_vars = new dummy_vars(e);
    outer_vals = new dep_vars(N_bins);
    
    inner_dummy_vars = new sub_dummy_vars*[N_bins];
    inner_vals = new dep_vars*[N_bins];
    
    p4_values = new sub_dummy_vars*[N_bins];    

    num_F = 2;
    F_values = new double**[4*num_F];
    for(int j = 0; j < 4*num_F; j++)
        F_values[j] = new double*[N_bins];

    int p4_len;
    for(int i = 0; i < N_bins; i++){
        inner_dummy_vars[i] = new sub_dummy_vars(e);
        inner_vals[i] = new dep_vars(N_bins);
        
        p4_len = 0;
        for(int j = 0; j < N_bins; j++){
            if(eps_value + outer_dummy_vars->get_value(i) - inner_dummy_vars[i]->get_value(j) >= 0)
                p4_len++;
            else
                break;
        }
        
        if(p4_len > 0){
            p4_values[i] = new sub_dummy_vars(e, p4_len);
            for(int j = 0; j < p4_len; j++)
                p4_values[i]->set_value(j, eps_value + outer_dummy_vars->get_value(i) - inner_dummy_vars[i]->get_value(j));
            p4_values[i]->set_interp();
        }
        else
            p4_values[i] = nullptr;
        
        for(int j = 0; j < 4*num_F; j++)
            F_values[j][i] = new double[N_bins]();
    }
}

nu_nu_collision::~nu_nu_collision(){
    if (eps_value == 0.)
         return;
    for(int i = 0; i < N_bins; i++){
        if(!p4_values[i])
            delete p4_values[i];
    }
}

int nu_nu_collision::estimate_load(){
    if (eps_value == 0)
        return 0;
        
    int iter = 0;
    
    for(int i = 0; i < outer_dummy_vars->get_length(); i++)
        iter += inner_dummy_vars[i]->bin_below(eps_value + outer_dummy_vars->get_value(i));
        
    return iter * load_factor;

}

void nu_nu_collision::get_p4_matrix(density* dens, int p2, int p3, bool nu, matrix* result, bool use_matrix){
//    cout << "* " << p4_values[p2]->get_value(p3) << ", " << p4_values[p2]->get_interp_index(p3) << "*";
    get_inner_matrix(dens, p4_values[p2]->get_value(p3), p4_values[p2], p3, nu, result, use_matrix);
}

/**********************************

NOTE: populate_F is built upon p1, p2, p3 all being at the bins (no interpolation needed)
*** 9/26: addressed: using get_inner_matrix() to generalize p3, (interpolation or not)

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

    get_inner_matrix(dens, inner_dummy_vars[p2]->get_value(p3), p2, p3, neutrino, p_3, true);
    get_p4_matrix(dens, p2, p3, neutrino, p_4, true);
    
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
    
    get_inner_matrix(dens, inner_dummy_vars[p2]->get_value(p3), p2, p3, neutrino, p_3, false);

    get_p4_matrix(dens, p2, p3, neutrino, p_4, false);
    

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
    
    get_inner_matrix(dens, inner_dummy_vars[p2]->get_value(p3), p2, p3, not neutrino, p_3, true);

    get_p4_matrix(dens, p2, p3, neutrino, p_4, true);
    
    
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
    
    get_inner_matrix(dens, inner_dummy_vars[p2]->get_value(p3), p2, p3, not neutrino, p_3, false);

    get_p4_matrix(dens, p2, p3, neutrino, p_4, false);
    

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
    if (eps_value == 0.)
         return;
    num_F = 4;
    F_values = new double**[4*num_F];

    int N_outer = 50;
    outer_dummy_vars = new gl_dummy_vars(N_outer, 0.);
    outer_vals = new dep_vars(N_outer);
    
    inner_dummy_vars = new sub_dummy_vars*[N_outer];
    inner_vals = new dep_vars*[N_outer];
    
    R_elec_values = new double***[2];
    epslim_values = new double**[2];
    
    for(int i = 0; i < 2; i++)
        R_elec_values[i] = new double**[N_outer];
    
    double* epslim;
    double eps2, eps3, eps_low_lim, p4_low, p4_high;
    
    int count_min, count_max, p4_len;
    
    int bot_shift, top_shift;
        
    epslim_values[0] = new double*[N_outer];
    for(int i = 0; i < N_outer; i++){
        epslim_values[0][i] = new double[5];
        epslim_R1(outer_dummy_vars->get_value(i), epslim_values[0][i]);
        epslim = epslim_values[0][i];
        
        eps2 = mom_to_eps(outer_dummy_vars->get_value(i));
        
        eps_low_lim = scaled_me;
        if(eps_value < 0.5 * scaled_me && eps2 > epslim[EPS_CUT_1])
            eps_low_lim = epslim[EPS_LIM_2];
            
        p4_low = eps_value + eps2 - epslim[EPS_LIM_1];
        p4_high = eps_value + eps2 - eps_low_lim;
        
//        cout << "nu_e_collision, " << eps_value << ", " << p4_low << ", " << p4_high << endl;
        
        inner_dummy_vars[i] = new sub_dummy_vars(eps, p4_low, p4_high, eps->get_num_gl());     
                
//        cout << "nu_e_collision, " << inner_dummy_vars[i]->length() << endl;   

        inner_vals[i] = new dep_vars(inner_dummy_vars[i]->get_length());
        
        R_elec_values[0][i] = new double*[inner_dummy_vars[i]->get_length()];
        for(int k = 0; k < inner_dummy_vars[i]->get_length(); k++){
            eps3 = eps_value + eps2 - inner_dummy_vars[i]->get_value(k);
            
            if(eps3 < scaled_me)
                eps3 = scaled_me;
        
            R_elec_values[0][i][k] = new double[4];
            
            R_elec_values[0][i][k][NUE_Q2] = outer_dummy_vars->get_value(i);
            R_elec_values[0][i][k][NUE_E2] = eps2;
            R_elec_values[0][i][k][NUE_E3] = eps3;
            R_elec_values[0][i][k][NUE_Q3] = eps_to_mom(eps3);
        }
        for(int j = 0; j < 8; j++){
            F_values[j] = new double*[N_outer];
            for(int k = 0; k < N_outer; k++)
                F_values[j][k] = new double[inner_vals[i]->get_length()]();
        }
    }
    
//    cout << "nu_e_collision, switch to R2" << endl;
    
    outer_dummy_vars_2 = new gl_dummy_vars(N_outer, 0.);
    outer_vals_2 = new dep_vars(N_outer);
    
    inner_dummy_vars_2 = new sub_dummy_vars*[N_outer];
    inner_vals_2 = new dep_vars*[N_outer];
    
    double* epslim_2;
    
    epslim_values[1] = new double*[N_outer];
    for(int i = 0; i < N_outer; i++){
        epslim_values[1][i] = new double[6];
        epslim_R2(outer_dummy_vars_2->get_value(i), epslim_values[1][i]);
        epslim_2 = epslim_values[1][i];
        
        eps3 = mom_to_eps(outer_dummy_vars_2->get_value(i));
        
        eps_low_lim = scaled_me;
        if(eps3 > epslim_2[EPS_CUT_2])
            eps_low_lim = epslim_2[EPS_LIM_2];

        p4_low = eps_value + eps_low_lim - eps3;
        
        if(eps_value < 0.5 * scaled_me && eps3 < epslim_2[EPS_CUT_1])
            p4_high = eps_value + epslim_2[EPS_LIM_1] - eps3;
        else
            p4_high = INNER_INTEGRAL_INFINITY;

        if(p4_high == INNER_INTEGRAL_INFINITY)
            inner_dummy_vars_2[i] = new sub_dummy_vars(eps, p4_low, INNER_INTEGRAL_INFINITY, eps->get_num_gl());   
        else
            inner_dummy_vars_2[i] = new sub_dummy_vars(eps, p4_low, p4_high, eps->get_num_gl());   
            
//        cout << "nu_e_collision, " << inner_dummy_vars_2[i]->get_length() << endl;

        inner_vals_2[i] = new dep_vars(inner_dummy_vars_2[i]->get_length());
        
        R_elec_values[1][i] = new double*[inner_dummy_vars_2[i]->get_length()];
        for(int k = 0; k < inner_dummy_vars_2[i]->get_length(); k++){
            eps2 = inner_dummy_vars_2[i]->get_value(k) + eps3 - eps_value;
        
            R_elec_values[1][i][k] = new double[4];
            
            R_elec_values[1][i][k][NUE_Q2] = eps_to_mom(eps2);
            R_elec_values[1][i][k][NUE_E2] = eps2;
            R_elec_values[1][i][k][NUE_E3] = eps3;
            R_elec_values[1][i][k][NUE_Q3] = outer_dummy_vars_2->get_value(i);
                if (R_elec_values[1][i][k][NUE_Q2] == -1.){
                    cout << "***" << p4_low << ", " << p4_high << "," << eps3 << ", " <<eps_value << endl;
                    cout << "*** " << i << ", " << k << ", " << inner_dummy_vars_2[i]->get_value(k) << endl;
                    }
        }
        for(int j = 8; j < 16; j++){
            F_values[j] = new double*[N_outer];
            for(int k = 0; k < N_outer; k++)
                F_values[j][k] = new double[inner_vals_2[i]->get_length()]();
        }
    }
    
}

nu_e_collision::~nu_e_collision(){
    if (eps_value == 0.)
         return;
        
    for(int i = 0; i < outer_dummy_vars->get_length(); i++){
        for(int k = 0; k < inner_dummy_vars[i]->get_length(); k++)
            delete[] R_elec_values[0][i][k];
        delete[] R_elec_values[0][i];
        
        delete[] epslim_values[0][i];
    }
    
    for(int i = 0; i < outer_dummy_vars_2->get_length(); i++){
        for(int k = 0; k < inner_dummy_vars_2[i]->get_length(); k++)
            delete[] R_elec_values[1][i][k];
        delete[] R_elec_values[1][i];
        
        delete[] epslim_values[1][i];
    }
    
    delete[] epslim_values;
    delete[] R_elec_values;
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


int nu_e_collision::estimate_load(){
    if (eps_value == 0.)
         return 0;

    int iter = 0;
    
    for(int i = 0; i < outer_dummy_vars->get_length(); i++)
        iter += inner_dummy_vars[i]->get_length() * 11;
        
    for(int i = 0; i < outer_dummy_vars_2->get_length(); i++)
        iter += inner_dummy_vars_2[i]->get_length() * 11;
        
    return iter;
}

void nu_e_collision::populate_F(density* dens, bool net){
//    cout << "populate_F, R1" << endl;
    F_R1_for_p1(dens, net);
//    cout << "populate_F, R2" << endl;
    F_R2_for_p1(dens, net);
//    cout << "populate_F, done" << endl;
    return;
}

// Interior integrals of R1
double nu_e_collision::interior_integral(int q2, int which_term){
   //limits of integration from BURST C5a
   //coefficients on integrand come from adjusting factor needed when using Froustey's statistical factor in BURST integral
   
    inner_vals[q2]->zeros();
   
    double q2_momentum = R_elec_values[0][q2][0][NUE_Q2];
    double E2 = R_elec_values[0][q2][0][NUE_E2];
    double q3_momentum, E3;
    int term = 0;
    int outer_branch = 0;
    
    if(eps_value < 0.5 * scaled_me){
       //case 1a
        if(E2 < epslim_values[0][q2][EPS_CUT_3])
            outer_branch = 0;
       //case 1b
        else if(E2 < epslim_values[0][q2][EPS_CUT_1])
            outer_branch = 1;
       //case 1c
        else
            outer_branch = 2;
    }
    else{
       //case 2a
        if(E2 < epslim_values[0][q2][EPS_CUT_3])
            outer_branch = 3;
       //case 2b
        else
            outer_branch = 4;
    }
    
    for(int p4 = 0; p4 < inner_dummy_vars[q2]->get_length(); p4++){
        q3_momentum = R_elec_values[0][q2][p4][NUE_Q3];
        E3 = R_elec_values[0][q2][p4][NUE_E3];
        switch(outer_branch){
       //case 1a
            case(0):
               //case 1ai
                if(E3 < E2)
                    term = 1;
               //case 1aii
                else if(E3 < epslim_values[0][q2][EPS_TRANS_2])
                    term = 2;
               //case 1aiii
                else
                    term = 3;
                break;
       //case 1b
            case(1):
               //case 1bi
                if(E3 < epslim_values[0][q2][EPS_TRANS_2])
                    term = 1;
               //case 1bii
                else if(E3 < E2)
                    term = 4;
               //case 1biii
                else
                    term = 3;
                break;
       //case 1c
            case(2):
               //case 1ci
                if(E3 < E2)
                    term = 4;
               //case 1cii
                else
                    term = 3;
                break;
       //case 2a
            case(3):
               //case 2ai
                if(E3 < E2)
                    term = 1;
               //case 2aii
                else if (E3 < epslim_values[0][q2][EPS_TRANS_2])
                    term = 2;
               //case 2aiii
                else
                    term = 3;
                break;
       //case 2b
            case(4):
               //case 2bi
                if(E3 < epslim_values[0][q2][EPS_TRANS_2])
                    term = 1;
               //case 2bii
                else if(E3 < E2)
                    term = 4;
               //case 2biii
                else
                    term = 3;
                break;
        }
        inner_vals[q2]->set_value(p4, 2 * F_values[which_term][q2][p4] * M_12(term, R_elec_values[0][q2][p4]) - 0.5 * me_squared * F_values[4+which_term][q2][p4] * M_11(term, R_elec_values[0][q2][p4]));
    }
    
    return inner_dummy_vars[q2]->integrate(inner_vals[q2]);
}

// Interior integrals of R2
double nu_e_collision::interior_integral_2(int q3, int which_term){
   //integration limits found in BURST C5b
   //coefficients on integrand are adjusting factors from using Froustey's statistical factor with BURST integral
   
    inner_vals_2[q3]->zeros();
   
    double q3_momentum = R_elec_values[1][q3][0][NUE_Q3];
    double E3 = R_elec_values[1][q3][0][NUE_E3];
    double q2_momentum, E2;
    
    int term = 0;
    int outer_branch = 0;

   //case 1    if(p1_me < (sqrt(5)-1)/4.)
    if(eps_value < (sqrt(5)-1)/4. * scaled_me){
       //case 1a if(q3_momentum < q_cut_3)
        if(E3 < epslim_values[1][q3][EPS_CUT_3])
            outer_branch = 0;
       //case 1b else if(q3_momentum < q_cut_2_R2)
        else if(E3 < epslim_values[1][q3][EPS_CUT_2])
            outer_branch = 1;
       //case 1c else if(q3_momentum < q_cut_1_R2)
        else if(E3 < epslim_values[1][q3][EPS_CUT_1])
            outer_branch = 2;
       //case 1d
        else
            outer_branch = 3; 
    }
   //case 2 else if(p1_me < 1./(2 * sqrt(2)))
    else if(eps_value < 1./(2 * sqrt(2)) * scaled_me){
       //case 2a if(q3_momentum<q_cut_3)
        if(E3 < epslim_values[1][q3][EPS_CUT_3])
            outer_branch = 4;
       //case 2b  else if(q3_momentum<q_cut_1_R2)
        else if(E3 < epslim_values[1][q3][EPS_CUT_1])
            outer_branch = 5;
       //case 2c else if(q3_momentum<q_cut_2_R2)
        else if(E3 < epslim_values[1][q3][EPS_CUT_2])
            outer_branch = 6;
       //case 2d
        else
            outer_branch = 7;        
    }
   //case 3 else if(p1_me < 0.5)
    else if (eps_value < 0.5 * scaled_me){
       //case 3a if(q3_momentum<q_cut_1_R2)
        if(E3 < epslim_values[1][q3][EPS_CUT_1])
            outer_branch = 8;
       //case 3b else if(q3_momentum<q_cut_3)
        else if(E3 < epslim_values[1][q3][EPS_CUT_3])
            outer_branch = 9;
       //case 3c else if(q3_momentum<q_cut_2_R2)
        else if(E3 < epslim_values[1][q3][EPS_CUT_2])
            outer_branch = 10;
       //case 3d
        else
            outer_branch = 11; 
    }
   //case 4
    else{
        //case 4a if(q3_momentum<q_cut_1_R2)
        if(E3 < epslim_values[1][q3][EPS_CUT_1])
            outer_branch = 12;
       //case 4b  else if(q3_momentum<q_cut_3)
        else if(E3 < epslim_values[1][q3][EPS_CUT_3])
            outer_branch = 13;
       //case 4c else if(q3_momentum<q_cut_2_R2)
        else if(E3 < epslim_values[1][q3][EPS_CUT_2])
            outer_branch = 14;
       //case 4d
        else
            outer_branch = 15; 
    }
    
    
    for(int p4 = 0; p4 < inner_dummy_vars_2[q3]->get_length(); p4++){
        q2_momentum = R_elec_values[1][q3][p4][NUE_Q2];
        E2 = R_elec_values[1][q3][p4][NUE_E2];
        switch(outer_branch){
        //CASE 1
            case(0):
                if(E2 < E3)
                    term = 1;
                else if(E2 < epslim_values[1][q3][EPS_TRANS_2])
                    term = 2;
                else
                    term = 3;
                break;
            case(1):
                if(E2 < epslim_values[1][q3][EPS_TRANS_2])
                    term = 1;
                else if(E2 < E3)
                    term = 4;
                else
                    term = 3;
                break;
            case(2):
                if(E2 < E3)
                    term = 4;
                else
                    term = 3;
                break;
            case(3):
                if(E2 < E3)
                    term = 4;
                else
                    term = 3;
                break;
        //CASE 2
            case(4):
                if(E2 < E3)
                    term = 1;
                else if(E2 < epslim_values[1][q3][EPS_TRANS_2])
                    term = 2;
                else
                    term = 3;
                break;
            case(5):
                if(E2 < epslim_values[1][q3][EPS_TRANS_2])
                    term = 1;
                else if(E2 < E3)
                    term = 4;
                else
                    term = 3;
                break;
            case(6):
                if(E2 < epslim_values[1][q3][EPS_TRANS_2])
                    term = 1;
                else if(E2 < E3)
                    term = 4;
                else
                    term = 3;
                break;
            case(7):
                if(E2 < E3)
                    term = 4;
                else
                    term = 3;
                break;
        //CASE 3
            case(8):
                if(E2 < E3)
                    term = 1;
                else if(E2 < epslim_values[1][q3][EPS_TRANS_2])
                    term = 2;
                else
                    term = 3;
                break;
            case(9):
                if(E2 < E3)
                    term = 1;
                else if(E2 < epslim_values[1][q3][EPS_TRANS_2])
                    term = 2;
                else
                    term = 3;
                break;
            case(10):
                if(E2 < epslim_values[1][q3][EPS_TRANS_2])
                    term = 1;
                else if(E2 < E3)
                    term = 4;
                else
                    term = 3;
                break;
            case(11):
                if(E2 < E3)
                    term = 4;
                else
                    term = 3;
                break;
        //CASE 4
            case(12):
                if(E2 < E3)
                    term = 1;
                else
                    term = 2;
                break;
            case(13):
                if(E2 < E3)
                    term = 1;
                else if(E2 < epslim_values[1][q3][EPS_TRANS_2])
                    term = 2;
                else
                    term = 3;
                break;
            case(14):
                if(E2 < epslim_values[1][q3][EPS_TRANS_2])
                    term = 1;
                else if(E2 < E3)
                    term = 4;
                else
                    term = 3;
                break;
            case(15):
                if(E2 < E3)
                    term = 4;
                else
                    term = 3;
                break;
        }
        inner_vals_2[q3]->set_value(p4, 0.5 * F_values[8+which_term][q3][p4] * M_22(term, R_elec_values[1][q3][p4]) + 0.5 * me_squared * F_values[12+which_term][q3][p4] * M_21(term, R_elec_values[1][q3][p4]));
    }
    return inner_dummy_vars_2[q3]->integrate(inner_vals_2[q3]);
}
void nu_e_collision::whole_integral(density* dens, double* results, bool net){
    if (eps_value == 0)
        for(int j = 0; j < 4; j++)
            results[j] = 0.;
    else{
//        cout << "whole integral, populate F" << endl;
        populate_F(dens, net);
        
        double coeff = pow(Tcm, 5) * pow(_GF_, 2) / (16 * pow(2*_PI_, 3) * pow(eps_value,2 ));
        
        for(int i = 0; i < 4; i++){
//            cout << "whole integral, R1" << endl;
            for(int q2 = 0; q2 < outer_dummy_vars->get_length(); q2++)
                outer_vals->set_value(q2, R_elec_values[0][q2][0][NUE_Q2] / R_elec_values[0][q2][0][NUE_E2] * interior_integral(q2, i));
            results[i] = outer_dummy_vars->integrate(outer_vals); //******* only doing R2 for now
            
//            cout << "whole integral, R2" << endl;
            
            for(int q3 = 0; q3 < outer_dummy_vars_2->get_length(); q3++)
                outer_vals_2->set_value(q3, R_elec_values[1][q3][0][NUE_Q3] / R_elec_values[1][q3][0][NUE_E3] * interior_integral_2(q3, i));
                
            results[i] += outer_dummy_vars_2->integrate(outer_vals_2);
            
            results[i] *= coeff;
        }
    }
    return;
}

void nu_e_collision::compute_R(double Tcm, double T, double* results){
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

    return;

}


void nu_e_collision::F_R1_for_p1(density* dens, bool net){
    double F0 = 0;
    three_vector* Fxyz = new three_vector();
        
    for(int q2 = 0; q2 < outer_dummy_vars->get_length(); q2++){
        for(int p4 = 0; p4 < inner_dummy_vars[q2]->get_length(); p4++){            
            F_LL_F_RR(dens, true, R_elec_values[0][q2][p4][NUE_E2], R_elec_values[0][q2][p4][NUE_E3], q2, p4, net, &F0, Fxyz);
            
            F_values[0][q2][p4] = F0;
            F_values[1][q2][p4] = Fxyz->get_value(0);
            F_values[2][q2][p4] = Fxyz->get_value(1);
            F_values[3][q2][p4] = Fxyz->get_value(2);
            
            F_LR_F_RL(dens, true, R_elec_values[0][q2][p4][NUE_E2], R_elec_values[0][q2][p4][NUE_E3], q2, p4, net, &F0, Fxyz);
            
            F_values[4][q2][p4] = F0;
            F_values[5][q2][p4] = Fxyz->get_value(0);
            F_values[6][q2][p4] = Fxyz->get_value(1);
            F_values[7][q2][p4] = Fxyz->get_value(2);
        }
    }
}

void nu_e_collision::F_R2_for_p1(density* dens, bool net){
    double F0 = 0;
    three_vector* Fxyz = new three_vector();
    
    for(int q3 = 0; q3 < outer_dummy_vars_2->get_length(); q3++){
        for(int p4 = 0; p4 < inner_dummy_vars_2[q3]->get_length(); p4++){
//            cout << "F_R2, " << p4 << " / " << inner_dummy_vars_2[q3]->get_length() << endl;
//            cout << "*F_R2, " << R_elec_values[1][q3][p4][NUE_E2] << ", " <<  R_elec_values[1][q3][p4][NUE_E3] << endl;
            
            F_LL_F_RR(dens, false, R_elec_values[1][q3][p4][NUE_E2], R_elec_values[1][q3][p4][NUE_E3], q3, p4, net, &F0, Fxyz);
//            cout << "F_R2, F_LLRR done" << endl;
            
            F_values[8][q3][p4] = F0;
            F_values[9][q3][p4] = Fxyz->get_value(0);
            F_values[10][q3][p4] = Fxyz->get_value(1);
            F_values[11][q3][p4] = Fxyz->get_value(2);
            
            
            F_LR_F_RL(dens, false, R_elec_values[1][q3][p4][NUE_E2], R_elec_values[1][q3][p4][NUE_E3], q3, p4, net, &F0, Fxyz);
            
//            cout << "F_R2, F_LRRL done" << endl;

            F_values[12][q3][p4] = 2 * F0;
            F_values[13][q3][p4] = 2 * Fxyz->get_value(0);
            F_values[14][q3][p4] = 2 * Fxyz->get_value(1);
            F_values[15][q3][p4] = 2 * Fxyz->get_value(2);
        }
    }
}



void nu_e_collision::F_LL_F_RR(density* dens, bool is_R1, double E2, double E3, int outer_index, int inner_index, bool net, double* F0, three_vector* F){
    double T__Tcm = dens->get_T() / Tcm;

    double f2 = 1 / (exp(E2/T__Tcm)+1);
    double f3 = 1 / (exp(E3/T__Tcm)+1);
    
    matrix* p_1 = new matrix();
    matrix* minus_p_1 = new matrix();
    matrix* p_4 = new matrix();
    
    p_1->convert_p_to_matrix(dens, neutrino, p1);
    minus_p_1->convert_p_to_identity_minus_matrix(dens, neutrino, p1);
    
    double p4_energy;
    if(is_R1)
        p4_energy = inner_dummy_vars[outer_index]->get_value(inner_index);
    else
        p4_energy = inner_dummy_vars_2[outer_index]->get_value(inner_index);
        
/*    if (!is_R1 && eps_value + E2 != (E3 + p4_energy))
        cout << "FLL, Delta E = " << eps_value + E2 - (E3 + p4_energy) << endl; */
    
    get_inner_matrix(dens, p4_energy, outer_index, inner_index, neutrino, p_4, true, is_R1);
    
    matrix* minus_p_4 = new matrix(p_4);
    minus_p_4->convert_this_to_identity_minus_this();
    
    double fwd = real(p_1->get_A0()) * f2 * (1-f3) * real(minus_p_4->get_A0());
    double rev = real(p_4->get_A0()) * f3 * (1-f2) * real(minus_p_1->get_A0());
    
    if (!is_R1 && abs(1-fwd/rev) > 1.e-13 && (fwd > 1.e-10 || rev > 1.e-1) && false){
        cout << E2 << "\t" << E3 << "\t" << p4_energy << "\t" << fwd << "\t" << rev << "\t" << 1 - fwd/rev << "; " << f2 << ",\t" << f3 << "\t" << p_4->get_A0() << endl;
        cout << "**" << p4_energy << "\t" << inner_dummy_vars_2[outer_index]->get_value(0) << "\t" << inner_dummy_vars_2[outer_index]->get_value(inner_dummy_vars_2[outer_index]->get_length()-1) << endl;
//        cout << inner_dummy_vars_2[outer_index]->get_need_interp(inner_index) << "\t" << inner_dummy_vars_2[outer_index]->get_interp_index(inner_index) << endl;

        for(int j = 0; j < inner_dummy_vars_2[outer_index]->get_length(); j++)
            cout << inner_dummy_vars_2[outer_index]->get_interp_index(j) << ", ";
        cout << endl;
    }
    
    matrix* F_dummy1 = new matrix();
    matrix* F_dummy2 = new matrix();
    matrix* F_dummy3 = new matrix();
    matrix* F_dummy4 = new matrix();
    matrix* F_dummy5 = new matrix();
    matrix* F_dummy6 = new matrix();
    matrix* F_dummy7 = new matrix();
    matrix* F_dummy8 = new matrix();
    matrix* F_dummy9 = new matrix();
    matrix* F_dummy10 = new matrix();
    matrix* F_dummy11 = new matrix();
    
    F_dummy1->matrix_multiply(G_L, p_4);
    F_dummy2->matrix_multiply(G_L, minus_p_1);
    F_dummy3->matrix_multiply(F_dummy1, F_dummy2);
    F_dummy3->multiply_by(f3 * (1-f2));
    F_dummy4->matrix_multiply(G_L, minus_p_4);
    F_dummy5->matrix_multiply(G_L, p_1);
    F_dummy6->matrix_multiply(F_dummy4, F_dummy5);
    F_dummy6->multiply_by(f2 * (1-f3));
    if(net==true){
        F_dummy6->multiply_by(complex<double> (-1,0));
    }
    F_dummy7->matrix_add(F_dummy3, F_dummy6);
    
    F_dummy8->matrix_multiply(p_4, minus_p_1);
    F_dummy8->multiply_by(f3 * (1-f2));
    F_dummy9->matrix_multiply(minus_p_4, p_1);
    F_dummy9->multiply_by(f2 * (1-f3));
    if(net==true){
        F_dummy9->multiply_by(complex<double> (-1,0));
    }
    F_dummy10->matrix_add(F_dummy8, F_dummy9);
    
    F_dummy10->multiply_by(pow(_sin_squared_theta_W_,2));
     
    F_dummy11->matrix_add(F_dummy7, F_dummy10);
     
    complex<double> comp_F0 = F_dummy11->get_A0();
    complex_three_vector* comp_F = F_dummy11->get_A();
    
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
    delete F_dummy10;
    delete F_dummy11;
    delete p_4;
    delete minus_p_4;
    delete p_1;
    delete minus_p_1;  
}

void nu_e_collision::F_LR_F_RL(density* dens, bool is_R1, double E2, double E3, int outer_index, int inner_index, bool net, double* F0, three_vector* F){
    double T__Tcm = dens->get_T() / Tcm;

    double f2 = 1 / (exp(E2/T__Tcm)+1);
    double f3 = 1 / (exp(E3/T__Tcm)+1);
    
    matrix* p_1 = new matrix();
    matrix* minus_p_1 = new matrix();
    matrix* p_4 = new matrix();
    
    p_1->convert_p_to_matrix(dens, neutrino, p1);
    minus_p_1->convert_p_to_identity_minus_matrix(dens, neutrino, p1);
    
    double p4_energy;
    if(is_R1)
        p4_energy = inner_dummy_vars[outer_index]->get_value(inner_index);
    else
        p4_energy = inner_dummy_vars_2[outer_index]->get_value(inner_index);

    get_inner_matrix(dens, p4_energy, outer_index, inner_index, neutrino, p_4, true, is_R1);
    
    matrix* minus_p_4 = new matrix(p_4);
    
    minus_p_4->convert_this_to_identity_minus_this();
    
/*    if (!is_R1 && eps_value + E2 != (E3 + p4_energy))
        cout << "FLR, Delta E = " << eps_value + E2 - (E3 + p4_energy) << endl;*/
    
//    cout << "FLR, " << minus_p_4->get_A0() << ", " << p_4->get_A0() << endl;
    
    /*
    F_dummy1 = G_L * rho_4
    F_dummy2 = F_dummy1 * (1-rho_1)
    F_dummy3 = rho_4 * G_L
    F_dummy4 = F_dummy3 * (1-rho_1)
    
    F_dummy5 = G_L * (1-rho_4)
    F_dummy6 = F_dummy5 * rho_1
    F_dummy7 = (1-rho_4) * G_L
    F_dummy8 = F_dummy7 * rho_1
    
    F_dummy9 = F_dummy2 + F_dummy4 => F_dummy9
    F_dummy10 = F_dummy6 + F_dummy8 => F_dummy10
    F_dummy11 = F_dummy9 - F_dummy10 => F_dummy11 = F_LR + F_RL
    
    F_LR = F_dummy2 - F_dummy_6
    F_RL = F_dummy4 - F_dummy_8
    see Froustey C.19
    */
    
    matrix* F_dummy1 = new matrix();
    matrix* F_dummy2 = new matrix();
    matrix* F_dummy3 = new matrix();
    matrix* F_dummy4 = new matrix();
    matrix* F_dummy5 = new matrix();
    matrix* F_dummy6 = new matrix();
    matrix* F_dummy7 = new matrix();
    matrix* F_dummy8 = new matrix();
    matrix* F_dummy9 = new matrix();
    matrix* F_dummy10 = new matrix();
    matrix* F_dummy11 = new matrix();
    
    F_dummy1->matrix_multiply(G_L, p_4);
    F_dummy2->matrix_multiply(F_dummy1, minus_p_1);
    F_dummy3->matrix_multiply(p_4, G_L);
    F_dummy4->matrix_multiply(F_dummy3, minus_p_1);
    
    F_dummy5->matrix_multiply(G_L, minus_p_4);
    F_dummy6->matrix_multiply(F_dummy5, p_1);
    F_dummy7->matrix_multiply(minus_p_4, G_L);
    F_dummy8->matrix_multiply(F_dummy7, p_1);
    
    F_dummy9->matrix_add(F_dummy2, F_dummy4);
    F_dummy9->multiply_by(f3 * (1-f2));
    F_dummy10->matrix_add(F_dummy6, F_dummy8);
    F_dummy10->multiply_by(f2 * (1-f3));

//    cout << p_4->get_A0() << ", " << F_dummy9->get_A0() << ", " << F_dummy10->get_A0() << ", " << F_dummy9->get_A0() - F_dummy10->get_A0() << endl;

     
    if(net==true){
        F_dummy10->multiply_by(complex<double> (-1,0));
    }
    F_dummy11->matrix_add(F_dummy9, F_dummy10);
    
    F_dummy11->multiply_by(_sin_squared_theta_W_);
    
    complex<double> comp_F0 = F_dummy11->get_A0();
    complex_three_vector* comp_F = F_dummy11->get_A();
    
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
    delete F_dummy10;
    delete F_dummy11;
    delete p_4;
    delete minus_p_4;
    delete p_1;
    delete minus_p_1;  
}

double nu_e_collision::M_11(int which, double* kinematic){
    double E2 = kinematic[NUE_E2];
    double E3 = kinematic[NUE_E3];
    double q2_momentum = kinematic[NUE_Q2];
    double q3_momentum = kinematic[NUE_Q3];

    double b = 0;
    double a = 0;
    double C1 = pow(eps_value + E2,2) - me_squared;
    
    if(which==1){
        a = eps_value + E2 - E3 - q3_momentum;
        b = eps_value + E2 - E3 + q3_momentum;
    }
    else if(which==2){
        a = eps_value - q2_momentum;
        b = eps_value + q2_momentum;
    }
    else if(which==3){
        a = E3 + q3_momentum - eps_value - E2;
        b = eps_value + q2_momentum;
    }
    else{
        a = q2_momentum - eps_value;
        b = eps_value + E2 - E3 + q3_momentum;
    }
    
    return 0.5 * (C1 * b - 1./3 * pow(b,3)) - 0.5 * (C1 * a - 1./3 * pow(a,3));
}

double nu_e_collision::M_12(int which, double* kinematic){
    double E2 = kinematic[NUE_E2];
    double E3 = kinematic[NUE_E3];
    double q2_momentum = kinematic[NUE_Q2];
    double q3_momentum = kinematic[NUE_Q3];

    double b = 0;
    double a = 0;
    double C1 = pow(eps_value + E2,2) - me_squared;
    
    if(which==1){
        a = eps_value + E2 - E3 - q3_momentum;
        b = eps_value + E2 - E3 + q3_momentum;
    }
    else if(which==2){
        a = eps_value - q2_momentum;
        b = eps_value + q2_momentum;
    }
    else if(which==3){
        a = E3 + q3_momentum - eps_value - E2;
        b = eps_value + q2_momentum;
    }
    else{
        a = q2_momentum - eps_value;
        b = eps_value + E2 - E3 + q3_momentum;
    }
    
    return 0.25 * (pow(C1,2) * b - 2./3 * C1 * pow(b,3) + 1./5 * pow(b,5)) - 0.25 * (pow(C1,2) * a - 2./3 * C1 * pow(a,3) + 1./5 * pow(a,5));
   
}

double nu_e_collision::M_21(int which, double* kinematic){
    double E2 = kinematic[NUE_E2];
    double E3 = kinematic[NUE_E3];
    double q2_momentum = kinematic[NUE_Q2];
    double q3_momentum = kinematic[NUE_Q3];
    
    
    double b = 0;
    double a = 0;
    double C2 = pow(eps_value - E3,2) - me_squared;
    
    if(which==1){
        a = eps_value + E2 - E3 - q2_momentum;
        b = eps_value + E2 - E3 + q2_momentum;
    }
    else if(which==2){
        a = eps_value - q3_momentum;
        b = eps_value + q3_momentum;
    }
    else if(which==3){
        a = E3 - eps_value - E2 + q2_momentum;
        b = eps_value + q3_momentum;
    }
    else{
        a = q3_momentum - eps_value;
        b = eps_value + E2 - E3 + q2_momentum;
    }
   
    return 0.5 * (1./3 * pow(b,3) - C2 * b) - 0.5 * (1./3 * pow(a,3) - C2 * a);
}

double nu_e_collision::M_22(int which, double* kinematic){
    double E2 = kinematic[NUE_E2];
    double E3 = kinematic[NUE_E3];
    double q2_momentum = kinematic[NUE_Q2];
    double q3_momentum = kinematic[NUE_Q3];

    double b = 0;
    double a = 0;
    double C2 = pow(eps_value - E3,2) - me_squared;
    
    if(which==1){
        a = eps_value + E2 - E3 - q2_momentum;
        b = eps_value + E2 - E3 + q2_momentum;
    }
    else if(which==2){
        a = eps_value - q3_momentum;
        b = eps_value + q3_momentum;
    }
    else if(which==3){
        a = E3 - eps_value - E2 + q2_momentum;
        b = eps_value + q3_momentum;
    }
    else{
        a = q3_momentum - eps_value;
        b = eps_value + E2 - E3 + q2_momentum;
    }
    
    return 0.25 * (1./5 * pow(b,5) - 2./3 * C2 * pow(b,3) + pow(C2,2) * b) - 0.25 * (1./5 * pow(a,5) - 2./3 * C2 * pow(a,3) + pow(C2,2) * a);
}

