#include "../test_extrapolation.hh"

#include "../code/include.hh"

int main(int argc, char* argv[]){

    ofstream file;
    file.open(argv[1]);

    ofstream orig;
    orig.open(argv[2]);
    
    linspace_and_gl* eps = new linspace_and_gl(0., EPS_MAX_LINSPACE, EPS_LINSPACE_POINTS, 5);
    density* model = new density(eps->get_length(), eps);

    linspace_for_trap* bins = new linspace_for_trap(0., 2 * eps->get_max_value(), 1000);
    density* interp_model = new density(bins->get_length(), bins);

    bins->print_csv(file);
    file << endl;

    eps->print_csv(orig);
    file << endl;
    
    int index;
    three_vector* P0P = new three_vector();
    
    for(int i = 0; i < Num_cases; i++){
        for(int j = 0; j < 8*eps->get_length()+2; j++)
            model->set_value(j, density_RK[i][j]);
        for(int k = 0; k < bins->get_length(); k++){
            index = eps->index_below_for_interpolation(bins->get_value(k));
            interp_model->set_value(k, true, 0, model->interpolated_matrix(true, index, bins->get_value(k), P0P));
            for(int jj=0; jj < 3; jj++)
                interp_model->set_value(k, true, jj+1, P0P->get_value(jj)/interp_model->p0(k, true));

            interp_model->set_value(k, false, 0, model->interpolated_matrix(false, index, bins->get_value(k), P0P));
            for(int jj=0; jj < 3; jj++)
                interp_model->set_value(k, false, jj+1, P0P->get_value(jj)/interp_model->p0(k, false));
        }

        file << i / 6 << ", " << i % 6 << ", ";
        interp_model->print_csv(file);
        file << endl;

        orig << i/6 << ", " << i % 6 << ", ";
        model->print_csv(orig);
        orig << endl;
    }

    file.close();
    orig.close();

    delete P0P;
    delete interp_model;
    delete bins;
    delete model;
    delete eps;
    
        
    return 0;
}
