#include "shared.hpp"
#include "coils.cpp"
#include "mag_field.cpp"
#include "read_write.cpp"

int main() {
   
    const char* input_coils = import_parameters("input.txt", "input_coils");

    int num_of_scan = stoi(import_parameters("input.txt", "num_of_scan"));

    int resolution = stoi(import_parameters("input.txt", "resolution"));

    double tension = stod(import_parameters("input.txt", "tension"));

    double init_point[3];
    init_point[0] = stod(import_parameters("input.txt", "init_point[0]"));
    init_point[1] = stod(import_parameters("input.txt", "init_point[1]"));
    init_point[2] = stod(import_parameters("input.txt", "init_point[2]"));

    int num_of_tracing = stoi(import_parameters("input.txt", "num_of_tracing"));

    double step = stod(import_parameters("input.txt", "step"));

    int num_of_section = stoi(import_parameters("input.txt", "num_of_section"));
    double zeta_of_poincare[num_of_section];
  
    for (int i = 0; i < num_of_section; i++) {

        char* str1 = nullptr;
        str1 = new char;
        sprintf(str1,"zeta_of_poincare[%d]\0",i);

        zeta_of_poincare[i] = stod(import_parameters("input.txt", str1));
 
        delete[] str1;  
    }

    int order_of_fitting = stoi(import_parameters("input.txt", "order_of_fitting"));

    double axis_guess_R = stod(import_parameters("input.txt", "axis_guess_R"));

    double axis_guess_zeta = stod(import_parameters("input.txt", "axis_guess_zeta"));

    double axis_guess_Z = stod(import_parameters("input.txt", "axis_guess_Z"));
    
    //////////////////////////////////////////////////////////////////////////

    int num_of_rows = num_of_row(input_coils);

    double** coor_init = nullptr;
    coor_init = new double* [num_of_rows];
    for (int i = 0; i < num_of_rows; i++) {
        coor_init[i] = new double [4];
    }
    coor_init = coordinates_initial(input_coils, num_of_rows);

    coils object_coils(coor_init, num_of_rows, resolution, tension);

    mag_field object_field(object_coils.ex_cart_coor, object_coils.ex_cyl_coor, object_coils.ex_current, \
    object_coils.num_of_rows + resolution * (object_coils.num_of_rows - object_coils.num_of_coils), object_coils.ex_separation, \
    order_of_fitting);

    write_coils_cyl(object_field.ex_cyl_coor, object_field.num_of_rows);
     
    object_field.coef_of_axis = mag_axis_cal(axis_guess_R, axis_guess_zeta, axis_guess_Z, step, num_of_tracing/10, order_of_fitting, \
    object_field.ex_cart_coor, object_field.ex_current, object_field.num_of_rows, object_field.ex_separation, object_coils.num_of_coils);

    //////////////////////////////////////////////////////////////////////////

    double axis_r, axis_z, scanned_point[3];
    for (int j = 1; j < 2 * order_of_fitting + 1; j += 2) {
        axis_r += (object_field.coef_of_axis[0][j] * cos((j+1)/2.0 * init_point[1]) + object_field.coef_of_axis[0][j+1] * sin((j+1)/2.0 * init_point[1]));
        axis_z += (object_field.coef_of_axis[1][j] * cos((j+1)/2.0 * init_point[1]) + object_field.coef_of_axis[1][j+1] * sin((j+1)/2.0 * init_point[1]));
    }
    axis_r += object_field.coef_of_axis[0][0];
    axis_z += object_field.coef_of_axis[1][0];
    
    double flux_t[num_of_scan], iota[num_of_scan];
    int num_of_poincare;
    printf("scanning from outside to inside...\n");
    for (int ns = 0; ns < num_of_scan; ns++) {
        printf("ns = %d/%d...\n",ns,num_of_scan);
        printf("\n");
        scanned_point[0] = ( (num_of_scan - ns) * init_point[0] + ns * axis_r ) / num_of_scan;
        scanned_point[1] = init_point[1];
        scanned_point[2] = ( (num_of_scan - ns) * init_point[2] + ns * axis_z ) / num_of_scan;

        double** trace_points = nullptr;
        trace_points = new double* [num_of_tracing];
        for (int i = 0; i < num_of_tracing; i++) {
        trace_points[i] = new double [3];
        }
        trace_points = adams_bashforth_moulton_cyl(scanned_point, step, num_of_tracing, \
        object_field.ex_cart_coor, object_field.ex_current, object_field.num_of_rows, object_field.ex_separation, object_coils.num_of_coils);

        write_FLT(trace_points, scanned_point, num_of_tracing);

        double** poincare_points = nullptr;
        for (int k = 0; k < num_of_section; k++) {

            num_of_poincare = poincare_w(trace_points, scanned_point, num_of_tracing, zeta_of_poincare[k]);

            poincare_points = new double* [num_of_poincare];
            for (int i = 0; i < num_of_poincare; i++) {
                poincare_points[i] = new double [3];
            }

            poincare_points = poincare(trace_points, num_of_tracing, zeta_of_poincare[k], num_of_poincare);

            if (k < num_of_section - 1) {
                for (int i = 0; i < num_of_poincare; i++) {
                    delete[] poincare_points[i];
                }
            }
        }

        double** coef_of_surface = nullptr;
        coef_of_surface = new double* [2];
        for (int i = 0; i < 2; i++) {
            coef_of_surface[i] = new double [2*order_of_fitting+1];
        }
        coef_of_surface = poincare_fitting(poincare_points, num_of_poincare, order_of_fitting, object_field.coef_of_axis);
        flux_t[ns] = flux_tor_cal(coef_of_surface, order_of_fitting, zeta_of_poincare[num_of_section-1], \
        object_field.ex_cart_coor, object_field.ex_current, object_field.num_of_rows, object_field.ex_separation, object_coils.num_of_coils);

        iota[ns] = iota_cal(trace_points, num_of_tracing, object_field.coef_of_axis, order_of_fitting);

        write_flux_iota(flux_t, iota, num_of_scan);

        for (int i = 0; i < 2; i++) {
            delete[] coef_of_surface[i];
        }
        delete[] coef_of_surface;

        for (int i = 0; i < num_of_poincare; i++) {
            delete[] poincare_points[i];
        }
        delete[] poincare_points;

        for (int i = 0; i < num_of_tracing; i++) {
            delete[] trace_points[i];
        }
        delete[] trace_points;

        printf("finished\n");
        printf("\n");    
    }
    

    for (int i = 0; i < num_of_rows; i++) {
        delete[] coor_init[i];
    }
    delete[] coor_init;
   
    object_field.~mag_field();
    object_coils.~coils();

    return 0 ;
}
