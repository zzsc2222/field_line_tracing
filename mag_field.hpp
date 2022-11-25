#include "shared.hpp"


class mag_field {
    public:
    
    double** ex_cart_coor = nullptr;
    double** ex_cyl_coor = nullptr;
    double** coef_of_axis = nullptr;
    double* ex_current = nullptr;
    int* ex_separation = nullptr;
    int num_of_rows, num_of_tracing, separation;
   
    mag_field(double**, double**, double*, int, int*, int);
    // mag_field(const coils &);
    ~mag_field();
};

double* biot_savart_cyl(double, double, double, double**, double*, int, int*, int);
double* biot_savart_cart(double, double, double, double**, double*, int, int*, int);
double** adams_bashforth_moulton_cyl(double* , double , int , double** , double* , int , int*, int);
double** poincare_fitting(double**, int, int, double**);
double** poincare(double**, int, double, int);
double** mag_axis_cal(double, double, double, double, int, int, double**, double*, int, int*, int);
double flux_tor_cal(double**, int, double, double**, double*, int, int*, int);
double iota_cal(double**, int, double**, int);