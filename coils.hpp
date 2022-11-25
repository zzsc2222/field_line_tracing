#include "shared.hpp"

class coils {
    public:

    int num_of_coils, num_of_rows, resolution;
    int* separation = nullptr;
    double* current = nullptr;
    double** cart_coor = nullptr;
    double** cyl_coor = nullptr; // R zeta Z
    double** ex_cart_coor = nullptr;
    double** ex_cyl_coor = nullptr;
    double* ex_current = nullptr;
    int* ex_separation = nullptr;

    coils(double**, int, int, double);
    // coils(const coils &);
    ~coils();

    private:
    
    double** catmull_rom = nullptr;
    
};

double** cart_to_cyl(double**, int);
double* cart_to_cyl(double*);
double** cyl_to_cart(double**, int);
double* cyl_to_cart(double*);
double** catmull_rom_spline(double**, int, double, int*);
double** cr_fitting(double**, double**, int, int, int, int*);
double* cr_fitting(double*, int, int, int, int*);