#include "shared.hpp"

char* import_parameters (const char*, const char*);
int num_of_row(const char*);
double** coordinates_initial(const char*, int);
void write_coils_cyl(double**, int);
void write_FLT(double**, double*, int);
void write_flux_iota(double*, double*, int);
int poincare_w(double**, double*, int, double);
