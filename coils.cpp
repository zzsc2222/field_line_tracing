#include "coils.hpp"

double** cart_to_cyl(double** var1, int var2) {

    double** array1 = nullptr;
    array1 = new double* [var2];
    for (int i = 0; i < var2; i++) {
        array1[i] = new double [3];
    }

    for (int i = 0; i < var2; i++) {
        array1[i][0] = sqrt( pow(var1[i][0],2) + pow(var1[i][1],2) );
        array1[i][1] = atan2(var1[i][1], var1[i][0]);
        array1[i][2] = var1[i][2];
    }

    return array1;
}
double* cart_to_cyl(double* var1) {

    double* array1 = nullptr;
    array1 = new double [3];

    array1[0] = sqrt( pow(var1[0],2) + pow(var1[1],2) );
    array1[1] = atan2(var1[1], var1[0]);
    array1[2] = var1[2];
    
    return array1;
}

double** cyl_to_cart(double** var1, int var2) {

    double** array1 = nullptr;
    array1 = new double* [var2];
    for (int i = 0; i < var2; i++) {
        array1[i] = new double [3];
    }

    for (int i = 0; i < var2; i++) {
        array1[i][0] = var1[i][0] * cos(var1[i][1]);
        array1[i][1] = var1[i][0] * sin(var1[i][1]);
        array1[i][2] = var1[i][2];
    }

    return array1;
}
double* cyl_to_cart(double* var1) {

    double* array1 = nullptr;
    array1 = new double [3];

    array1[0] = var1[0] * cos(var1[1]);
    array1[1] = var1[0] * sin(var1[1]);
    array1[2] = var1[2];

    return array1;
}

double** catmull_rom_spline(double** var1, int var2, double var3, int* var4) {

    double** array1 = nullptr;
    array1 = new double* [var2];
    for (int i = 0; i < var2; i++) {
        array1[i] = new double [12];
    }

    double dx0 = 0, dx1 = 0, dy0 = 0, dy1 = 0, dz0 = 0, dz1 = 0; 
    int num1 = 0;

    for (int i = 0; i < var2 - 1; i++) {

        if (i == var4[num1] - 1) {
            num1++;
            continue;
        }

        else if (i == 0 || i == var4[num1-1]) { 
            dx0 = var3 * (var1[i][0] - var1[var4[num1]-2][0]);
            dx1 = var3 * (var1[i+2][0] - var1[i+1][0]);
            dy0 = var3 * (var1[i][1] - var1[var4[num1]-2][1]);
            dy1 = var3 * (var1[i+2][1] - var1[i+1][1]);
            dz0 = var3 * (var1[i][2] - var1[var4[num1]-2][2]);
            dz1 = var3 * (var1[i+2][2] - var1[i+1][2]);
        }
        else if (i == var4[0] - 2) {       
            dx0 = var3 * (var1[i][0] - var1[i-1][0]);
            dx1 = var3 * (var1[1][0] - var1[i+1][0]);
            dy0 = var3 * (var1[i][1] - var1[i-1][1]);
            dy1 = var3 * (var1[1][1] - var1[i+1][1]);
            dz0 = var3 * (var1[i][2] - var1[i-1][2]);
            dz1 = var3 * (var1[1][2] - var1[i+1][2]);

        }
        else if (i == var4[num1] - 2) {
            dx0 = var3 * (var1[i][0] - var1[i-1][0]);
            dx1 = var3 * (var1[var4[num1-1]+1][0] - var1[var4[num1-1]][0]);
            dy0 = var3 * (var1[i][1] - var1[i-1][1]);
            dy1 = var3 * (var1[var4[num1-1]+1][1] - var1[var4[num1-1]][1]);
            dz0 = var3 * (var1[i][2] - var1[i-1][2]);
            dz1 = var3 * (var1[var4[num1-1]+1][2] - var1[var4[num1-1]][2]);
        }
        else {
            dx0 = var3 * (var1[i][0] - var1[i-1][0]);
            dx1 = var3 * (var1[i+2][0] - var1[i+1][0]);
            dy0 = var3 * (var1[i][1] - var1[i-1][1]);
            dy1 = var3 * (var1[i+2][1] - var1[i+1][1]);
            dz0 = var3 * (var1[i][2] - var1[i-1][2]);
            dz1 = var3 * (var1[i+2][2] - var1[i+1][2]);
        }
  
        array1[i][0] = var1[i][0];
        array1[i][1] = dx0;
        array1[i][2] = 3*(var1[i+1][0] - var1[i][0]) - (dx1 + 2*dx0);
        array1[i][3] = 2*(var1[i][0] - var1[i+1][0]) + dx1 +dx0;
        array1[i][4] = var1[i][1];
        array1[i][5] = dy0;
        array1[i][6] = 3*(var1[i+1][1] - var1[i][1]) - (dy1 + 2*dy0);
        array1[i][7] = 2*(var1[i][1] - var1[i+1][1]) + dy1 +dy0;
        array1[i][8] = var1[i][2];
        array1[i][9] = dz0;
        array1[i][10] = 3*(var1[i+1][2] - var1[i][2]) - (dz1 + 2*dz0);
        array1[i][11] = 2*(var1[i][2] - var1[i+1][2]) + dz1 +dz0;
    }
    
    return array1;
}

double** cr_fitting(double** var1, double** var2, int var3, int var4, int var5, int* var6) {

    // cart_coor/cyl_coor catmull_rom num_of_rows num_of_coils resolution separation
    double** array1 = nullptr;
    array1 = new double* [var3 + (var3 - var4) * var5];
    for (int i = 0; i < var3 + (var3 - var4) * var5; i++) {
        array1[i] = new double [3];
    }
    
    int num1 = 0, num2 = 0, num4 = 0;
    const double num3 = (const double) 1 / (var5 + 1);
    for (int i = 0; i < var3; i++) {
        if (i == var6[num4] - 1) {
            num4++;
            array1[num2][0] = var1[i][0];
            array1[num2][1] = var1[i][1];
            array1[num2][2] = var1[i][2];
            num2++;
        }
        else {
            for (int j = 0; j < var5 + 1; j++) {
                array1[num2][0] = var2[i][0] + var2[i][1] * pow(num3*j,1) + var2[i][2] * pow(num3*j,2) + var2[i][3] * pow(num3*j,3);
                array1[num2][1] = var2[i][4] + var2[i][5] * pow(num3*j,1) + var2[i][6] * pow(num3*j,2) + var2[i][7] * pow(num3*j,3);
                array1[num2][2] = var2[i][8] + var2[i][9] * pow(num3*j,1) + var2[i][10] * pow(num3*j,2) + var2[i][11] * pow(num3*j,3);
                num2++;
            } 
        }        
    }

    return array1;
}
double* cr_fitting(double* var1, int var2, int var3, int var4, int* var5) {

    // current num_of_rows num_of_coils resolution
    double* array1 = nullptr;
    array1 = new double [var2 + (var2 - var3) * var4];
    
    int num1 = 0, num2 =0;
    for (int i = 0; i < var2; i++) {
        if (i == var5[num2] - 1) {
            num2++;
            array1[num1] = var1[i];
            num1++;
        }
        else {
            for (int j = 0; j < var4 + 1; j++) {
                array1[num1] = var1[i];
                num1++;
            } 
        }   
    }

    return array1;
}

coils::coils(double** var1, int var2, int var3, double var4) {

    printf("constructing and refining coils...\n");

    num_of_coils = (int) var1[var2 + 1][0];
    num_of_rows = var2;
    resolution = var3;

    separation = new int [num_of_coils];
    for (int i = 0; i < num_of_coils; i++) {
        separation[i] = (int) var1[num_of_rows][i];
    }
    
    current = new double [num_of_rows];
    for (int i = 0; i < num_of_rows; i++) {
        current[i] = var1[i][3];
    }
  
    cart_coor = new double* [num_of_rows];
    for (int i = 0; i < num_of_rows; i++) {
        cart_coor[i] = new double [3];  
    }
    for (int i = 0; i < num_of_rows; i++) {
        for (int j = 0; j < 3; j++) {
            cart_coor[i][j] = var1[i][j];
        }
    }
 
    cyl_coor = cart_to_cyl(cart_coor, num_of_rows);
  
    catmull_rom = catmull_rom_spline(cart_coor, num_of_rows, var4, separation);
   
    ex_cart_coor = cr_fitting(cart_coor, catmull_rom, num_of_rows, num_of_coils, resolution, separation);
 
    ex_current = cr_fitting(current, num_of_rows, num_of_coils, resolution, separation);
  
    ex_cyl_coor = cart_to_cyl(ex_cart_coor, num_of_rows + (num_of_rows - num_of_coils) * resolution);
  
    ex_separation = new int [num_of_coils];
    for (int i = 0; i < num_of_coils; i++) {
        ex_separation[i] = separation[i] + (separation[i] - i - 1) * resolution;     
    }

    printf("...\n");
}

coils::~coils() {

    printf("deleting an object of class <coils>...\n");
    
    for (int i = 0; i < num_of_rows; i++) {        
        delete[] cart_coor[i];
        delete[] cyl_coor[i];
        delete[] catmull_rom[i];
    }
    // for (int i = 0; i < num_of_rows + (num_of_rows - num_of_coils) * resolution; i++) {
    //     delete[] ex_cart_coor[i];
    //     delete[] ex_cyl_coor[i];
    // }
    delete[] separation;
    delete[] current; 
    delete[] cart_coor;
    delete[] cyl_coor;
    delete[] catmull_rom;
    // delete[] ex_cart_coor;
    // delete[] ex_cyl_coor;
    // delete[] ex_current;

    printf("done...\n");
}
