#include "mag_field.hpp"

double* biot_savart_cart(double var1_0, double var1_1, double var1_2, double** var2, double* var3, int var4, int* var5, int var6) { 

    double* array1 = nullptr;
    array1 = new double [3];

    double sum0, sum1, sum2, array2[var4-1];
    double dx[var4-1], dy[var4-1], dz[var4-1];
    double xx[var4-1], yy[var4-1], zz[var4-1];

    for (int i = 0; i < var6; i++) {
        var3[var5[i]-1] = 0;
    }

    #ifdef num_thread
    #pragma omp parallel for reduction(+:sum0) reduction(+:sum1) reduction(+:sum2) num_threads(num_thread)
    #else
    #pragma omp parallel for reduction(+:sum0) reduction(+:sum1) reduction(+:sum2)
    #endif
    for (int i = 0; i < var4 - 1; i++) {    
       
        dx[i] = var2[i+1][0]-var2[i][0];
        dy[i] = var2[i+1][1]-var2[i][1];
        dz[i] = var2[i+1][2]-var2[i][2];

        xx[i] = var1_0-var2[i][0];
        yy[i] = var1_1-var2[i][1];
        zz[i] = var1_2-var2[i][2];

        array2[i] = bs * (var3[i]*pow(10,6)) / pow( pow(xx[i],2) + pow(yy[i],2) + pow(zz[i],2), 1.5);

        sum0 += ( (zz[i]*dy[i] - yy[i]*dz[i]) * array2[i] );
        sum1 += ( (xx[i]*dz[i] - zz[i]*dx[i]) * array2[i] );
        sum2 += ( (yy[i]*dx[i] - xx[i]*dy[i]) * array2[i] );

    }

    array1[0] = sum0;
    array1[1] = sum1;
    array1[2] = sum2;

    return array1;  
}

double* biot_savart_cyl(double var1_0, double var1_1, double var1_2, double** var2, double* var3, int var4, int* var5, int var6) {

    extern double* cyl_to_cart(double*);

    double* array1 = nullptr; 
    double* array2 = nullptr; 
    double* array3 = nullptr; 
    double* array4 = nullptr;
    array1 = new double [3];
    array2 = new double [3];
    array3 = new double [3];
    array4 = new double [3];

    array1[0] = var1_0;
    array1[1] = var1_1;
    array1[2] = var1_2;
    array2 = cyl_to_cart(array1);
    array3 = biot_savart_cart(array2[0], array2[1], array2[2], var2, var3, var4, var5, var6);

    array4[0] =  array3[0] * cos(var1_1) + array3[1] * sin(var1_1);
    array4[1] = -array3[0] * sin(var1_1) + array3[1] * cos(var1_1);
    array4[2] =  array3[2];
    
    delete[] array1;
    delete[] array2;
    delete[] array3;

    return array4;
}

double** adams_bashforth_moulton_cyl(double* var1, double var2, int var3, double** var4, double* var5, int var6, int* var7, int var8) {

    printf("field line tracing...\n");
    printf("starting from R/zeta/Z = %0.6f/%0.6f/%0.6f\n", var1[0], var1[1], var1[2]);

    time_t time;
    
    double** array1 = nullptr;
    double** array2 = nullptr;
    double** array3 = nullptr;
    array1 = new double* [var3];
    array2 = new double* [4];
    array3 = new double* [var3];
    for (int i = 0; i < var3; i++) {
        array1[i] = new double [3];
    }
    for (int i = 0; i < 4; i++) {
        array2[i] = new double [3];
    }
    
    double num1_r, num2_r, num3_r, num4_r, num5_r, num6_r, num7_r;
    double num1_t, num2_t, num3_t, num4_t, num5_t, num6_t, num7_t;
    double num1_z, num2_z, num3_z, num4_z, num5_z, num6_z, num7_z;
    double num1;
    double num1_rr, num2_rr, num3_rr, num4_rr, num5_rr, num6_rr;
    double num1_tt, num2_tt, num3_tt, num4_tt, num5_tt, num6_tt;
    double num1_zz, num2_zz, num3_zz, num4_zz, num5_zz, num6_zz;

    array1[0][0] = var1[0];
    array1[0][1] = var1[1];
    array1[0][2] = var1[2];
    array3[0] = biot_savart_cyl(var1[0], var1[1], var1[2], var4, var5, var6, var7, var8);
    
    int num2 = 1;
    for (int i = 1; i < var3; i++) {
        time = clock();

        if (i == num2) {
            printf("tracing... i = %d; time = %d ms\n",i,(int) time);
            num2 += ceil(var3/10);
        }
        else if (i == var3 - 1) {
            printf("tracing... i = %d; time = %d ms\n",i,(int) time);
        }

        if (i < 5) {
            array2[0] = biot_savart_cyl(array1[i-1][0], \
            array1[i-1][1], \
            array1[i-1][2], var4, var5, var6, var7, var8);
            num1 = sqrt( pow(array2[0][0],2) + pow(array2[0][1],2) + pow(array2[0][2],2) );
            num1_r = var2 * array2[0][0] / num1;
            num1_t = var2 * array2[0][1] / num1 / (array1[i-1][0]);
            num1_z = var2 * array2[0][2] / num1;

            array2[1] = biot_savart_cyl(array1[i-1][0] + (num1_r/2), \
            array1[i-1][1] + (num1_t/2), \
            array1[i-1][2] + (num1_z/2), var4, var5, var6, var7, var8);
            num1 = sqrt( pow(array2[1][0],2) + pow(array2[1][1],2) + pow(array2[1][2],2) );
            num2_r = var2 * array2[1][0] / num1;
            num2_t = var2 * array2[1][1] / num1 / (array1[i-1][0] + (num1_r/2));
            num2_z = var2 * array2[1][2] / num1;

            array2[2] = biot_savart_cyl(array1[i-1][0] + (num1_r*2/9) + (num2_r*4/9), \
            array1[i-1][1] + (num1_t*2/9) + (num2_t*4/9), \
            array1[i-1][2] + (num1_z*2/9) + (num2_z*4/9), var4, var5, var6, var7, var8);
            num1 = sqrt( pow(array2[2][0],2) + pow(array2[2][1],2) + pow(array2[2][2],2) );
            num3_r = var2 * array2[2][0] / num1;
            num3_t = var2 * array2[2][1] / num1 / (array1[i-1][0] + (num1_r*2/9) + (num2_r*4/9));
            num3_z = var2 * array2[2][2] / num1;

            array2[3] = biot_savart_cyl(array1[i-1][0] + (num1_r*7/36) + (num2_r*2/9) - (num3_r/12), \
            array1[i-1][1] + (num1_t*7/36) + (num2_t*2/9) - (num3_t/12), \
            array1[i-1][2] + (num1_z*7/36) + (num2_z*2/9) - (num3_z/12), var4, var5, var6, var7, var8);
            num1 = sqrt( pow(array2[3][0],2) + pow(array2[3][1],2) + pow(array2[3][2],2) );
            num4_r = var2 * array2[3][0] / num1;
            num4_t = var2 * array2[3][1] / num1 / (array1[i-1][0] + (num1_r*7/36) + (num2_r*2/9) - (num3_r/12));
            num4_z = var2 * array2[3][2] / num1;

            array2[3] = biot_savart_cyl(array1[i-1][0] - (num1_r*35/144) - (num2_r*55/36) + (num3_r*35/48) + (num4_r*15/8), \
            array1[i-1][1] - (num1_t*35/144) - (num2_t*55/36) + (num3_t*35/48) + (num4_t*15/8), \
            array1[i-1][2] - (num1_z*35/144) - (num2_z*55/36) + (num3_z*35/48) + (num4_z*15/8), var4, var5, var6, var7, var8);
            num1 = sqrt( pow(array2[3][0],2) + pow(array2[3][1],2) + pow(array2[3][2],2) );
            num5_r = var2 * array2[3][0] / num1;
            num5_t = var2 * array2[3][1] / num1 / (array1[i-1][0] - (num1_r*35/144) - (num2_r*55/36) + (num3_r*35/48) + (num4_r*15/8));
            num5_z = var2 * array2[3][2] / num1;

            array2[3] = biot_savart_cyl(array1[i-1][0] - (num1_r/360) - (num2_r*11/36) - (num3_r/8) + (num4_r/2) + (num5_r/10), \
            array1[i-1][1] - (num1_t/360) - (num2_t*11/36) - (num3_t/8) + (num4_t/2) + (num5_t/10), \
            array1[i-1][2] - (num1_z/360) - (num2_z*11/36) - (num3_z/8) + (num4_z/2) + (num5_z/10), var4, var5, var6, var7, var8);
            num1 = sqrt( pow(array2[3][0],2) + pow(array2[3][1],2) + pow(array2[3][2],2) );
            num6_r = var2 * array2[3][0] / num1;
            num6_t = var2 * array2[3][1] / num1 / (array1[i-1][0] - (num1_r/360) - (num2_r*11/36) - (num3_r/8) + (num4_r/2) + (num5_r/10));
            num6_z = var2 * array2[3][2] / num1;

            array2[3] = biot_savart_cyl(array1[i-1][0] - (num1_r*41/260) + (num2_r*23/13) + (num3_r*43/156) - (num4_r*118/39) + (num5_r*32/195) + (num6_r*80/39), \
            array1[i-1][1] - (num1_t*41/260) + (num2_t*23/13) + (num3_t*43/156) - (num4_t*118/39) + (num5_t*32/195) + (num6_t*80/39), \
            array1[i-1][2] - (num1_z*41/260) + (num2_z*23/13) + (num3_z*43/156) - (num4_z*118/39) + (num5_z*32/195) + (num6_z*80/39), var4, var5, var6, var7, var8);
            num1 = sqrt( pow(array2[3][0],2) + pow(array2[3][1],2) + pow(array2[3][2],2) );
            num7_r = var2 * array2[3][0] / num1;
            num7_t = var2 * array2[3][1] / num1 / (array1[i-1][0] - (num1_r*41/260) + (num2_r*23/13) + (num3_r*43/156) - (num4_r*118/39) + (num5_r*32/195) + (num6_r*80/39));
            num7_z = var2 * array2[3][2] / num1;


            array1[i][0] = array1[i-1][0] + (num1_r*13/200 + num3_r*11/40 + num4_r*11/40 + num5_r*4/25 + num6_r*4/25 + num7_r*13/200) / 6;
            array1[i][1] = array1[i-1][1] + (num1_t*13/200 + num3_t*11/40 + num4_t*11/40 + num5_t*4/25 + num6_t*4/25 + num7_t*13/200) / 6;
            array1[i][2] = array1[i-1][2] + (num1_z*13/200 + num3_z*11/40 + num4_z*11/40 + num5_z*4/25 + num6_z*4/25 + num7_z*13/200) / 6; 
            array3[i] = biot_savart_cyl(array1[i][0], array1[i][1], array1[i][2], var4, var5, var6, var7, var8);
        }
        else {
            num1 = sqrt( pow(array3[i-5][0],2) + pow(array3[i-5][1],2) + pow(array3[i-5][2],2) );
            num1_rr = array3[i-5][0] / num1;
            num1_tt = array3[i-5][1] / num1 / (array1[i-5][0]);
            num1_zz = array3[i-5][2] / num1;

            num1 = sqrt( pow(array3[i-4][0],2) + pow(array3[i-4][1],2) + pow(array3[i-4][2],2) );
            num2_rr = array3[i-4][0] / num1;
            num2_tt = array3[i-4][1] / num1 / (array1[i-4][0]);
            num2_zz = array3[i-4][2] / num1;
            
            num1 = sqrt( pow(array3[i-3][0],2) + pow(array3[i-3][1],2) + pow(array3[i-3][2],2) );
            num3_rr = array3[i-3][0] / num1;
            num3_tt = array3[i-3][1] / num1 / (array1[i-3][0]);
            num3_zz = array3[i-3][2] / num1;
            
            num1 = sqrt( pow(array3[i-2][0],2) + pow(array3[i-2][1],2) + pow(array3[i-2][2],2) );
            num4_rr = array3[i-2][0] / num1;
            num4_tt = array3[i-2][1] / num1 / (array1[i-2][0]);
            num4_zz = array3[i-2][2] / num1;
 
            num1 = sqrt( pow(array3[i-1][0],2) + pow(array3[i-1][1],2) + pow(array3[i-1][2],2) );
            num5_rr = array3[i-1][0] / num1;
            num5_tt = array3[i-1][1] / num1 / (array1[i-1][0]);
            num5_zz = array3[i-1][2] / num1;

            array1[i][0] = array1[i-1][0] + var2 * ((1901.0/720.0)*num5_rr + (-2774.0/720.0)*num4_rr + (2616.0/720.0)*num3_rr + (-1274.0/720.0)*num2_rr + (251.0/720.0)*num1_rr);
            array1[i][1] = array1[i-1][1] + var2 * ((1901.0/720.0)*num5_tt + (-2774.0/720.0)*num4_tt + (2616.0/720.0)*num3_tt + (-1274.0/720.0)*num2_tt + (251.0/720.0)*num1_tt);
            array1[i][2] = array1[i-1][2] + var2 * ((1901.0/720.0)*num5_zz + (-2774.0/720.0)*num4_zz + (2616.0/720.0)*num3_zz + (-1274.0/720.0)*num2_zz + (251.0/720.0)*num1_zz);

            array3[i] = biot_savart_cyl(array1[i][0], array1[i][1], array1[i][2], var4, var5, var6, var7, var8);
            num1 = sqrt( pow(array3[i][0],2) + pow(array3[i][1],2) + pow(array3[i][2],2) );
            num6_rr = array3[i][0] / num1;
            num6_tt = array3[i][1] / num1 / (array1[i][0]);
            num6_zz = array3[i][2] / num1;
            
            array1[i][0] = array1[i-1][0] + var2 * ((251.0/720.0)*num6_rr + (646.0/720.0)*num5_rr + (-264.0/720.0)*num4_rr + (106.0/720.0)*num3_rr + (-19.0/720.0)*num2_rr);
            array1[i][1] = array1[i-1][1] + var2 * ((251.0/720.0)*num6_tt + (646.0/720.0)*num5_tt + (-264.0/720.0)*num4_tt + (106.0/720.0)*num3_tt + (-19.0/720.0)*num2_tt);
            array1[i][2] = array1[i-1][2] + var2 * ((251.0/720.0)*num6_zz + (646.0/720.0)*num5_zz + (-264.0/720.0)*num4_zz + (106.0/720.0)*num3_zz + (-19.0/720.0)*num2_zz);  
   
            array3[i] = biot_savart_cyl(array1[i][0], array1[i][1], array1[i][2], var4, var5, var6, var7, var8);
        } 
    }
    
    delete[] array2;
    delete[] array3;

    printf("done... \n");
    printf("\n");

    return array1;
}

double** poincare_fitting(double** var1, int var2, int var3, double** var4) {

    printf("2D curve fitting prepared for flux calculation...\n");

    Matrix<double, Dynamic, Dynamic> matrix1(var2, 2*var3+1);
    Matrix<double, Dynamic, 1> matrix2_r(var2,1);
    Matrix<double, Dynamic, 1> matrix2_z(var2,1);
    Matrix<double, Dynamic, 1> matrix3_r(2*var3+1,1);
    Matrix<double, Dynamic, 1> matrix3_z(2*var3+1,1);

    double num1 = 0, num2 = 0;

    for (int i = 1; i < 2 * var3; i += 2) {
        num1 += var4[0][i] * cos( (i+1)/2 * var1[0][1]) + var4[0][i+1] * sin((i+1)/2 * var1[0][1] );
        num2 += var4[1][i] * cos( (i+1)/2 * var1[0][1]) + var4[1][i+1] * sin((i+1)/2 * var1[0][1] );
    }
    num1 += var4[0][0];
    num2 += var4[1][0];

    // printf("R0 = %f\n",num1);
    // printf("Z0 = %f\n",num2);

    for (int i = 0; i < var2; i++) {
        for (int j = 1; j < 2 * var3; j += 2) {
            matrix1(i,j  ) = cos( (j+1)/2 * atan2(var1[i][2]-num2, var1[i][0]-num1) );
            matrix1(i,j+1) = sin( (j+1)/2 * atan2(var1[i][2]-num2, var1[i][0]-num1) );
        }
        matrix1(i,0) = 1;

        matrix2_r(i,0) = var1[i][0];
        matrix2_z(i,0) = var1[i][2];
    }
    matrix3_r = matrix1.bdcSvd(ComputeThinU | ComputeThinV).solve(matrix2_r);
    matrix3_z = matrix1.bdcSvd(ComputeThinU | ComputeThinV).solve(matrix2_z);

    // matrix3_r = matrix1.completeOrthogonalDecomposition().solve(matrix2_r);
    // matrix3_z = matrix1.completeOrthogonalDecomposition().solve(matrix2_z);

    // cout << "R = " << matrix3_r << endl;
    // cout << "Z = " << matrix3_z << endl;

    double** array1 = nullptr;
    array1 = new double* [2];
    array1[0] = new double[2*var3+1];
    array1[1] = new double[2*var3+1];
    for (int i = 0; i < 2*var3+1;i++) {
        array1[0][i] = matrix3_r(i,0);
        array1[1][i] = matrix3_z(i,0);
    }
    
    return array1;    
}

double** poincare(double** var1, int var2, double var3, int var4) {

    printf("...\n");

    double** array1 = nullptr;
    array1 = new double* [var2];
    for (int i = 0; i < var2; i++) {
        array1[i] = new double [3];
        for (int j = 0; j < 3; j++) {
            array1[i][j] = var1[i][j];
        }
    }
    
    for (int i = 0; i < var2; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; ; k++) {
                if (k > var2) {
                    printf("loop error in poincare plot\n");
                    abort();
                }
                if (array1[i][1] > 2*pi) {
                    array1[i][1] -= 2*pi;
                }
                else {
                    break;
                }
            }  
        }
    }

    double** array2 = nullptr;
    array2 = new double* [var4];
    for (int i = 0; i < var4; i++) {
        array2[i] = new double [3];
    }

    int num1 = 0;
    for (int i = 0; i < var2; i++) {
        if(fabs(array1[i][1]-var3) < 2*pi/720) {
            for (int j = 0; j < 3; j++) {
                array2[num1][j] = array1[i][j];   
            }
            num1++;
        }
    }

    for (int i = 0; i < var2; i++) {
        delete[] array1[i];
    }
    delete[] array1;

    printf("done...\n");
    printf("\n");

    return array2;
}

double** mag_axis_cal(double var1, double var2, double var3, double var4, int var5, int var6,    double** var7, double* var8, int var9, int* var10, int var11) {
    
    printf("magnetic axis calculating...\n");

    time_t time;
    FILE saved_stdout = *stdout;

    double** array1 = nullptr;
    array1 = new double* [var5];
    for (int i = 0; i < var5; i++) {
        array1[i] = new double [3];
    }
    double array2[3];
    array2[0] = var1;
    array2[1] = var2;
    array2[2] = var3;

    int num1;
    double num2 = -var4, num3;
    double max_r, min_r;
    double max_z, min_z;
    
    time = clock();

    for (int iter = 0; iter < 10; iter++) {

        printf("iteration = %d/10\n",iter+1);
        *stdout = *fopen("log.txt","w");
        num1 = 0;

        for (int i = 0; ; i++) {

            array2[0] += num2;
            max_r = INT_MIN;
            min_r = INT_MAX;
            array1 = adams_bashforth_moulton_cyl(array2, var4, var5, var7, var8, var9 ,var10 ,var11);

            for (int j = 0; j < var5; j++) {
                for (int k = 0; ; k++) {
                    if ( k > 99999) {
                        abort();
                    }
                    if (array1[j][1] > 2*pi) {
                        array1[j][1] -= 2*pi;
                    }
                    else {
                        break;
                    }
                }
            }

            for (int j = 0; j < var5; j++) {
                if ( fabs(array1[j][1]-var2) < 2*pi/720 ) {
                    if ( array1[j][0] > max_r ) {
                        max_r = array1[j][0];
                    }
                    if (array1[j][0] < min_r) {
                        min_r = array1[j][0];
                    }
                }
            }
        
            if (i == 0) {
                num3 = max_r - min_r;
            }
            else {
                if ( (max_r - min_r) < num3 ) {
                    num3 = max_r - min_r;
                }
                else {
                    num2 = -num2;
                    num1++;
                }
            }

            if (num1 == 2) {
                array2[0] -= num2/2;
                break;
            }
    
        }
        num2 = -fabs(num2);

        num1 = 0;
        for (int i = 0; ; i++) {

            array2[2] += num2;
            max_z = INT_MIN;
            min_z = INT_MAX;
            array1 = adams_bashforth_moulton_cyl(array2, var4, var5, var7, var8, var9 ,var10 ,var11);

            for (int j = 0; j < var5; j++) {
                for (int k = 0; ; k++) {
                    if ( k > 99999) {
                        abort();
                    }
                    if (array1[j][1] > 2*pi) {
                        array1[j][1] -= 2*pi;
                    }
                    else {
                        break;
                    }
                }
            }

            for (int j = 0; j < var5; j++) {
                if ( fabs(array1[j][1]-var2) < 2*pi/720 ) {
                    if ( array1[j][2] > max_z ) {
                        max_z = array1[j][2];
                    }
                    if (array1[j][2] < min_z) {
                        min_z = array1[j][2];
                    }
                }
            }
        
            if (i == 0) {
                num3 = max_z - min_z;
            }
            else {
                if ( (max_z - min_z) < num3 ) {
                    num3 = max_z - min_z;
                }
                else {
                    num2 = -num2;
                    num1++;
                }
            }
    
            if (num1 == 2) {
                array2[2] -= num2/2;
                break;
            }
        }
        num2 = -fabs(num2);

        *stdout = saved_stdout;
        if (iter > 0 && iter%2 == 0) {
            num2 /= 2;
        }
        
        time = clock();

        printf(" time = %d ms\n",(int) time);
    }
    
    printf("R/zetz/Z = %f/%f/%f\n",array2[0],array2[1],array2[2]);
    
    *stdout = *fopen("log.txt","w");
    array1 = adams_bashforth_moulton_cyl(array2, var4, var5, var7, var8, var9 ,var10 ,var11);
    *stdout = saved_stdout;

    Matrix<double, Dynamic, Dynamic> matrix1(var5, 2*var6+1);
    Matrix<double, Dynamic, 1> matrix2_r(var5,1);
    Matrix<double, Dynamic, 1> matrix2_z(var5,1);
    Matrix<double, Dynamic, 1> matrix3_r(2*var6+1,1);
    Matrix<double, Dynamic, 1> matrix3_z(2*var6+1,1);

    for (int i = 0; i < var5; i++) {
        for (int j = 1; j < 2 * var6; j += 2) {
            matrix1(i,j) = cos( (j+1)/2.0 * array1[i][1] );
            matrix1(i,j+1) = sin( (j+1)/2.0 * array1[i][1] );
        }
        matrix1(i,0) = 1;

        matrix2_r(i,0) = array1[i][0];
        matrix2_z(i,0) = array1[i][2];
    }
 
    matrix3_r = matrix1.bdcSvd(ComputeThinU | ComputeThinV).solve(matrix2_r);
    matrix3_z = matrix1.bdcSvd(ComputeThinU | ComputeThinV).solve(matrix2_z);
    // matrix3_r = matrix1.jacobiSvd(ComputeThinU | ComputeThinV).solve(matrix2_r);
    // matrix3_z = matrix1.jacobiSvd(ComputeThinU | ComputeThinV).solve(matrix2_z);
    // matrix3_r = matrix1.completeOrthogonalDecomposition().solve(matrix2_r);
    // matrix3_z = matrix1.completeOrthogonalDecomposition().solve(matrix2_z);
    // matrix3_r = matrix1.fullPivHouseholderQr().solve(matrix2_r);
    // matrix3_z = matrix1.fullPivHouseholderQr().solve(matrix2_z);
    // matrix3_r = matrix1.colPivHouseholderQr().solve(matrix2_r);
    // matrix3_z = matrix1.colPivHouseholderQr().solve(matrix2_z);
    // matrix3_r = matrix1.householderQr().solve(matrix2_r);
    // matrix3_z = matrix1.householderQr().solve(matrix2_z);
    // matrix3_r = (matrix1.transpose() * matrix1).ldlt().solve(matrix1.transpose() * matrix2_r);
    // matrix3_z = (matrix1.transpose() * matrix1).ldlt().solve(matrix1.transpose() * matrix2_z);

    double** array3 = nullptr;
    array3 = new double* [2];
    array3[0] = new double[2*var6+1];
    array3[1] = new double[2*var6+1];
    for (int i = 0; i < 2*var6+1;i++) {
        array3[0][i] = matrix3_r(i,0);
        array3[1][i] = matrix3_z(i,0);
    }
    // cout << "R = " << matrix3_r << endl;
    // cout << "Z = " << matrix3_z << endl;

    for (int i = 0; i < var5; i++) {
        delete[] array1[i];
    }
    delete[] array1;

    printf("done...\n");
    printf("...\n");

    return array3;
}

double flux_tor_cal(double** var1, int var2, double var3, double** var4, double* var5, int var6, int* var7, int var8) {

    printf("flux calculating...\n");

    time_t time;
    time = clock();

    double num1 = 0, num2 = 0;
    int row = 360, col = 360;
    printf("rho x theta = %d x %d\n", row, col);

    double** array1 = nullptr;
    double** array2 = nullptr; 
    double*** array3 = nullptr;
    double** r = nullptr; 
    double** z = nullptr; 
    array1 = new double* [row];
    array2 = new double* [row];
    array3 = new double** [row];
    r = new double* [row];
    z = new double* [row];
    for (int i = 0; i < row; i++) {
        array1[i] = new double [col];
        array2[i] = new double [col];
        array3[i] = new double* [col];
        for (int j = 0; j < col; j++) {
            array3[i][j] = new double [5];

        }
   }
   for (int i = 0; i < row; i++) {
        r[i] = new double [col];
        z[i] = new double [col];
   }

    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            for (int k = 1; k < 2 * var2 + 1; k += 2) {
                r[i][j] += ((double)row-i)/row * (cos((k+1)/2*(2*pi/col*j)) * var1[0][k] + sin((k+1)/2*(2*pi/col*j)) * var1[0][k+1]); 
                z[i][j] += ((double)row-i)/row * (cos((k+1)/2*(2*pi/col*j)) * var1[1][k] + sin((k+1)/2*(2*pi/col*j)) * var1[1][k+1]); 
                if (k == 2 * var2 - 1) {
                    r[i][j] += var1[0][0];
                    z[i][j] += var1[1][0];
                }
            }                  
        }
    }

    #ifdef num_thread
    #pragma omp parallel for num_threads(num_thread)
    #else
    #pragma omp parallel for 
    #endif
    for (int i = 0; i < row - 1; i++) {
        for (int j = 0; j < col; j++) {
            if (j < col -1) {
                array3[i][j][0] = sqrt( pow(r[i  ][j+1]-r[i  ][j  ],2) + pow(z[i  ][j+1]-z[i  ][j  ],2) );
                array3[i][j][1] = sqrt( pow(r[i+1][j  ]-r[i  ][j+1],2) + pow(z[i+1][j  ]-z[i  ][j+1],2) );
                array3[i][j][2] = sqrt( pow(r[i+1][j  ]-r[i  ][j  ],2) + pow(z[i+1][j  ]-z[i  ][j  ],2) );
                array3[i][j][3] = sqrt( pow(r[i+1][j+1]-r[i+1][j  ],2) + pow(z[i+1][j+1]-z[i+1][j  ],2) );
                array3[i][j][4] = sqrt( pow(r[i  ][j+1]-r[i+1][j+1],2) + pow(z[i  ][j+1]-z[i+1][j+1],2) );
            }
            else {
                array3[i][j][0] = sqrt( pow(r[i  ][0  ]-r[i  ][j  ],2) + pow(z[i  ][0  ]-z[i  ][j  ],2) );
                array3[i][j][1] = sqrt( pow(r[i+1][j  ]-r[i  ][0  ],2) + pow(z[i+1][j  ]-z[i  ][0  ],2) );
                array3[i][j][2] = sqrt( pow(r[i+1][j  ]-r[i  ][j  ],2) + pow(z[i+1][j  ]-z[i  ][j  ],2) );
                array3[i][j][3] = sqrt( pow(r[i+1][0  ]-r[i+1][j  ],2) + pow(z[i+1][0  ]-z[i+1][j  ],2) );
                array3[i][j][4] = sqrt( pow(r[i  ][0  ]-r[i+1][0  ],2) + pow(z[i  ][0  ]-z[i+1][0  ],2) );
            }
            array1[i][j] = sqrt( \
            (fabs(+array3[i][j][0]+array3[i][j][1]+array3[i][j][2])) * \
            (fabs(-array3[i][j][0]+array3[i][j][1]+array3[i][j][2])) * \
            (fabs(+array3[i][j][0]-array3[i][j][1]+array3[i][j][2])) * \
            (fabs(+array3[i][j][0]+array3[i][j][1]-array3[i][j][2])) ) / 4 + sqrt( \
            (fabs(+array3[i][j][1]+array3[i][j][3]+array3[i][j][4])) * \
            (fabs(-array3[i][j][1]+array3[i][j][3]+array3[i][j][4])) * \
            (fabs(+array3[i][j][1]-array3[i][j][3]+array3[i][j][4])) * \
            (fabs(+array3[i][j][1]+array3[i][j][3]-array3[i][j][4])) ) / 4;

            if (j < col -1) {
                array2[i][j] = *(biot_savart_cyl((r[i][j]+r[i+1][j]+r[i][j+1]+r[i+1][j+1])/4, var3, (z[i][j]+z[i+1][j]+z[i][j+1]+z[i+1][j+1])/4, var4, var5, var6, var7, var8)+1);
            }
            else {
                array2[i][j] = *(biot_savart_cyl((r[i][j]+r[i+1][j]+r[i][0]+r[i+1][0])/4, var3, (z[i][j]+z[i+1][j]+z[i][0]+z[i+1][0])/4, var4, var5, var6, var7, var8)+1);
            }
        }
    }

    for (int i = 0; i < row - 1; i++) {
        for (int j = 0; j < col; j++) {
            num2 += array1[i][j];        
            num1 += array1[i][j] * array2[i][j];
        }
    }
    num1 += (*(biot_savart_cyl(var1[0][0], var3, var1[1][0], var4, var5, var6, var7, var8)+1)) * num2 / (pow(row,2)-1);
 
    for (int i = 0; i < row; i++) {
        delete[] array1[i];
        delete[] array2[i];
        delete[] r[i];
        delete[] z[i];     
        for (int j = 0; j < col; j++) {

           delete[] array3[i][j];
        }
        delete[] array3[i];
   }
   delete[] array1;
   delete[] array2;
   delete[] r;
   delete[] z;
   delete[] array3;

  printf("PSIzeta = %f Wb\n",num1);
  printf("done...  time = %d ms\n",(int) (time));
  printf("\n");

  return num1;
}

double iota_cal(double** var1, int var2, double** var3, int var4) {

    printf("iota calculating...\n");

    double num1 = 0;
    double* array1 = nullptr;
    double* array2 = nullptr;
    double* array3 = nullptr;
    array1 = new double [var2];
    array2 = new double [var2];
    array3 = new double [var2];
    for (int i = 0; i < var2; i++) {
        array1[i] = 0;
        array2[i] = 0;
        array3[i] = 0;
    }
    
    for (int i = 0; i < var2; i++) {
        for (int j = 1; j < 2 * var4 + 1; j += 2) {
            array2[i] += (var3[0][j] * cos((j+1)/2.0 * var1[i][1]) + var3[0][j+1] * sin((j+1)/2.0 * var1[i][1]));
            array3[i] += (var3[1][j] * cos((j+1)/2.0 * var1[i][1]) + var3[1][j+1] * sin((j+1)/2.0 * var1[i][1]));
        }
        array2[i] += var3[0][0];
        array2[i] += var3[1][0];
    }

    for (int i = 0; i < var2 - 1; i++) {
 
        array1[i] = atan2(var1[i+1][2] - array3[i+1], var1[i+1][0] - array2[i+1]) - atan2(var1[i][2] - array3[i], var1[i][0] - array2[i]);
        if (array1[i] > pi) {
            num1 += array1[i] - 2*pi;
        }
        else if (array1[i] < -pi) {
            num1 += array1[i] + 2*pi;
        }
        else {
            num1 += array1[i];
        }

    }
    num1 /= (var1[var2-1][1]-var1[0][1]);
    
    printf("iota = %f\n", num1);
    printf("done... \n");
    printf("\n");

    return num1;
}

mag_field::mag_field(double** var1, double** var2, double* var3, int var4, int* var5, int var6) {

    printf("...\n");
    ex_cart_coor = var1;
    ex_cyl_coor = var2;
    ex_current = var3;
    num_of_rows = var4;
    ex_separation = var5;
    
    coef_of_axis = new double* [2];
    coef_of_axis[0] = new double [2*var6+1];
    coef_of_axis[1] = new double [2*var6+1];

    printf("done...\n");
    printf("\n");
}

mag_field::~mag_field() {

    printf("deleting an object of class <mag_field>...\n");
    for (int i = 0; i < num_of_rows; i++) {
        delete[] ex_cart_coor[i];
        delete[] ex_cyl_coor[i];
        
    }
    delete[] ex_cart_coor;
    delete[] ex_cyl_coor;
    delete[] ex_current;

    printf("done...\n");
}
