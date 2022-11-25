#include "read_write.hpp"

char* import_parameters (const char* var1, const char* var2) {

    printf("importing parameters...\n");

    fstream file1;
    file1.open(var1,ios::in);
    if (!file1.is_open()) {
        printf("file does not exist or has improper access permission\n");
        abort();
    }

    file1.seekg(0,ios::end);
    int num1 = file1.tellg();
    file1.seekg(0);

    int num2 = 0;
    for (int i = 0; i < num1; i++) {
        if(file1.get() == '\n') {
            num2++;
        }
    }
    file1.close();

    file1.open(var1,ios::in);

    char char_temp;
    char str0[num1 - num2];
    for (int i = 0; i < num1 - num2; i++) {
        str0[i] = file1.get();
    }
    file1.close();

    char* str_var1_0 = strstr(str0,var2);
    char* str_var1_1 = strstr(str_var1_0,"=");
    
    int num3 = 0, num4 = 0;
    for (int i = 0; ;i++) {

        if ( i > 99999) {
            printf("error in importing parameter;\n");
            abort();
        }
        if (str_var1_1[i+1] == '\n' || str_var1_1[i+1] == '\0') {
            break;
        }
        if (str_var1_1[i+1] == ' ') {
            continue;
        }
        num4++;

    }

    char* str1 = nullptr;
    str1 = new char [num4+1]; 
    for (int i = 0; ;i++) {

        if ( i > 99999) {
            printf("error in importing parameter;\n");
            abort();
        }
        if (str_var1_1[i+1] == '\n' || str_var1_1[i+1] == '\0') {
            break;
        }
        if (str_var1_1[i+1] == ' ') {
            num3++;
            continue;
        }
        str1[i-num3] = str_var1_1[i+1];

    }
    str1[num4] = '\0';

    printf("done...\n");
    printf("\n");

    return str1;

}

int num_of_row(const char *var1) { 
    
    printf("importing the data of coils...\n");
    fstream file1(var1,ios::in);
    if (!file1.is_open()) {
        printf("file does not exist or has improper access permission\n");
        abort();
    }

    file1.seekg(0, ios::end);
    int num1 = file1.tellg();
    int num2 = 1;
    char char1;
    file1.seekg(0);
    for (int i = 0; i <= num1; i++ ) {
        file1.get(char1);
             if (char1 == '\n') {
                num2++;
        }
    }
    file1.close();
    
    printf("...\n");
    return num2;   
}

double** coordinates_initial(const char* var1, int var2) {
    printf("...\n");
    fstream file1(var1,ios::in);
    if (!file1.is_open()) {
        printf("file does not exist or has improper access permission\n");
        abort();
    }

    double** array1 = nullptr;
    double* array2 = nullptr;
    array1 = new double* [var2];
    array2 = new double [3];   
    int num1 = 0, num2 = 0;
     
    file1 >> *(array2 + 0); 
    file1 >> *(array2 + 1); 
    file1 >> *(array2 + 2); 
    file1.seekg(0);
    
    for (int i = 0; i < var2; i++) {
        *(array1+i) = new double [4];   
    }

    for (int i = 0; i < var2; i++) {
        for (int j = 0; j < 4; j++) {
            file1 >> *(*(array1 + i) + j);
            if (fabs(array1[i][j]) < pow(10,-12)) {
                array1[i][j] = 0;
            }           
        }
    }

    for (int i = 1; i < var2; i++) {
        if (*(array2 + 0) == *(*(array1 + i) + 0)) {
            if (*(array2 + 1) == *(*(array1 + i) + 1)) {
                if (*(array2 + 2) == *(*(array1 + i) + 2)) { 
                    num1++;
                    if (i < var2 -1) {
                        *(array2 + 0) = *(*(array1 + i + 1) + 0);
                        *(array2 + 1) = *(*(array1 + i + 1) + 1);
                        *(array2 + 2) = *(*(array1 + i + 1) + 2);
                        i++;
                    }  
                }
            }
        }       
    }

    file1.seekg(0);
    file1 >> *(array2 + 0); 
    file1 >> *(array2 + 1); 
    file1 >> *(array2 + 2); 

    double** array3 = nullptr;
    array3 = new double* [var2 + 2];
    for (int i = 0; i < var2 + 2; i++) {
        if (i < var2) {
            *(array3+i) = new double [4];
        }
        else if (i == var2) {
            *(array3+i) = new double [num1];
        }
        else if (i == var2 + 1) {
            *(array3+i) = new double [1];
        }          
    }

    for (int i = 0; i < var2; i++) {
        for (int j = 0; j < 4; j++) {
            *(*(array3 + i) + j) = *(*(array1 + i) + j); 
        }
    }
    for (int i = 1; i < var2; i++) {
        if (*(array2 + 0) == *(*(array1 + i) + 0)) {
            if (*(array2 + 1) == *(*(array1 + i) + 1)) {
                if (*(array2 + 2) == *(*(array1 + i) + 2)) {
                    num2++; 
                    if (i < var2 - 1) {
                        *(array2 + 0) = *(*(array1 + i + 1) + 0);
                        *(array2 + 1) = *(*(array1 + i + 1) + 1);
                        *(array2 + 2) = *(*(array1 + i + 1) + 2);    
                        i++;    
                    }
                    else if (i < var2) {
                        i++;
                    }
                    *(*(array3 + var2) + num2 - 1) = i;  
                }
            }
        }       
    }
    *(*(array3 + var2 + 1) + 0) = num2;

    for (int i = 0; i < var2; i++) {
        delete [] *(array1 + i);   
    }
    delete [] array1;
    delete [] array2;
    file1.close();
    
    printf("done...\n");
    printf("\n");
    return array3;  
}

void write_coils_cyl(double** var1, int var3) {
    printf("exporting the cylindrical coordinates of the refined coils ...\n");
    if (!_mkdir("output")) {
        _mkdir("output");
    }
    char array1[99];
    sprintf(array1, ".\\output\\refined_coils_cyl.txt");
    fstream file1(array1,ios::out);
    if (!file1.is_open()) {
        printf("file does not exist or has improper access permission\n");
        abort();
    }
    for (int i = 0; i < var3; i++) {
        for (int j = 0; j < 3; j++) {
            file1 << var1[i][j] << " ";
            if (j == 2) {
                file1 << "\n";
            }
        }
    }
    file1.close();
    printf("done ...\n");
    printf("\n");
}

void write_FLT(double** var1, double* var2, int var3) {
    printf("exporting the cylindrical coordinates of the points of FLT ...\n");
    if (!_mkdir("output")) {
        _mkdir("output");
    }
    char array1[99];
    sprintf(array1, ".\\output\\FLT_R_%0.4f_zeta_%0.4f_Z_%0.4f.txt", var2[0], var2[1], var2[2]);
    fstream file1(array1,ios::out);
    if (!file1.is_open()) {
        printf("file does not exist or has improper access permission\n");
        abort();
    }
    for (int i = 0; i < var3; i++) {
        for (int j = 0; j < 3; j++) {
            file1 << var1[i][j] << " ";
            if (j == 2) {
                file1 << "\n";
            }
        }
    }
    file1.close();
    printf("done ...\n");
    printf("\n");
}

void write_flux_iota(double* var1, double* var2, int var3) {
    printf("exporting toroidal flux and iota ...\n");
    if (!_mkdir("output")) {
        _mkdir("output");
    }
    char array1[99];
    sprintf(array1, ".\\output\\flux_iotaf.txt", var3);
    fstream file1(array1,ios::out);
    if (!file1.is_open()) {
        printf("file does not exist or has improper access permission\n");
        abort();
    }
    for (int i = 0; i < var3; i++) {
        file1 << sqrt(var1[i]/var1[0]) << " " << var1[i] << " " << var2[i] << "\n";
    }
    file1.close();
    printf("done ...\n");
    printf("\n");

}

int poincare_w(double** var1, double* var2, int var3, double var4) {
    printf("poincare section; zeta = %f*pi\n",var4/pi);
    double** array1 = nullptr;
    array1 = new double* [var3];
    for (int i = 0; i < var3; i++) {
        array1[i] = new double [3];
        for (int j = 0; j < 3; j++) {
            array1[i][j] = var1[i][j];
        }
    }

    for (int i = 0; i < var3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; ; k++) {
                if (k > var3) {
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
    
    printf("exporting the points in cylindrical coordinates...\n");
    int num1 = 0;
    if (!_mkdir("output")) {
        _mkdir("output");
    }
    char array2[99];
    sprintf(array2,".\\output\\poincare_R_%0.4f_zeta_%0.4f_Z_%0.4f_torSec_%0.4f.txt",var2[0], var2[1], var2[2], var4);
    fstream file1(array2,ios::out);
    if (!file1.is_open()) {
        printf("file does not exist or has improper access permission\n");
        abort();
    }

    for (int i = 0; i < var3; i++) {
        if(fabs(array1[i][1]-var4) < 2*pi/720) {
            num1++;
            file1 << array1[i][0] << " ";
            file1 << array1[i][1] << " ";
            file1 << array1[i][2] << " ";
            file1 << "\n";
        }
    }
    file1.close(); 

    for (int i = 0; i < var3; i++) {
        delete[] array1[i];
    }
    delete[] array1;
    printf("continue...\n");
    return num1;
}

