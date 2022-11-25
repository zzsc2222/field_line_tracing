#pragma once
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <ctime>
#include <omp.h>
#include <eigen-3.4.0/Eigen/Dense>
#include <direct.h>
// #include <iostream>
// #include <memory>

// #define num_thread 5
#define pi  3.14159265358979
#define mu0 0.000001256637061435917
#define bs mu0/(4*pi)
#define bz pow(1.380649, -23) 
#define ev pow(1.6021766208,-19)
#define me pow(9.10956, -31)
#define mi pow(1.6726231, -27)


// using std::cout;
// using std::cin;
// using std::endl;
using std::stod;
using std::stoi;
using std::fstream;
using std::ios;
using std::max;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::ComputeThinU;
using Eigen::ComputeThinV;
// using std::unique_ptr;
// using std::make_unique;
