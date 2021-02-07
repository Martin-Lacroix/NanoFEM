#include <bits/stdc++.h>
#include <algorithm>
#include "linalg.h"
#include "stdafx.h"

#ifndef MATH_H
#define MATH_H

typedef std::vector<int> ivector;
typedef std::vector<double> dvector;
typedef alglib::sparsematrix sparse;
typedef alglib::real_2d_array matrix;
typedef alglib::real_1d_array darray;
typedef alglib::integer_1d_array iarray;
typedef std::vector<std::string> svector;
typedef std::pair<std::string,double> sdpair;

// ----------------------------------------|
// Structure storing the quadrature rule   |
// ----------------------------------------|

struct quadStruct{

    dvector weight;
    std::vector<dvector> gRST;
};

// -------------------------------------------|
// Functions available in the math.cpp file   |
// -------------------------------------------|

namespace math{
    
    quadStruct legendre(int dim,int order);
    std::vector<dvector> to2D(std::vector<dvector> &nXYZ);

    // Sparse or dense matrix-matrix addition

    void add(double k1,double k2,matrix &M1,matrix &M2);
    void add(double k1,double k2,sparse &M1,sparse &M2);

    // Matrix-matrix or matrix-array product

    darray prod(double k,matrix &M1,darray &V1,int t1=0);
    matrix prod(double k,matrix &M1,matrix &M2,int t1=0,int t2=0);

    // Fills a matrix or an array with zeros

    void zero(darray &V);
    void zero(matrix &M);

    // Non-zero indices of a sparse matrix and vector cross product

    std::vector<ivector> sparsemap(sparse &K);
    darray cross(darray &V1,darray &V2);

    // Kernel function and stiffness tensor in Voigh notation

    double kernel(dvector x,dvector y);
    matrix stiffness(double E,double v);

    // Symmetric sparse matrix operations only

    double symget(sparse &M,int row,int col);
    void symset(sparse &M,int row,int col,double val);
    void symadd(sparse &M,int row,int col,double val);
}

#endif