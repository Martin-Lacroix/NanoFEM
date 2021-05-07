#include <algorithm>
#include "linalg.h"
#include "stdafx.h"
#include <array>

#ifndef MATH_H
#define MATH_H

// Definition of usual variable types

typedef std::vector<int> ivector;
typedef std::vector<double> dvector;
typedef alglib::sparsematrix sparse;
typedef alglib::real_2d_array matrix;
typedef alglib::real_1d_array darray;
typedef std::array<double,3> array3d;

// ----------------------------------------|
// Structure storing the quadrature rule   |
// ----------------------------------------|

struct quadStruct{

    int gLen;
    dvector weight;
    std::vector<dvector> gRST;
};

// ---------------------------------------|
// Functions available in the math file   |
// ---------------------------------------|

namespace math{

    matrix stiffness(array3d LmX);
    quadStruct legendre(int dim,int order);
    dvector lagrange(int var,dvector node,dvector val);

    // Matrix-matrix or matrix-array operations

    void add(double k1,double k2,matrix &M1,matrix &M2);
    darray prod(double k,matrix &M1,darray &V1,int t1=0);
    matrix prod(double k,matrix &M1,matrix &M2,int t1=0,int t2=0);
    
    // Fills a matrix or an array with zeros

    void zero(darray &V);
    void zero(matrix &M);

    // Computes the L2 norm of a vector

    double norm(darray &V);
    double norm(array3d &V);

    // Other general fiunctions for sparse or dense matrices
    
    matrix projection(array3d norm);
    matrix invert(matrix &M,double det);
    array3d cross(array3d &V1,array3d &V2);
    array3d dotsub(array3d &V1,array3d &V2);
    std::vector<ivector> sparsemap(sparse &M);

    // Symmetric sparse matrix operations only

    double get(sparse &M,int row,int col);
    void symset(sparse &M,int row,int col,double val);
    void symadd(sparse &M,int row,int col,double val);
}

#endif