#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <array>

#ifndef MATH_H
#define MATH_H

// Definition of usual variable types

typedef std::vector<int> ivector;
typedef std::vector<double> dvector;
typedef std::array<double,3> array3d;
typedef Eigen::SparseMatrix<double> sparse;

typedef Eigen::Matrix3d matrix3d;
typedef Eigen::Vector3d darray3d;
typedef Eigen::MatrixXd matrix;
typedef Eigen::VectorXd darray;

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

    // Creates and fill a quadrature rule structure

    quadStruct makeQuad2D(darray &V, darray &val);
    quadStruct makeQuad3D(darray &V, darray &val);

    // Computes the L2 norm of a vector

    double norm(array3d &V);

    // Other general functions for sparse or dense matrices
    
    matrix projection(array3d norm);
    matrix3d invert(matrix3d &M,double det);
    array3d cross(array3d &V1,array3d &V2);
    array3d dotsub(array3d &V1,array3d &V2);
    std::vector<ivector> sparsemap(sparse &M);

    // Symmetric sparse matrix operations only

    double get(sparse &M,int row,int col);
    void symset(sparse &M,int row,int col,double val);
    void symadd(sparse &M,int row,int col,double val);
}

#endif