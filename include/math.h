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

struct quadStruct{

    dvector weight;
    std::vector<dvector> gRST;
};

namespace math{

    void zero(darray &V);
    void zero(matrix &M);
    double kernel(dvector x,dvector y);
    matrix stiffness(double E,double v);
    darray cross(darray &V1,darray &V2);
    quadStruct legendre(int dim,int order);
    void add(double k1,double k2,matrix &M1,matrix &M2);
    void add(double k1,double k2,sparse &M1,sparse &M2);
    darray prod(double k,matrix &M1,darray &V1,int t1=0);
    matrix prod(double k,matrix &M1,matrix &M2,int t1=0,int t2=0);
    std::vector<dvector> to2D(std::vector<dvector> &nXYZ);
}

#endif