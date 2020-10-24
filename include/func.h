#include "linalg.h"
#include "stdafx.h"
#ifndef FUNC_H
#define FUNC_H

typedef alglib::sparsematrix sparse;
typedef alglib::real_2d_array matrix;
typedef alglib::real_1d_array darray;
typedef alglib::integer_1d_array iarray;

struct quadStruct{

    std::vector<double> weight;
    std::vector<darray> gRST;
};

namespace math{

    void zero(darray &V);
    void zero(matrix &M);
    void clean(sparse &M);
    quadStruct quadrature(int dim);
    darray cross(darray &V1,darray &V2);
    void add(double k1,double k2,matrix &M1,matrix &M2);
    std::vector<darray> to2D(std::vector<darray> &nXYZ);
    darray prod(double k,matrix &M1,darray &V1,int t1=0);
    matrix prod(double k,matrix &M1,matrix &M2,int t1=0,int t2=0);
}

#endif