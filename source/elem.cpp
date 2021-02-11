#include "..\include\elem.h"
using namespace std;

// --------------------------------------------------|
// Class of 8-node hexahedron linear element in 3D   |
// --------------------------------------------------|

Elem::Elem(vector<dvector> nXYZ,shapeStruct shape){

    matrix J;
    matrix invJ;
    matrix coord;

    // Number of nodes and Gauss points

    nLen = shape.N.rows();
    int gLen = shape.N.cols();

    // Memory allocation and initialization

    dxN.setlength(nLen,gLen);
    dyN.setlength(nLen,gLen);
    dzN.setlength(nLen,gLen);
    coord.setlength(nLen,3);
    invJ.setlength(9,gLen);
    detJ.setlength(gLen);
    J.setlength(9,gLen);

    // Fills the matrices with zeros

    math::zero(J);
    math::zero(dxN);
    math::zero(dyN);
    math::zero(dzN);
    
    // Builds the Jacobian matrix and global coordinates

    for(int i=0; i<nLen; i++){

        coord(i,0) = nXYZ[i][0];
        coord(i,1) = nXYZ[i][1];
        coord(i,2) = nXYZ[i][2];

        for(int k=0; k<3; k++){

            // Jacobian matrix at Gauss points

            for(int j=0; j<gLen; j++){
                
                J(k,j) += shape.drN(i,j)*nXYZ[i][k];
                J(k+3,j) += shape.dsN(i,j)*nXYZ[i][k];
                J(k+6,j) += shape.dtN(i,j)*nXYZ[i][k];
            }
        }
    }

    matrix v = math::prod(1,shape.N,coord,1);
    for(int i=0; i<gLen; i++){

        // Global coordinates of Gauss points

        gXYZ.push_back({v(i,0),v(i,1),v(i,2)});

        // Determinant of the Jacobian matrix

        detJ(i) = J[0][i]*(J[4][i]*J[8][i]-J[5][i]*J[7][i]);
        detJ(i) += J[1][i]*(J[5][i]*J[6][i]-J[3][i]*J[8][i]);
        detJ(i) += J[2][i]*(J[3][i]*J[7][i]-J[4][i]*J[6][i]);

        // Inverse of the Jacobian matrix

        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){

                invJ(3*j+k,i) = J[(j+1)%3*3+(k+1)%3][i]*J[(j+2)%3*3+(k+2)%3][i];
                invJ(3*j+k,i) -= J[(j+2)%3*3+(k+1)%3][i]*J[(j+1)%3*3+(k+2)%3][i];
                invJ(3*j+k,i) /= detJ(i);
            }
        }

        // Global derivatives of shape functions

        for(int j=0; j<nLen; j++){

            dxN(j,i) += shape.drN(j,i)*invJ(0,i)+shape.dsN(j,i)*invJ(1,i)+shape.dtN(j,i)*invJ(2,i);
            dyN(j,i) += shape.drN(j,i)*invJ(3,i)+shape.dsN(j,i)*invJ(4,i)+shape.dtN(j,i)*invJ(5,i);
            dzN(j,i) += shape.drN(j,i)*invJ(6,i)+shape.dsN(j,i)*invJ(7,i)+shape.dtN(j,i)*invJ(8,i);
        }
    }
}

// ---------------------------------------------------------|
// Computes the elemental stiffness matrix K for local FEM  |
// ---------------------------------------------------------|

matrix Elem::selfK(quadStruct quad,matrix D){

    matrix B,K;
    B.setlength(3*nLen,6);
    K.setlength(3*nLen,3*nLen);
    int gLen = quad.weight.size();
    math::zero(B);
    math::zero(K);

    // Performs the numerical integration

    for(int i=0; i<gLen; i++){
        for(int j=0; j<nLen; j++){

            // Computes the shape functions derivative matrix B
            
            B(j,0) = B(j+nLen,3) = B(j+2*nLen,4) = dxN(j,i);
            B(j,3) = B(j+nLen,1) = B(j+2*nLen,5) = dyN(j,i);
            B(j,4) = B(j+nLen,5) = B(j+2*nLen,2) = dzN(j,i);
        }

        // Computes K by Gauss-Legendre quadrature

        matrix K1 = math::prod(quad.weight[i],B,D);
        matrix K2 = math::prod(detJ[i],K1,B,0,1);
        math::add(1,1,K2,K);
    }
    return K;
}

// -------------------------------------------------------|
// Computes the averaged Von Mises stress in the element  |
// -------------------------------------------------------|

double Elem::stress(quadStruct quad,matrix D,darray u){

    matrix B;
    B.setlength(3*nLen,6);
    int gLen = quad.weight.size();
    double VM=0,volume=0;
    math::zero(B);

    // Performs the numerical integration

    for(int i=0; i<gLen; i++){
        for(int j=0; j<nLen; j++){

            // Computes the shape functions derivative matrix B
            
            B(j,0) = B(j+nLen,3) = B(j+2*nLen,4) = dxN(j,i);
            B(j,3) = B(j+nLen,1) = B(j+2*nLen,5) = dyN(j,i);
            B(j,4) = B(j+nLen,5) = B(j+2*nLen,2) = dzN(j,i);
        }

        // Computes the stress field at Gauss points

        darray strain = math::prod(1,B,u,1);
        darray sigma = math::prod(1,D,strain,0);

        // Integrates the Von Mises stress in the element and the volume

        double VM1 = 3*(sigma[3]*sigma[3]+sigma[4]*sigma[4]+sigma[5]*sigma[5]);
        VM1 += (pow(sigma[0]-sigma[1],2)+pow(sigma[1]-sigma[2],2)+pow(sigma[2]-sigma[0],2))/2;
        VM += quad.weight[i]*sqrt(VM1)*detJ[i];
        volume += quad.weight[i]*detJ[i];
    }

    VM /= volume;
    return VM;
}

// ---------------------------------------------------|
// Computes the elemental S matrix for non-local FEM  |
// ---------------------------------------------------|
/*
matrix Elem::selfS(quadStruct quad,dvector xyz){

    matrix S;
    S.setlength(6,3*nLen);
    int gLen = quad.weight.size();
    math::zero(S);

    // Performs the numerical integration

    for(int i=0; i<gLen; i++){
        double k = math::kernel(xyz,gXYZ[i])*quad.weight[i]*detJ[i];

        for(int j=0; j<nLen; j++){

            // Computes the non-local S matrix

            S(0,j) = S(3,j+nLen) = S(4,j+2*nLen) += k*dxN(j,i);
            S(3,j) = S(1,j+nLen) = S(5,j+2*nLen) += k*dyN(j,i);
            S(4,j) = S(5,j+nLen) = S(2,j+2*nLen) += k*dzN(j,i);
        }   
    }
    return S;
}
*/
// --------------------------------------------------|
// Class of 4-node quadrangle linear element in 2D   |
// --------------------------------------------------|

Face::Face(vector<dvector> nXYZ,shapeStruct shape){

    int gLen = shape.N.cols();
    nLen = shape.N.rows();
    dJ2D.setlength(gLen);
    math::zero(dJ2D);

    // Performs a loop over the Gauss points

    for(int j=0; j<gLen; j++){

        dvector J2Dr(3,0);
        dvector J2Ds(3,0);

        // Computes the partial Jacobian vectors for 3D to 2D mapping

        for(int k=0; k<3; k++){
            for(int i=0; i<nLen; i++){
            
                J2Dr[k] += shape.drN(i,j)*nXYZ[i][k];
                J2Ds[k] += shape.dsN(i,j)*nXYZ[i][k];
            }
        }

        // Computes the cross product and takes the norm

        dvector dJ = math::cross(J2Dr,J2Ds);
        for(int i=0; i<3; i++){dJ2D(j) += dJ[i]*dJ[i];}
        dJ2D(j) = sqrt(dJ2D[j]);
    }
}

// --------------------------------------------------------------|
// Integrates the matrix of shape functions N over the element   |
// --------------------------------------------------------------|

matrix Face::selfM(shapeStruct shape,quadStruct quad){

    matrix M; M.setlength(3*nLen,3);
    int gLen = shape.N.cols();
    math::zero(M);

    // Integrates the matrix of local shape functions

    for(int i=0; i<nLen; i++){
        for(int j=0; j<gLen; j++){
            
            M(i,0) += shape.N(i,j)*quad.weight[j]*dJ2D(j);
            M(4+i,1) += shape.N(i,j)*quad.weight[j]*dJ2D(j);
            M(8+i,2) += shape.N(i,j)*quad.weight[j]*dJ2D(j);
        }
    }
    return M;
}