#include "..\include\elem.h"
using namespace std;

// Class of cubic linear 3D element

Elem::Elem(vector<dvector> nXYZ,shapeStruct shape){

    matrix J;
    matrix invJ;
    matrix coord;

    // Memory allocation and initialization

    nLen = shape.nLen;
    int gLen = shape.gLen;
    dxN.setlength(nLen,gLen);
    dyN.setlength(nLen,gLen);
    dzN.setlength(nLen,gLen);
    coord.setlength(nLen,3);
    invJ.setlength(9,gLen);
    detJ.setlength(gLen);
    J.setlength(9,gLen);

    math::zero(J);
    math::zero(dxN);
    math::zero(dyN);
    math::zero(dzN);
    
    // Builds the Jacobian matrix and global coordinates

    for(int i=0; i<nLen; i++){

        coord(i,0) = nXYZ[i][0];
        coord(i,1) = nXYZ[i][1];
        coord(i,2) = nXYZ[i][2];

        for(int j=0; j<gLen; j++){

            // Jacobian matrix at Gauss points

            J(0,j) += shape.drN(i,j)*nXYZ[i][0];
            J(1,j) += shape.drN(i,j)*nXYZ[i][1];
            J(2,j) += shape.drN(i,j)*nXYZ[i][2];
            J(3,j) += shape.dsN(i,j)*nXYZ[i][0];
            J(4,j) += shape.dsN(i,j)*nXYZ[i][1];
            J(5,j) += shape.dsN(i,j)*nXYZ[i][2];
            J(6,j) += shape.dtN(i,j)*nXYZ[i][0];
            J(7,j) += shape.dtN(i,j)*nXYZ[i][1];
            J(8,j) += shape.dtN(i,j)*nXYZ[i][2];
        }
    }

    matrix v1 = math::prod(1,shape.N,coord,1);

    for(int i=0; i<gLen; i++){

        // Global coordinates of Gauss points

        gXYZ.push_back({v1(i,0),v1(i,1),v1(i,2)});

        // Determinant of the Jacobian matrix

        detJ(i) = J[0][i]*(J[4][i]*J[8][i]-J[5][i]*J[7][i]);
        detJ(i) + J[1][i]*(J[5][i]*J[6][i]-J[3][i]*J[8][i]);
        detJ(i) + J[2][i]*(J[3][i]*J[7][i]-J[4][i]*J[6][i]);

        // Inverse of the Jacobian matrix

        invJ(0,i) = (J[4][i]*J[8][i]-J[5][i]*J[7][i])/detJ(i);
        invJ(1,i) = (J[2][i]*J[7][i]-J[1][i]*J[8][i])/detJ(i);
        invJ(2,i) = (J[1][i]*J[5][i]-J[2][i]*J[4][i])/detJ(i);
        invJ(3,i) = (J[5][i]*J[6][i]-J[3][i]*J[8][i])/detJ(i);
        invJ(4,i) = (J[0][i]*J[8][i]-J[2][i]*J[6][i])/detJ(i);
        invJ(5,i) = (J[2][i]*J[3][i]-J[0][i]*J[5][i])/detJ(i);
        invJ(6,i) = (J[3][i]*J[7][i]-J[4][i]*J[6][i])/detJ(i);
        invJ(7,i) = (J[1][i]*J[6][i]-J[0][i]*J[7][i])/detJ(i);
        invJ(8,i) = (J[0][i]*J[4][i]-J[1][i]*J[3][i])/detJ(i);

        // Global derivatives of shape functions

        for(int j=0; j<nLen; j++){

            dxN(j,i) += shape.drN(j,i)*invJ(0,i)+shape.dsN(j,i)*invJ(1,i)+shape.dtN(j,i)*invJ(2,i);
            dyN(j,i) += shape.drN(j,i)*invJ(3,i)+shape.dsN(j,i)*invJ(4,i)+shape.dtN(j,i)*invJ(5,i);
            dzN(j,i) += shape.drN(j,i)*invJ(6,i)+shape.dsN(j,i)*invJ(7,i)+shape.dtN(j,i)*invJ(8,i);
        }
    }
}

// Computes the elemental local stiffness matrix

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
            
            B(j,0) = B(j+nLen,3) = B(j+2*nLen,4) = dxN(j,i);
            B(j,3) = B(j+nLen,1) = B(j+2*nLen,5) = dyN(j,i);
            B(j,4) = B(j+nLen,5) = B(j+2*nLen,2) = dzN(j,i);
        }
        matrix K1 = math::prod(quad.weight[i],B,D);
        matrix K2 = math::prod(detJ[i],K1,B,0,1);
        math::add(1,1,K2,K);
    }
    return K;
}

// Computes the elemental local strain matrix

matrix Elem::selfS(quadStruct quad,dvector xyz){

    matrix S;
    S.setlength(6,3*nLen);
    int gLen = quad.weight.size();
    math::zero(S);

    // Performs the numerical integration

    for(int i=0; i<gLen; i++){
        double k = math::kernel(xyz,gXYZ[i])*quad.weight[i]*detJ[i];

        for(int j=0; j<nLen; j++){

            S(0,j) = S(3,j+nLen) = S(4,j+2*nLen) += k*dxN(j,i);
            S(3,j) = S(1,j+nLen) = S(5,j+2*nLen) += k*dyN(j,i);
            S(4,j) = S(5,j+nLen) = S(2,j+2*nLen) += k*dzN(j,i);
        }   
    }
    return S;
}

// Class of square linear 2D element

Face::Face(vector<dvector> nXYZ,shapeStruct shape){

    matrix J,invJ;
    vector<dvector> nXY = math::to2D(nXYZ);

    // Memory allocation and initialization

    nLen = shape.nLen;
    int gLen = shape.gLen;
    dxN.setlength(nLen,gLen);
    dyN.setlength(nLen,gLen);
    invJ.setlength(4,gLen);
    detJ.setlength(gLen);
    J.setlength(4,gLen);

    math::zero(J);
    math::zero(dxN);
    math::zero(dyN);
    
    // Builds the Jacobian matrix

    for(int i=0; i<nLen; i++){
        for(int j=0; j<gLen; j++){

            J(0,j) += shape.drN(i,j)*nXY[i][0];
            J(1,j) += shape.drN(i,j)*nXY[i][1];
            J(2,j) += shape.dsN(i,j)*nXY[i][0];
            J(3,j) += shape.dsN(i,j)*nXY[i][1];
        }
    }

    for(int i=0; i<gLen; i++){

        // Determinant and inverse of the Jacobian matrix

        detJ(i) = J[0][i]*J[3][i]-J[1][i]*J[2][i];
        invJ(1,i) = -J[1][i]/detJ(i);
        invJ(2,i) = -J[2][i]/detJ(i);
        invJ(0,i) = J[3][i]/detJ(i);
        invJ(3,i) = J[0][i]/detJ(i);

        // Global derivatives of shape functions

        for(int j=0; j<nLen; j++){

            dxN(j,i) += shape.drN(j,i)*invJ(0,i)+shape.dsN(j,i)*invJ(1,i);
            dyN(j,i) += shape.drN(j,i)*invJ(2,i)+shape.dsN(j,i)*invJ(3,i);
        }
    }
}

// Integrates the shape functions matrix

matrix Face::selfN(shapeStruct shape,quadStruct quad){

    matrix N; N.setlength(3*nLen,3);
    int gLen = quad.weight.size();
    math::zero(N);

    // Builds the shape function matrix

    for(int i=0; i<nLen; i++){
        for(int j=0; j<gLen; j++){
            
            N(i,0) += shape.N(i,j)*quad.weight[j]*detJ(j);
            N(4+i,1) += shape.N(i,j)*quad.weight[j]*detJ(j);
            N(8+i,2) += shape.N(i,j)*quad.weight[j]*detJ(j);
        }
    }
    return N;
}