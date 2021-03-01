#include "..\include\elem.h"
using namespace std;

// -----------------------------------------------------------|
// Class of hexahedron isoparametric lagrange element in 3D   |
// -----------------------------------------------------------|

Elem::Elem(vector<array3d> nXYZ,ivector surface){

    // Builds the Jacobian matrix and its inverse

    this->nXYZ = nXYZ;
    nLen = nXYZ.size();
    this->surface = surface;
}

// -------------------------------------------------------|
// Builds the Jacobian matrix evaluated at Gauss points   |
// -------------------------------------------------------|

vector<matrix> Elem::jacobian(shapeStruct shape){

    vector<matrix> J(shape.gLen);

    // Builds the Jacobian matrix and global coordinates

    for(int i=0; i<shape.gLen; i++){

        J[i].setlength(3,3);
        math::zero(J[i]);

        // Computes the Jacobian matrix at this Gauss point

        for(int j=0; j<nLen; j++){
            for(int k=0; k<3; k++){
                for(int n=0; n<3; n++){
                    J[i](n,k) += shape.dN[n](j,i)*nXYZ[j][k];
                }
            }
        }
    }
    return J;
}

// -------------------------------------------------------------|
// Determinant of J and global derivatives of shape functions   |
// -------------------------------------------------------------|

void Elem::derivative(shapeStruct shape,vector<matrix> J){

    detJ.resize(shape.gLen);
    vector<matrix> invJ(shape.gLen);

    // Resets the shape function derivative

    for(int i=0; i<3; i++){
        
        dN[i].setlength(nLen,shape.gLen);
        math::zero(dN[i]);
    }

    // Determinant and inverse of the Jacobian matrix

    for(int i=0; i<shape.gLen; i++){

        detJ[i] = alglib::rmatrixdet(J[i],3);
        invJ[i] = math::invert(J[i],detJ[i]);

        // Global derivatives of shape functions

        for(int j=0; j<nLen; j++){
            for(int k=0; k<3; k++){
                for(int n=0; n<3; n++){

                    dN[k](j,i) += shape.dN[n](j,i)*invJ[i](k,n);
                }
            }
        }
    }
}

// ---------------------------------------------------------|
// Computes the elemental stiffness matrix K for local FEM  |
// ---------------------------------------------------------|

matrix Elem::selfK(shapeStruct shape,matrix D){

    matrix B,K;
    B.setlength(3*nLen,6);
    K.setlength(3*nLen,3*nLen);
    math::zero(B);
    math::zero(K);

    // Computes the Jacobian and shape function derivatives

    vector<matrix> J = jacobian(shape);
    derivative(shape,J);

    // Performs the numerical integration

    for(int i=0; i<shape.gLen; i++){
        for(int j=0; j<nLen; j++){

            // Computes the shape functions derivative matrix B
            
            B(j,0) = B(j+nLen,5) = B(j+2*nLen,4) = dN[0](j,i);
            B(j,5) = B(j+nLen,1) = B(j+2*nLen,3) = dN[1](j,i);
            B(j,4) = B(j+nLen,3) = B(j+2*nLen,2) = dN[2](j,i);
        }

        // Computes K by Gauss-Legendre quadrature

        matrix K1 = math::prod(shape.weight[i],B,D);
        matrix K2 = math::prod(detJ[i],K1,B,0,1);
        math::add(1,1,K2,K);
    }
    return K;
}

// ----------------------------------------------------|
// Computes the elemental mass matrix M for local FEM  |
// ----------------------------------------------------|

matrix Elem::selfM(shapeStruct shape,double rho){

    matrix N,M;
    N.setlength(3*nLen,3);
    M.setlength(3*nLen,3*nLen);
    math::zero(N);
    math::zero(M);

    // Computes the Jacobian and shape function derivatives

    vector<matrix> J = jacobian(shape);
    derivative(shape,J);

    // Performs the numerical integration

    for(int i=0; i<shape.gLen; i++){
        for(int j=0; j<nLen; j++){
            N(j,0) = N(j+nLen,1) = N(j+2*nLen,2) = shape.N(j,i);
        }

        // Computes M by Gauss-Legendre quadrature

        double wdetJ = shape.weight[i]*detJ[i];
        matrix Me = math::prod(wdetJ,N,N,0,1);
        math::add(rho,1,Me,M);
    }
    return M;
}

// ----------------------------------------------------------------|
// Computes the elemental curl-displacemen matrix K for local FEM  |
// ----------------------------------------------------------------|

matrix Elem::selfG(shapeStruct shape){

    matrix B,K,N;
    N.setlength(3*nLen,3);
    B.setlength(3*nLen,3);
    K.setlength(3*nLen,3*nLen);
    math::zero(N);
    math::zero(B);
    math::zero(K);

    // Computes the Jacobian and shape function derivatives

    vector<matrix> J = jacobian(shape);
    derivative(shape,J);

    // Performs the numerical integration

    for(int i=0; i<shape.gLen; i++){
        for(int j=0; j<nLen; j++){

            // Computes the shape functions and derivative matrices
            
            B(j,1) = dN[2](j,i); B(j+nLen,0) = -dN[2](j,i);
            B(j,2) = -dN[1](j,i); B(j+2*nLen,0) = dN[1](j,i);
            B(j+nLen,2) = dN[0](j,i); B(j+2*nLen,1) = -dN[0](j,i);
            N(j,0) = N(j+nLen,1) = N(j+2*nLen,2) = shape.N(j,i);
        }

        // Computes K by Gauss-Legendre quadrature

        double wdetJ = shape.weight[i]*detJ[i];
        matrix Ke = math::prod(wdetJ,B,N,0,1);
        math::add(-0.5,1,Ke,K);
    }
    return K;
}

// ------------------------------------------------------------------|
// Computes the elemental surface stiffness matrix Ks for local FEM  |
// ------------------------------------------------------------------|

matrix Elem::selfKs(shapeStruct shape[6],matrix D){

    matrix B,K;
    array3d norm;
    B.setlength(3*nLen,6);
    K.setlength(3*nLen,3*nLen);
    math::zero(B);
    math::zero(K);

    // Computes the Jacobian and shape function derivatives




    for(int k:surface){

        vector<matrix> J = jacobian(shape[k]);
        derivative(shape[k],J);


        for(int i=0; i<shape[k].gLen; i++){
            for(int j=0; j<nLen; j++){

                // Computes the shape functions derivative matrix B
                
                B(j,0) = B(j+nLen,5) = B(j+2*nLen,4) = dN[0](j,i);
                B(j,5) = B(j+nLen,1) = B(j+2*nLen,3) = dN[1](j,i);
                B(j,4) = B(j+nLen,3) = B(j+2*nLen,2) = dN[2](j,i);
            }


            // Computes the shape derivative vectors

            array3d vr = {J[i](0,0),J[i](0,1),J[i](0,2)};
            array3d vs = {J[i](1,0),J[i](1,1),J[i](1,2)};
            array3d vt = {J[i](2,0),J[i](2,1),J[i](2,2)};
            array3d v;



            switch (k) {
            case 0:
                v = math::dotsub(vs,vt);
                norm = math::cross(v,vs,1);
                break;
            case 1:
                v = math::dotsub(vs,vt);
                norm = math::cross(vs,v,1);
                break;
            case 2:
                v = math::dotsub(vt,vr);
                norm = math::cross(vt,v,1);
                break;
            case 3:
                v = math::dotsub(vt,vr);
                norm = math::cross(v,vt,1);
                break;
            case 4:
                v = math::dotsub(vr,vs);
                norm = math::cross(v,vr,1);
                break;
            case 5:
                v = math::dotsub(vr,vs);
                norm = math::cross(vr,v,1);
                break;
            }

            matrix T,P;
            P.setlength(3,3);
            T.setlength(6,6);

            for(int j=0; j<3; j++){
                for(int k=0; k<3; k++){
                    
                    if(j==k){P(j,k) = 1-norm[j]*norm[k];}
                    else{P(j,k) = -norm[j]*norm[k];}
                }
            }

            T(0,0) = P(0,0)*P(0,0);
            T(0,1) = P(0,1)*P(0,1);
            T(0,2) = P(2,0)*P(2,0);
            T(1,0) = P(0,1)*P(0,1);
            T(1,1) = P(1,1)*P(1,1);
            T(1,2) = P(1,2)*P(1,2);
            T(2,0) = P(2,0)*P(2,0);
            T(2,1) = P(1,2)*P(1,2);
            T(2,2) = P(2,2)*P(2,2);

            T(0,3) = P(2,0)*P(0,1);
            T(0,4) = P(2,0)*P(0,0);
            T(0,5) = P(0,1)*P(0,0);
            T(1,3) = P(1,2)*P(1,1);
            T(1,4) = P(1,2)*P(0,1);
            T(1,5) = P(1,1)*P(0,1);
            T(2,3) = P(1,2)*P(2,2);
            T(2,4) = P(2,2)*P(2,0);
            T(2,5) = P(1,2)*P(2,0);

            T(3,0) = 2*P(2,0)*P(0,1);
            T(3,1) = 2*P(1,2)*P(1,1);
            T(3,2) = 2*P(2,2)*P(1,2);
            T(4,0) = 2*P(2,0)*P(0,0);
            T(4,1) = 2*P(1,2)*P(0,1);
            T(4,2) = 2*P(1,2)*P(2,2);
            T(5,0) = 2*P(0,1)*P(0,0);
            T(5,1) = 2*P(1,1)*P(0,1);
            T(5,2) = 2*P(2,0)*P(1,2);

            T(3,3) = P(2,2)*P(1,1)+P(1,2)*P(1,2);
            T(3,4) = P(2,2)*P(0,1)+P(2,0)*P(1,1);
            T(3,5) = P(2,0)*P(1,1)+P(0,1)*P(1,2);
            T(4,3) = P(1,2)*P(2,0)+P(2,2)*P(0,1);
            T(4,4) = P(2,2)*P(0,0)+P(2,0)*P(2,0);
            T(4,5) = P(2,0)*P(0,1)+P(1,2)*P(0,0);
            T(5,3) = P(1,2)*P(0,1)+P(2,0)*P(1,1);
            T(5,4) = P(2,0)*P(0,1)+P(0,0)*P(1,2);
            T(5,5) = P(0,0)*P(1,1)+P(0,1)*P(0,1);

            // Computes the surface stiffness tensor

            matrix S = math::prod(1,T,D,1);
            S = math::prod(1,S,T);

            // Computes K by Gauss-Legendre quadrature

            matrix BT = math::prod(1,B,T,0,1);
            matrix K1 = math::prod(shape[k].weight[i],BT,S);
            matrix K2 = math::prod(detJ[i],K1,BT,0,1);
            math::add(1,1,K2,K);
        }
    }
    return K;
}

// -------------------------------------------------------|
// Computes the averaged Von Mises stress in the element  |
// -------------------------------------------------------|

darray Elem::stress(shapeStruct shape,matrix D,darray u){

    matrix B;
    darray sigma;
    double volume = 0;
    sigma.setlength(6);
    B.setlength(3*nLen,6);
    math::zero(sigma);
    math::zero(B);

    // Computes the Jacobian and shape function derivatives

    vector<matrix> J = jacobian(shape);
    derivative(shape,J);

    // Performs the numerical integration

    for(int i=0; i<shape.gLen; i++){
        for(int j=0; j<nLen; j++){

            // Computes the shape functions derivative matrix B
            
            B(j,0) = B(j+nLen,5) = B(j+2*nLen,4) = dN[0](j,i);
            B(j,5) = B(j+nLen,1) = B(j+2*nLen,3) = dN[1](j,i);
            B(j,4) = B(j+nLen,3) = B(j+2*nLen,2) = dN[2](j,i);
        }

        // Computes the stress field at Gauss points

        darray strain = math::prod(1,B,u,1);
        darray stress = math::prod(1,D,strain,0);

        // Integrates the Von Mises stress in the element and the volume

        math::add(shape.weight[i]*detJ[i],1,stress,sigma);
        volume += shape.weight[i]*detJ[i];
    }

    for(int i=0; i<6; i++){sigma[i] /= volume;}
    return sigma;
}

// -----------------------------------------------------------|
// Class of quadrangle isoparametric lagrange element in 2D   |
// -----------------------------------------------------------|

Face::Face(vector<array3d> nXYZ,shapeStruct shape){

    nLen = shape.N.cols();
    dJ2D.setlength(shape.gLen);
    math::zero(dJ2D);

    // Performs a loop over the Gauss points

    for(int j=0; j<shape.gLen; j++){

        array3d J2Dr = {0,0,0};
        array3d J2Ds = {0,0,0};

        // Computes the partial Jacobian vectors for 3D to 2D mapping

        for(int k=0; k<3; k++){
            for(int i=0; i<nLen; i++){
            
                J2Dr[k] += shape.dN[0](i,j)*nXYZ[i][k];
                J2Ds[k] += shape.dN[1](i,j)*nXYZ[i][k];
            }
        }

        // Computes the cross product and takes the norm

        array3d dJ = math::cross(J2Dr,J2Ds);
        for(int i=0; i<3; i++){dJ2D(j) += dJ[i]*dJ[i];}
        dJ2D(j) = sqrt(dJ2D[j]);
    }
}

// --------------------------------------------------------------|
// Integrates the matrix of shape functions N over the element   |
// --------------------------------------------------------------|

darray Face::selfB(shapeStruct shape,darray F){

    matrix N;
    darray B;
    B.setlength(3*nLen);
    N.setlength(3*nLen,3);
    math::zero(N);
    math::zero(B);

    // Integrates the matrix of local shape functions

    for(int i=0; i<shape.gLen; i++){
        for(int j=0; j<nLen; j++){
            N(j,0) = N(j+nLen,1) = N(j+2*nLen,2) = shape.N(j,i);
        }

        // Computes M by Gauss-Legendre quadrature

        double wdetJ = shape.weight[i]*dJ2D(i);
        darray Be = math::prod(wdetJ,N,F,0);
        math::add(1,1,Be,B);
    }
    return B;
}