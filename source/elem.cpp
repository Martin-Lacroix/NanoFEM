#include "..\include\elem.h"
using namespace std;

// -----------------------------------------------------------|
// Class of hexahedron isoparametric Lagrange element in 3D   |
// -----------------------------------------------------------|

Elem::Elem(vector<array3d> inp1,ivector inp2) : nXYZ{inp1},surface{inp2}{

    // Sets the number of nodes in each dimension

    nLen = nXYZ.size();
    sLen = cbrt(nLen+1e-5);
}

// ---------------------------------------------------------|
// Jacobian and derivative of shape functions in the bulk   |
// ---------------------------------------------------------|

void Elem::updateJ(shapeStruct &shape){

    detJ.resize(shape.gLen);
    vector<matrix> J(shape.gLen);
    vector<matrix> invJ(shape.gLen);

    // Resets the shape function derivative

    for(int i=0; i<3; i++){
        
        dN[i].setlength(nLen,shape.gLen);
        math::zero(dN[i]);
    }

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

    // Determinant and inverse of the Jacobian matrix

    for(int i=0; i<shape.gLen; i++){

        detJ[i] = alglib::rmatrixdet(J[i],3);
        invJ[i] = math::invert(J[i],detJ[i]);

        // Global derivatives of shape functions

        for(int j=0; j<nLen; j++){
            for(int k=0; k<3; k++){
                for(int n=0; n<3; n++){
                    dN[k](j,i) += shape.dN[n](j,i)*invJ[i](n,k);
                }
            }
        }
    }
}

// ------------------------------------------------------------|
// Jacobian and derivative of shape functions at the surface   |
// ------------------------------------------------------------|

void Elem::updateS(shapeStruct (&shape)[6]){

    for(int s:surface){

        array3d v;
        norm[s].resize(shape[s].gLen);
        detJ2D[s].resize(shape[s].gLen);
        vector<matrix> J(shape[s].gLen);
        vector<matrix> invJ(shape[s].gLen);

        // Resets the shape function derivative

        for(int i=0; i<3; i++){
            
            dNs[s][i].setlength(nLen,shape[s].gLen);
            math::zero(dNs[s][i]);
        }

        // Builds the Jacobian matrix and global coordinates

        for(int i=0; i<shape[s].gLen; i++){

            J[i].setlength(3,3);
            math::zero(J[i]);

            // Computes the Jacobian matrix at this Gauss point

            for(int j=0; j<nLen; j++){
                for(int k=0; k<3; k++){
                    for(int n=0; n<3; n++){
                        J[i](n,k) += shape[s].dN[n](j,i)*nXYZ[j][k];
                    }
                }
            }

            // Normal to the s-th surface of the element

            array3d vr = {J[i](0,0),J[i](0,1),J[i](0,2)};
            array3d vs = {J[i](1,0),J[i](1,1),J[i](1,2)};
            array3d vt = {J[i](2,0),J[i](2,1),J[i](2,2)};

            switch(s){
            case 0:
                v = math::dotsub(vs,vr);
                norm[s][i] = math::cross(vs,v);
                break;
            case 1:
                v = math::dotsub(vs,vr);
                norm[s][i] = math::cross(v,vs);
                break;
            case 2:
                v = math::dotsub(vr,vt);
                norm[s][i] = math::cross(v,vr);
                break;
            case 3:
                v = math::dotsub(vr,vt);
                norm[s][i] = math::cross(vr,v);
                break;
            case 4:
                v = math::dotsub(vt,vs);
                norm[s][i] = math::cross(vt,v);
                break;
            case 5:
                v = math::dotsub(vt,vs);
                norm[s][i] = math::cross(v,vt);
                break;
            }
        }

        // Determinant and inverse of the Jacobian matrix

        for(int i=0; i<shape[s].gLen; i++){

            double detJ = alglib::rmatrixdet(J[i],3);
            invJ[i] = math::invert(J[i],detJ);

            // 2D Jacobian determinant and the surface normal

            detJ2D[s][i] = math::norm(norm[s][i]);
            for(int j=0; j<3; j++){norm[s][i][j] /= detJ2D[s][i];}

            // Global derivatives of shape functions

            for(int j=0; j<nLen; j++){
                for(int k=0; k<3; k++){
                    for(int n=0; n<3; n++){
                        dNs[s][k](j,i) += shape[s].dN[n](j,i)*invJ[i](n,k);
                    }
                }
            }
        }
    }
}

// --------------------------------------------------------|
// This should normally free the memory of those vectors   |
// --------------------------------------------------------|

void Elem::clean(){

    // Clear bulk Jacobian and parameters

    dvector().swap(detJ);
    vector<darray>().swap(E);
    vector<matrix>().swap(F);
    for(int i=0; i<3; i++){dN[i].setlength(0,0);}

    // Clear surface Jacobian and parameters

    for(int s:surface){
        
        dvector().swap(detJ2D[s]);
        vector<array3d>().swap(norm[s]);
        for(int i=0; i<3; i++){dNs[s][i].setlength(0,0);}
    }
}

// ----------------------------------------------------------|
// Computes the elemental stiffness matrix K for local FEM   |
// ----------------------------------------------------------|

matrix Elem::selfK(shapeStruct &shape,array3d LmR){

    matrix B,K;
    B.setlength(3*nLen,6);
    K.setlength(3*nLen,3*nLen);
    matrix D = math::stiffness(LmR);
    math::zero(B);
    math::zero(K);

    // Performs the numerical integration

    for(int i=0; i<shape.gLen; i++){
        for(int j=0; j<nLen; j++){

            // Computes the shape functions derivative matrix B
            
            B(j,0) = B(j+nLen,3) = B(j+2*nLen,5) = dN[0](j,i);
            B(j,3) = B(j+nLen,1) = B(j+2*nLen,4) = dN[1](j,i);
            B(j,5) = B(j+nLen,4) = B(j+2*nLen,2) = dN[2](j,i);
        }
        
        // Computes K by Gauss-Legendre quadrature

        matrix K1 = math::prod(shape.weight[i],B,D);
        matrix K2 = math::prod(detJ[i],K1,B,0,1);
        math::add(1,1,K2,K);
    }
    return K;
}

// -----------------------------------------------------|
// Computes the elemental mass matrix M for local FEM   |
// -----------------------------------------------------|

matrix Elem::selfM(shapeStruct &shape,double rho){

    matrix N,M;
    N.setlength(3*nLen,3);
    M.setlength(3*nLen,3*nLen);
    math::zero(N);
    math::zero(M);

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

// ---------------------------------------------------------|
// Computes the surface stiffness matrix Ks for local FEM   |
// ---------------------------------------------------------|

matrix Elem::selfKS(shapeStruct (&shape)[6],array3d LmS){

    matrix B,K;
    B.setlength(3*nLen,6);
    K.setlength(3*nLen,3*nLen);
    matrix D = math::stiffness(LmS);
    math::zero(K);

    // Update the Jacobian and numerical integration

    for(int k:surface){
        math::zero(B);

        for(int i=0; i<shape[k].gLen; i++){
            for(int j=0; j<nLen; j++){

                // Computes the shape functions derivative matrix B
                
                B(j,0) = B(j+nLen,3) = B(j+2*nLen,5) = dNs[k][0](j,i);
                B(j,3) = B(j+nLen,1) = B(j+2*nLen,4) = dNs[k][1](j,i);
                B(j,5) = B(j+nLen,4) = B(j+2*nLen,2) = dNs[k][2](j,i);
            }

            // Computes the isotropic surface stiffness tensor

            matrix T = math::projection(norm[k][i]);
            matrix S = math::prod(1,T,D);
            S = math::prod(1,S,T,0,1);

            // Computes K by Gauss-Legendre quadrature

            matrix BT = math::prod(1,B,T);
            matrix K1 = math::prod(shape[k].weight[i],BT,S);
            matrix K2 = math::prod(detJ2D[k][i],K1,BT,0,1);
            math::add(1,1,K2,K);
        }
    }
    return K;
}

// --------------------------------------------------------|
// Computes the surface boundary vector Fs for local FEM   |
// --------------------------------------------------------|

darray Elem::selfFS(shapeStruct (&shape)[6],array3d LmS){

    matrix B;
    darray F;
    F.setlength(3*nLen);
    B.setlength(3*nLen,6);
    math::zero(F);

    // Surface traction vector and stiffness matrix

    darray tau;
    double tauS[6] = {LmS[2],LmS[2],LmS[2],0,0,0};
    tau.setcontent(6,tauS);

    // Update the Jacobian and numerical integration

    for(int s:surface){
        math::zero(B);

        for(int i=0; i<shape[s].gLen; i++){
            for(int j=0; j<nLen; j++){

                // Computes the shape functions derivative matrix B
                
                B(j,0) = B(j+nLen,3) = B(j+2*nLen,5) = dNs[s][0](j,i);
                B(j,3) = B(j+nLen,1) = B(j+2*nLen,4) = dNs[s][1](j,i);
                B(j,5) = B(j+nLen,4) = B(j+2*nLen,2) = dNs[s][2](j,i);
            }

            // Computes the surface tension by Gauss-Legendre quadrature

            matrix T = math::projection(norm[s][i]);
            matrix BT = math::prod(1,B,T);
            darray Fe = math::prod(shape[s].weight[i],BT,tau);
            alglib::vsub(&F[0],&Fe[0],3*nLen,detJ2D[s][i]);
        }
    }
    return F;
}

// -----------------------------------------------------|
// Computes the elemental deformation gradient tensor   |
// -----------------------------------------------------|

void Elem::updateF(shapeStruct &shape,darray u){

    matrix B;
    B.setlength(3*nLen,9);
    math::zero(B);

    // Update the Jacobian and the deformation gradient

    F.resize(shape.gLen);
    E.resize(shape.gLen);

    for(int i=0; i<shape.gLen; i++){
        for(int j=0; j<nLen; j++){

            // Computes the linear shape functions derivative matrix
            
            B(j,0) = B(j+nLen,4) = B(j+2*nLen,7) = dN[0](j,i);
            B(j,3) = B(j+nLen,1) = B(j+2*nLen,6) = dN[1](j,i);
            B(j,8) = B(j+nLen,5) = B(j+2*nLen,2) = dN[2](j,i);
        }

        // Stores the deformation gradient tensor

        F[i].setlength(3,3);
        darray Fv =  math::prod(1,B,u,1);

        F[i](0,0) = Fv[0]+1;
        F[i](1,1) = Fv[1]+1;
        F[i](2,2) = Fv[2]+1;
        F[i](0,1) = Fv[3];
        F[i](1,0) = Fv[4];
        F[i](1,2) = Fv[5];
        F[i](2,1) = Fv[6];
        F[i](2,0) = Fv[7];
        F[i](0,2) = Fv[8];

        // Stores the Green-Lagrange strain tensor

        E[i].setlength(6);
        matrix Ev = math::prod(1,F[i],F[i],1,0);

        E[i](0) = (Ev[0][0]-1)/2;
        E[i](1) = (Ev[1][1]-1)/2;
        E[i](2) = (Ev[2][2]-1)/2;
        E[i](3) = Ev[0][1]/2;
        E[i](4) = Ev[1][2]/2;
        E[i](5) = Ev[2][0]/2;
    }
}

// -------------------------------------------------------------|
// Computes the non-linear stiffness matrix for finite strain   |
// -------------------------------------------------------------|

matrix Elem::selfKN(shapeStruct &shape,array3d LmR){

    matrix K,B;
    B.setlength(3*nLen,6);
    K.setlength(3*nLen,3*nLen);
    matrix D = math::stiffness(LmR);
    math::zero(K);
    math::zero(B);

    // Performs the numerical integration

    for(int i=0; i<shape.gLen; i++){
        for(int j=0; j<nLen; j++){

            // Computes the non-linear shape functions derivative matrix

            for(int k=0; k<3; k++){
            
                B(j+k*nLen,0) = F[i](k,0)*dN[0](j,i);
                B(j+k*nLen,1) = F[i](k,1)*dN[1](j,i);
                B(j+k*nLen,2) = F[i](k,2)*dN[2](j,i);

                B(j+k*nLen,3) = F[i](k,0)*dN[1](j,i)+F[i](k,1)*dN[0](j,i);
                B(j+k*nLen,4) = F[i](k,1)*dN[2](j,i)+F[i](k,2)*dN[1](j,i);
                B(j+k*nLen,5) = F[i](k,2)*dN[0](j,i)+F[i](k,0)*dN[2](j,i);
            }
        }

        // Computes K by Gauss-Legendre quadrature

        matrix K1 = math::prod(shape.weight[i],B,D);
        matrix K2 = math::prod(detJ[i],K1,B,0,1);
        math::add(1,1,K2,K);
    }
    return K;
}

// ---------------------------------------------------------|
// Computes the linear stiffness matrix for finite strain   |
// ---------------------------------------------------------|

matrix Elem::selfKL(shapeStruct &shape,array3d LmR){

    matrix K,B,S;
    S.setlength(9,9);
    B.setlength(3*nLen,9);
    K.setlength(3*nLen,3*nLen);
    matrix D = math::stiffness(LmR);
    math::zero(K);
    math::zero(B);
    math::zero(S);

    // Performs the numerical integration

    for(int i=0; i<shape.gLen; i++){
        for(int j=0; j<nLen; j++){

            // Computes the linear shape functions derivative matrix
            
            B(j,0) = B(j+nLen,4) = B(j+2*nLen,7) = dN[0](j,i);
            B(j,3) = B(j+nLen,1) = B(j+2*nLen,6) = dN[1](j,i);
            B(j,8) = B(j+nLen,5) = B(j+2*nLen,2) = dN[2](j,i);
        }

        // Computes the matrix of second Piola-Kirchhoff components

        darray Sv = math::prod(1,D,E[i]);

        S(0,0) = S(4,4) = S(7,7) = Sv[0];
        S(1,1) = S(3,3) = S(6,6) = Sv[1];
        S(2,2) = S(5,5) = S(8,8) = Sv[2];
        S(0,3) = S(1,4) = S(6,7) = Sv[3];
        S(1,5) = S(2,6) = S(3,8) = Sv[4];
        S(0,8) = S(2,7) = S(4,5) = Sv[5];

        // Computes K by Gauss-Legendre quadrature

        alglib::rmatrixenforcesymmetricity(S,9,1);
        matrix K1 = math::prod(shape.weight[i],B,S);
        matrix K2 = math::prod(detJ[i],K1,B,0,1);
        math::add(1,1,K2,K);
    }
    return K;
}

// --------------------------------------------------------------|
// Computes the non-linear equilibrium force for finite strain   |
// --------------------------------------------------------------|

darray Elem::selfFX(shapeStruct &shape,array3d LmR){

    matrix B;
    darray Fx;
    Fx.setlength(3*nLen);
    B.setlength(3*nLen,6);
    matrix D = math::stiffness(LmR);
    math::zero(Fx);
    math::zero(B);

    // Performs the numerical integration

    for(int i=0; i<shape.gLen; i++){
        for(int j=0; j<nLen; j++){

            // Computes the non-linear shape functions derivative matrix

            for(int k=0; k<3; k++){
            
                B(j+k*nLen,0) = F[i](k,0)*dN[0](j,i);
                B(j+k*nLen,1) = F[i](k,1)*dN[1](j,i);
                B(j+k*nLen,2) = F[i](k,2)*dN[2](j,i);

                B(j+k*nLen,3) = F[i](k,0)*dN[1](j,i)+F[i](k,1)*dN[0](j,i);
                B(j+k*nLen,4) = F[i](k,1)*dN[2](j,i)+F[i](k,2)*dN[1](j,i);
                B(j+k*nLen,5) = F[i](k,2)*dN[0](j,i)+F[i](k,0)*dN[2](j,i);
            }
        }

        // Computes K by Gauss-Legendre quadrature

        darray F1 = math::prod(shape.weight[i],D,E[i]);
        darray F2 = math::prod(detJ[i],B,F1);
        alglib::vadd(&Fx[0],&F2[0],3*nLen);
    }
    return Fx;
}

// --------------------------------------------------------|
// Computes the averaged Vin Mises stress in the element   |
// --------------------------------------------------------|

double Elem::stress(shapeStruct &shape,array3d LmR,darray u){

    matrix S;
    S.setlength(3,3);
    double VM = 0,volume=0;
    matrix D = math::stiffness(LmR);

    // Volume and stress field at Gauss points

    for(int i=0; i<shape.gLen; i++){

        darray Sv = math::prod(1,D,E[i]);

        // Biuilds the second Pila Kirchhoff stress tensor

        S(0,0) = Sv(0);
        S(1,1) = Sv(1);
        S(2,2) = Sv(2);
        S(0,1) = S(1,0) = Sv(3);
        S(1,2) = S(2,1) = Sv(4);
        S(2,0) = S(0,2) = Sv(5);

        // Builds the Cauchy stress tensor

        double detF = alglib::rmatrixdet(F[i],3);
        matrix C = math::prod(1/detF,F[i],S);
        C = math::prod(1,S,F[i],0,1);

        // Computes the square of the Von Mises stress

        double VMe = (pow(C(0,0)-C(1,1),2)+pow(C(1,1)-C(2,2),2)+pow(C(2,2)-C(0,0),2))/2;
        VMe += 3*(C[0][1]*C[0][1]+C[1][2]*C[1][2]+C[2][0]*C[2][0]);

        // Integration of the volume and the stress

        volume += shape.weight[i]*detJ[i];
        VM += shape.weight[i]*detJ[i]*sqrt(VMe);
    }

    VM /= volume;
    return VM;
}

// -----------------------------------------------------------|
// Class of quadrangle isoparametric lagrange element in 2D   |
// -----------------------------------------------------------|

Face::Face(vector<array3d> nXYZ,shapeStruct &shape){

    nLen = shape.N.cols();
    detJ2D.resize(shape.gLen,0);

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
        for(int i=0; i<3; i++){detJ2D[j] += dJ[i]*dJ[i];}
        detJ2D[j] = sqrt(detJ2D[j]);
    }
}

// --------------------------------------------------------------|
// Integrates the matrix of shape functions N over the element   |
// --------------------------------------------------------------|

darray Face::selfFT(shapeStruct &shape,darray F){

    matrix N;
    darray FT;
    FT.setlength(3*nLen);
    N.setlength(3*nLen,3);
    math::zero(FT);
    math::zero(N);

    // Integrates the matrix of local shape functions

    for(int i=0; i<shape.gLen; i++){
        for(int j=0; j<nLen; j++){
            N(j,0) = N(j+nLen,1) = N(j+2*nLen,2) = shape.N(j,i);
        }

        // Computes M by Gauss-Legendre quadrature

        double wdetJ = shape.weight[i]*detJ2D[i];
        darray Fe = math::prod(wdetJ,N,F,0);
        alglib::vadd(&FT[0],&Fe[0],3*nLen);
    }
    return FT;
}