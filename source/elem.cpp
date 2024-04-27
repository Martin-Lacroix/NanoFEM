#include "elem.h"
using namespace std;

// -----------------------------------------------------------|
// Class of hexahedron isoparametric Lagrange element in 3D   |
// -----------------------------------------------------------|

Elem::Elem(vector<array3d> inp1,ivector inp2) : nXYZ{inp1},surface{inp2}{
    nLen = nXYZ.size();
}

// ---------------------------------------------------------|
// Jacobian and derivative of shape functions in the bulk   |
// ---------------------------------------------------------|

void Elem::updateJ(shapeStruct &shape){

    detJ.resize(shape.gLen);
    vector<matrix3d> J(shape.gLen);
    vector<matrix3d> invJ(shape.gLen);

    // Resets the shape function derivative

    for(int i = 0; i < 3; i++){
        
        dN[i].resize(nLen,shape.gLen);
        dN[i].setZero();
    }

    // Builds the Jacobian matrix and global coordinates

    for(int i = 0; i < shape.gLen; i++){

        J[i].setZero();

        // Computes the Jacobian matrix at this Gauss point

        for(int j = 0; j < nLen; j++){
            for(int k = 0; k < 3; k++){
                for(int n = 0; n < 3; n++){
                    J[i](n,k) += shape.dN[n](j,i)*nXYZ[j][k];
                }
            }
        }
    }

    // Determinant and inverse of the Jacobian matrix

    for(int i = 0; i < shape.gLen; i++){

        invJ[i] = J[i].inverse();
        detJ[i] = J[i].determinant();

        // Global derivatives of shape functions

        for(int j = 0; j < nLen; j++){
            for(int k = 0; k < 3; k++){
                for(int n = 0; n < 3; n++){
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

    for(int s : surface){

        darray3d v;
        norm[s].resize(shape[s].gLen);
        detJ2D[s].resize(shape[s].gLen);
        vector<matrix3d> J(shape[s].gLen);
        vector<matrix3d> invJ(shape[s].gLen);

        // Resets the shape function derivative

        for(int i = 0; i < 3; i++){
            
            dNs[s][i].resize(nLen,shape[s].gLen);
            dNs[s][i].setZero();
        }

        // Builds the Jacobian matrix and global coordinates

        for(int i = 0; i < shape[s].gLen; i++){

            J[i].setZero();

            // Computes the Jacobian matrix at this Gauss point

            for(int j = 0; j < nLen; j++){
                for(int k = 0; k < 3; k++){
                    for(int n = 0; n < 3; n++){
                        J[i](n,k) += shape[s].dN[n](j,i)*nXYZ[j][k];
                    }
                }
            }

            // Normal to the s-th surface of the element

            darray3d vr = {J[i](0,0),J[i](0,1),J[i](0,2)};
            darray3d vs = {J[i](1,0),J[i](1,1),J[i](1,2)};
            darray3d vt = {J[i](2,0),J[i](2,1),J[i](2,2)};

            switch(s){
            case 0:
                v = vs-vs.dot(vr)*vr;
                norm[s][i] = vs.cross(v);
                break;
            case 1:
                v = vs-vs.dot(vr)*vr;
                norm[s][i] = v.cross(vs);
                break;
            case 2:
                v = vr-vr.dot(vt)*vt;
                norm[s][i] = v.cross(vr);
                break;
            case 3:
                v = vr-vr.dot(vt)*vt;
                norm[s][i] = vr.cross(v);
                break;
            case 4:
                v = vt-vt.dot(vs)*vs;
                norm[s][i] = vt.cross(v);
                break;
            case 5:
                v = vt-vt.dot(vs)*vs;
                norm[s][i] = v.cross(vt);
                break;
            }
        }

        // Determinant and inverse of the Jacobian matrix

        for(int i = 0; i < shape[s].gLen; i++){
            
            invJ[i] = J[i].inverse();
            detJ2D[s][i] = norm[s][i].norm();
            for(int j = 0; j < 3; j++) norm[s][i][j] /= detJ2D[s][i];

            // Global derivatives of shape functions

            for(int j = 0; j < nLen; j++){
                for(int k = 0; k < 3; k++){
                    for(int n = 0; n < 3; n++){
                        dNs[s][k](j,i) += shape[s].dN[n](j,i)*invJ[i](n,k);
                    }
                }
            }
        }
    }
}

// -----------------------------------------------------|
// Computes the elemental deformation gradient tensor   |
// -----------------------------------------------------|

void Elem::updateF(shapeStruct &shape,darray u){

    matrix B(3*nLen,9);
    B.setZero();

    // Update the Jacobian and the deformation gradient

    F.resize(shape.gLen);
    E.resize(shape.gLen);

    for(int i = 0; i < shape.gLen; i++){
        for(int j = 0; j < nLen; j++){

            // Computes the linear shape functions derivative matrix
            
            B(j,0) = B(j+nLen,4) = B(j+2*nLen,7) = dN[0](j,i);
            B(j,3) = B(j+nLen,1) = B(j+2*nLen,6) = dN[1](j,i);
            B(j,8) = B(j+nLen,5) = B(j+2*nLen,2) = dN[2](j,i);
        }

        // Stores the deformation gradient tensor

        darray Fv =  B.transpose()*u;

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

        E[i].resize(6);
        matrix Ev = F[i].transpose()*F[i];

        E[i](0) = (Ev(0,0)-1)/2;
        E[i](1) = (Ev(1,1)-1)/2;
        E[i](2) = (Ev(2,2)-1)/2;
        E[i](3) = Ev(0,1)/2;
        E[i](4) = Ev(1,2)/2;
        E[i](5) = Ev(2,0)/2;
    }
}

// --------------------------------------------------------|
// This should normally free the memory of those vectors   |
// --------------------------------------------------------|

void Elem::clean(){

    // Clear bulk Jacobian and parameters

    dvector().swap(detJ);
    vector<darray>().swap(E);
    vector<matrix3d>().swap(F);
    for(int i = 0; i < 3; i++) dN[i].resize(0,0);

    // Clear surface Jacobian and parameters

    for(int s : surface){
        
        dvector().swap(detJ2D[s]);
        vector<darray3d>().swap(norm[s]);
        for(int i = 0; i < 3; i++) dNs[s][i].resize(0,0);
    }
}

// ----------------------------------------------------------|
// Computes the elemental stiffness matrix K for local FEM   |
// ----------------------------------------------------------|

matrix Elem::selfK(shapeStruct &shape,array3d LmR){

    matrix B(3*nLen,6);
    matrix K(3*nLen,3*nLen);
    matrix D = math::stiffness(LmR);
    B.setZero();
    K.setZero();

    // Performs the numerical integration

    for(int i = 0; i < shape.gLen; i++){
        for(int j = 0; j < nLen; j++){

            // Computes the shape functions derivative matrix B
            
            B(j,0) = B(j+nLen,3) = B(j+2*nLen,5) = dN[0](j,i);
            B(j,3) = B(j+nLen,1) = B(j+2*nLen,4) = dN[1](j,i);
            B(j,5) = B(j+nLen,4) = B(j+2*nLen,2) = dN[2](j,i);
        }

        // Computes K by Gauss-Legendre quadrature

        double wdetJ = shape.weight[i]*detJ[i];
        K += wdetJ*B*D*B.transpose();
    }
    return K;
}

// -----------------------------------------------------|
// Computes the elemental mass matrix M for local FEM   |
// -----------------------------------------------------|

matrix Elem::selfM(shapeStruct &shape,double rho){

    matrix N(3*nLen,3);
    matrix M(3*nLen,3*nLen);
    N.setZero();
    M.setZero();

    // Performs the numerical integration

    for(int i = 0; i < shape.gLen; i++){
        for(int j = 0; j < nLen; j++){
            N(j,0) = N(j+nLen,1) = N(j+2*nLen,2) = shape.N(j,i);
        }

        // Computes M by Gauss-Legendre quadrature

        double wdetJ = shape.weight[i]*detJ[i];
        M += wdetJ*rho*N*N.transpose();
    }
    return M;
}

// ------------------------------------------------------|
// Computes the surface stiffness matrix Ks for SST/SET  |
// ------------------------------------------------------|

matrix Elem::selfKS(shapeStruct (&shape)[6],array3d LmS){

    matrix B(3*nLen,6);
    matrix K(3*nLen,3*nLen);
    matrix D = math::stiffness(LmS);
    K.setZero();

    // Update the Jacobian and numerical integration

    for(int k : surface){
        B.setZero();

        for(int i = 0; i < shape[k].gLen; i++){
            for(int j = 0; j < nLen; j++){

                // Computes the shape functions derivative matrix B
                
                B(j,0) = B(j+nLen,3) = B(j+2*nLen,5) = dNs[k][0](j,i);
                B(j,3) = B(j+nLen,1) = B(j+2*nLen,4) = dNs[k][1](j,i);
                B(j,5) = B(j+nLen,4) = B(j+2*nLen,2) = dNs[k][2](j,i);
            }

            // Computes the isotropic surface stiffness tensor

            matrix T = math::projection(norm[k][i]);
            matrix S = T*D*T.transpose();
            matrix BT = B*T;

            // Computes K by Gauss-Legendre quadrature

            double wdetJ = shape[k].weight[i]*detJ2D[k][i];
            K += wdetJ*BT*S*BT.transpose();
        }
    }
    return K;
}

// ------------------------------------------------------|
// Computes the surface boundary vector Fs for SST/SET   |
// ------------------------------------------------------|

darray Elem::selfFS(shapeStruct (&shape)[6],array3d LmS){

    matrix B(3*nLen,6);
    darray F(3*nLen);
    F.setZero();

    // Surface traction vector and stiffness matrix

    darray tau(6);
    tau << LmS[2],LmS[2],LmS[2],0,0,0;

    // Update the Jacobian and numerical integration

    for(int s : surface){
        B.setZero();

        for(int i = 0; i < shape[s].gLen; i++){
            for(int j = 0; j < nLen; j++){

                // Computes the shape functions derivative matrix B
                
                B(j,0) = B(j+nLen,3) = B(j+2*nLen,5) = dNs[s][0](j,i);
                B(j,3) = B(j+nLen,1) = B(j+2*nLen,4) = dNs[s][1](j,i);
                B(j,5) = B(j+nLen,4) = B(j+2*nLen,2) = dNs[s][2](j,i);
            }

            // Computes the surface tension by Gauss-Legendre quadrature

            double wdetJ = shape[s].weight[i]*detJ2D[s][i];
            matrix T = math::projection(norm[s][i]);
            F -= wdetJ*B*T*tau;
        }
    }
    return F;
}

// ---------------------------------------------------|
// Computes the non-linear stiffness matrix for SVK   |
// ---------------------------------------------------|

matrix Elem::selfKN(shapeStruct &shape,array3d LmR){

    matrix B(3*nLen,6);
    matrix K(3*nLen,3*nLen);
    matrix D = math::stiffness(LmR);
    K.setZero();
    B.setZero();

    // Performs the numerical integration

    for(int i = 0; i < shape.gLen; i++){
        for(int j = 0; j < nLen; j++){

            // Computes the non-linear shape functions derivative matrix

            for(int k = 0; k < 3; k++){
            
                B(j+k*nLen,0) = F[i](k,0)*dN[0](j,i);
                B(j+k*nLen,1) = F[i](k,1)*dN[1](j,i);
                B(j+k*nLen,2) = F[i](k,2)*dN[2](j,i);

                B(j+k*nLen,3) = F[i](k,0)*dN[1](j,i)+F[i](k,1)*dN[0](j,i);
                B(j+k*nLen,4) = F[i](k,1)*dN[2](j,i)+F[i](k,2)*dN[1](j,i);
                B(j+k*nLen,5) = F[i](k,2)*dN[0](j,i)+F[i](k,0)*dN[2](j,i);
            }
        }

        // Computes K by Gauss-Legendre quadrature
        
        double wdetJ = shape.weight[i]*detJ[i];
        K += wdetJ*B*D*B.transpose();
    }
    return K;
}

// -----------------------------------------------|
// Computes the linear stiffness matrix for SVK   |
// -----------------------------------------------|

matrix Elem::selfKL(shapeStruct &shape,array3d LmR){

    matrix S(9,9);
    matrix B(3*nLen,9);
    matrix K(3*nLen,3*nLen);
    matrix D = math::stiffness(LmR);
    K.setZero();
    B.setZero();
    S.setZero();

    // Performs the numerical integration

    for(int i = 0; i < shape.gLen; i++){
        for(int j = 0; j < nLen; j++){

            // Computes the linear shape functions derivative matrix
            
            B(j,0) = B(j+nLen,4) = B(j+2*nLen,7) = dN[0](j,i);
            B(j,3) = B(j+nLen,1) = B(j+2*nLen,6) = dN[1](j,i);
            B(j,8) = B(j+nLen,5) = B(j+2*nLen,2) = dN[2](j,i);
        }

        // Computes the matrix of second Piola-Kirchhoff components

        darray Sv = D*E[i];

        S(0,0) = S(4,4) = S(7,7) = Sv[0];
        S(1,1) = S(3,3) = S(6,6) = Sv[1];
        S(2,2) = S(5,5) = S(8,8) = Sv[2];

        S(0,3) = S(1,4) = S(6,7) = Sv[3];
        S(1,5) = S(2,6) = S(3,8) = Sv[4];
        S(0,8) = S(2,7) = S(4,5) = Sv[5];

        S(3,0) = S(4,1) = S(7,6) = Sv[3];
        S(5,1) = S(6,2) = S(8,3) = Sv[4];
        S(8,0) = S(7,2) = S(5,4) = Sv[5];

        // Computes K by Gauss-Legendre quadrature

        double wdetJ = shape.weight[i]*detJ[i];
        K += wdetJ*B*S*B.transpose();
    }
    return K;
}

// ----------------------------------------------------|
// Computes the non-linear equilibrium force for SVK   |
// ----------------------------------------------------|

darray Elem::selfFX(shapeStruct &shape,array3d LmR){

    darray Fx(3*nLen);
    matrix B(3*nLen,6);
    matrix D = math::stiffness(LmR);
    Fx.setZero();
    B.setZero();

    // Performs the numerical integration

    for(int i = 0; i < shape.gLen; i++){
        for(int j = 0; j < nLen; j++){

            // Computes the non-linear shape functions derivative matrix

            for(int k = 0; k < 3; k++){
            
                B(j+k*nLen,0) = F[i](k,0)*dN[0](j,i);
                B(j+k*nLen,1) = F[i](k,1)*dN[1](j,i);
                B(j+k*nLen,2) = F[i](k,2)*dN[2](j,i);

                B(j+k*nLen,3) = F[i](k,0)*dN[1](j,i)+F[i](k,1)*dN[0](j,i);
                B(j+k*nLen,4) = F[i](k,1)*dN[2](j,i)+F[i](k,2)*dN[1](j,i);
                B(j+k*nLen,5) = F[i](k,2)*dN[0](j,i)+F[i](k,0)*dN[2](j,i);
            }
        }

        // Computes K by Gauss-Legendre quadrature

        double wdetJ = shape.weight[i]*detJ[i];
        Fx += wdetJ*B*D*E[i];
    }
    return Fx;
}

// --------------------------------------------------|
// Averaged Von Mises stress in the element in SVK   |
// --------------------------------------------------|

double Elem::stress(shapeStruct &shape,array3d LmR){

    matrix3d S;
    double VM = 0,volume=0;
    matrix D = math::stiffness(LmR);

    // Volume and stress field at Gauss points

    for(int i = 0; i < shape.gLen; i++){

        darray Sv = D*E[i];

        // Biuilds the second Pila Kirchhoff stress tensor

        S(0,0) = Sv(0);
        S(1,1) = Sv(1);
        S(2,2) = Sv(2);
        S(0,1) = S(1,0) = Sv(3);
        S(1,2) = S(2,1) = Sv(4);
        S(2,0) = S(0,2) = Sv(5);

        // Builds the Cauchy stress tensor

        double detF = F[i].determinant();
        S = F[i]*S*F[i].transpose()/detF;

        // Computes the square of the Von Mises stress

        double VMe = (pow(S(0,0)-S(1,1),2)+pow(S(1,1)-S(2,2),2)+pow(S(2,2)-S(0,0),2))/2;
        VMe += 3*(S(0,1)*S(0,1)+S(1,2)*S(1,2)+S(2,0)*S(2,0));

        // Integration of the volume and the stress

        volume += shape.weight[i]*detJ[i];
        VM += shape.weight[i]*detJ[i]*sqrt(VMe);
    }

    VM /= volume;
    return VM;
}

// ------------------------------------------------------|
// Averaged Von Mises stress in the element in SST/SET   |
// ------------------------------------------------------|

double Elem::stress(shapeStruct &shape,array3d LmR,darray u){

    matrix B(3*nLen,6);
    double VM = 0,volume=0;
    matrix D = math::stiffness(LmR);
    B.setZero();

    // Performs the numerical integration

    for(int i = 0; i < shape.gLen; i++){
        for(int j = 0; j < nLen; j++){

            // Computes the shape functions derivative matrix B
            
            B(j,0) = B(j+nLen,3) = B(j+2*nLen,5) = dN[0](j,i);
            B(j,3) = B(j+nLen,1) = B(j+2*nLen,4) = dN[1](j,i);
            B(j,5) = B(j+nLen,4) = B(j+2*nLen,2) = dN[2](j,i);
        }
        
        // Builds the Cauchy stress tensor

        darray S = B.transpose()*u;
        S = D*S;

        // Computes the square of the Von Mises stress

        double VMe = (pow(S[0]-S[1],2)+pow(S[1]-S[2],2)+pow(S[2]-S[0],2))/2;
        VMe += 3*(S[3]*S[3]+S[4]*S[4]+S[5]*S[5]);

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

    for(int j = 0; j < shape.gLen; j++){

        darray3d J2Dr = darray3d::Zero();
        darray3d J2Ds = darray3d::Zero();

        // Computes the partial Jacobian vectors for 3D to 2D mapping

        for(int k = 0; k < 3; k++){
            for(int i = 0; i < nLen; i++){
            
                J2Dr[k] += shape.dN[0](i,j)*nXYZ[i][k];
                J2Ds[k] += shape.dN[1](i,j)*nXYZ[i][k];
            }
        }

        // Computes the cross product and takes the norm

        darray3d dJ = J2Dr.cross(J2Ds);
        for(int i = 0; i < 3; i++) detJ2D[j] += dJ[i]*dJ[i];
        detJ2D[j] = sqrt(detJ2D[j]);
    }
}

// --------------------------------------------------------------|
// Integrates the matrix of shape functions N over the element   |
// --------------------------------------------------------------|

darray Face::selfFT(shapeStruct &shape,darray F){

    darray FT(3*nLen);
    matrix N(3*nLen,3);
    FT.setZero();
    N.setZero();

    // Integrates the matrix of local shape functions

    for(int i = 0; i < shape.gLen; i++){
        for(int j = 0; j < nLen; j++){
            N(j,0) = N(j+nLen,1) = N(j+2*nLen,2) = shape.N(j,i);
        }

        // Computes M by Gauss-Legendre quadrature

        double wdetJ = shape.weight[i]*detJ2D[i];
        FT += wdetJ*N*F;
    }
    return FT;
}