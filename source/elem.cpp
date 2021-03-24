#include "..\include\elem.h"
using namespace std;

// -----------------------------------------------------------|
// Class of hexahedron isoparametric lagrange element in 3D   |
// -----------------------------------------------------------|

Elem::Elem(vector<array3d> inp1,ivector inp2) : nXYZ{inp1},surface{inp2}{

    // Sets the number of nodes in each dimension

    nLen = nXYZ.size();
    sLen = cbrt(nLen+1e-5);
}

// -----------------------------------------------------------|
// Evaluates the Jacobian and derivative of shape functions   |
// -----------------------------------------------------------|

void Elem::updateJ(shapeStruct &shape){

    matrix coord;
    J.resize(shape.gLen);
    detJ.resize(shape.gLen);
    gXYZ.resize(shape.gLen);
    coord.setlength(nLen,3);
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

    // Global coordinates of the Gauss points

    for(int i=0; i<nLen; i++){

        coord(i,0) = nXYZ[i][0];
        coord(i,1) = nXYZ[i][1];
        coord(i,2) = nXYZ[i][2];
    }

    matrix M = math::prod(1,shape.N,coord,1);
    for(int i=0; i<shape.gLen; i++){gXYZ[i] = {M(i,0),M(i,1),M(i,2)};}
}

// --------------------------------------------------------|
// This should normally free the memory of those vectors   |
// --------------------------------------------------------|

void Elem::freeJdN(){

    dvector().swap(detJ);
    vector<matrix>().swap(J);
    vector<array3d>().swap(gXYZ);
    for(int i=0; i<3; i++){dN[i].setlength(0,0);}
}

// ----------------------------------------------------------|
// Computes the elemental stiffness matrix K for local FEM   |
// ----------------------------------------------------------|

matrix Elem::selfK(shapeStruct &shape,array3d LmR){

    matrix B,K;
    B.setlength(3*nLen,6);
    K.setlength(3*nLen,3*nLen);
    math::zero(B);
    math::zero(K);

    // Update the Jacobian and builds the stiffness matrix

    updateJ(shape);
    matrix D = math::stiffness(LmR);

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
    updateJ(shape);
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

// -----------------------------------------------------------|
// Computes the elemental strain matrix S for non-local FEM   |
// -----------------------------------------------------------|

matrix Elem::selfS(shapeStruct &shape,array3d xyz,double range){

    matrix S;
    S.setlength(6,3*nLen);
    math::zero(S);

    // Performs the numerical integration

    for(int i=0; i<shape.gLen; i++){
        double k = math::kernel(xyz,gXYZ[i],range)*shape.weight[i]*detJ[i];

        for(int j=0; j<nLen; j++){

            S(0,j) = S(3,j+nLen) = S(5,j+2*nLen) += k*dN[0](j,i);
            S(3,j) = S(1,j+nLen) = S(4,j+2*nLen) += k*dN[1](j,i);
            S(5,j) = S(4,j+nLen) = S(2,j+2*nLen) += k*dN[2](j,i);
        }
    }
    return S;
}

// -------------------------------------------------------------------|
// Computes the elemental surface stiffness matrix Ks for local FEM   |
// -------------------------------------------------------------------|

pair<matrix,darray> Elem::selfKB(shapeStruct &shape,shapeStruct (&shapeS)[6],array3d LmS){

    darray BS;
    matrix B,K;
    ivector node;
    array3d v,norm;
    BS.setlength(3*nLen);
    B.setlength(3*nLen,6);
    K.setlength(3*nLen,3*nLen);
    math::zero(BS);
    math::zero(K);

    // Surface traction vector and stiffness matrix

    darray tau;
    matrix D = math::stiffness(LmS);
    double tauS[6] = {LmS[2],LmS[2],LmS[2],0,0,0};
    tau.setcontent(6,tauS);

    // Update the Jacobian and numerical integration

    for(int k:surface){

        math::zero(B);
        updateJ(shapeS[k]);
        dvector detJ2D = surfaceJ(shape,node,k);

        for(int i=0; i<shapeS[k].gLen; i++){
            for(int j:node){

                // Computes the shape functions derivative matrix B
                
                B(j,0) = B(j+nLen,3) = B(j+2*nLen,5) = dN[0](j,i);
                B(j,3) = B(j+nLen,1) = B(j+2*nLen,4) = dN[1](j,i);
                B(j,5) = B(j+nLen,4) = B(j+2*nLen,2) = dN[2](j,i);
            }

            // Computes the exterior normal to the surface

            array3d vr = {J[i](0,0),J[i](0,1),J[i](0,2)};
            array3d vs = {J[i](1,0),J[i](1,1),J[i](1,2)};
            array3d vt = {J[i](2,0),J[i](2,1),J[i](2,2)};

            switch(k){
            case 0:
                v = math::dotsub(vs,vr);
                norm = math::cross(vs,v,1);
                break;
            case 1:
                v = math::dotsub(vs,vr);
                norm = math::cross(v,vs,1);
                break;
            case 2:
                v = math::dotsub(vr,vt);
                norm = math::cross(v,vr,1);
                break;
            case 3:
                v = math::dotsub(vr,vt);
                norm = math::cross(vr,v,1);
                break;
            case 4:
                v = math::dotsub(vt,vs);
                norm = math::cross(vt,v,1);
                break;
            case 5:
                v = math::dotsub(vt,vs);
                norm = math::cross(v,vt,1);
                break;
            }

            // Computes the isotropic surface stiffness tensor

            matrix T = math::projection(norm);
            matrix S = math::prod(1,T,D);
            S = math::prod(1,S,T,0,1);

            // Computes K by Gauss-Legendre quadrature

            matrix BT = math::prod(1,B,T);
            matrix K1 = math::prod(shape.weight[i],BT,S);
            matrix K2 = math::prod(detJ2D[i],K1,BT,0,1);
            math::add(1,1,K2,K);

            // Computes B by Gauss-Legendre quadrature

            darray Be = math::prod(shape.weight[i],BT,tau);
            math::add(-detJ2D[i],1,Be,BS);
        }
    }

    pair<matrix,darray> Kb = {K,BS};
    return Kb;
}

// --------------------------------------------------------|
// Computes the averaged Von Mises stress in the element   |
// --------------------------------------------------------|

darray Elem::stress(shapeStruct &shape,array3d EvR,darray u){

    matrix B;
    darray sigma;
    double volume = 0;
    sigma.setlength(6);
    B.setlength(3*nLen,6);
    math::zero(sigma);
    math::zero(B);

    // Update the Jacobian and builds the stiffness matrix

    updateJ(shape);
    matrix D = math::stiffness(EvR);

    // Performs the numerical integration

    for(int i=0; i<shape.gLen; i++){
        for(int j=0; j<nLen; j++){

            // Computes the shape functions derivative matrix B
            
            B(j,0) = B(j+nLen,3) = B(j+2*nLen,5) = dN[0](j,i);
            B(j,3) = B(j+nLen,1) = B(j+2*nLen,4) = dN[1](j,i);
            B(j,5) = B(j+nLen,4) = B(j+2*nLen,2) = dN[2](j,i);
        }

        // Volume and stress field at Gauss points

        darray strain = math::prod(1,B,u,1);
        darray stress = math::prod(1,D,strain,0);
        math::add(shape.weight[i]*detJ[i],1,stress,sigma);
        volume += shape.weight[i]*detJ[i];
    }

    for(int i=0; i<6; i++){sigma[i] /= volume;}
    return sigma;
}

// -----------------------------------------------------------|
// Determinant of the 2D face Jacobian in the 3D hexahedron   |
// -----------------------------------------------------------|

dvector Elem::surfaceJ(shapeStruct &shape,ivector &node,int index){

    node.clear();
    vector<array3d> fXYZ;

    // Creates the node and coordinates of the face

    switch(index){
    case 0:
        for(int i=0; i<sLen*sLen; i++){
            fXYZ.push_back(nXYZ[i*sLen]);
            node.push_back(i*sLen);
        }
        break;
    case 1:
        for(int i=0; i<sLen*sLen; i++){
            node.push_back(i*sLen+sLen-1);
            fXYZ.push_back(nXYZ[i*sLen+sLen-1]);
        }
        break;
    case 2:
        for(int i=0; i<sLen; i++){
            for(int j=0; j<sLen; j++){
                node.push_back(i*sLen*sLen+j);
                fXYZ.push_back(nXYZ[i*sLen*sLen+j]);
            }
        }
        break;
    case 3:
        for(int i=0; i<sLen; i++){
            for(int j=0; j<sLen; j++){
                node.push_back(i*sLen*sLen+j+(sLen-1)*sLen);
                fXYZ.push_back(nXYZ[i*sLen*sLen+j+(sLen-1)*sLen]);
            }
        }
        break;
    case 4:
        for(int i=0; i<sLen*sLen; i++){
            fXYZ.push_back(nXYZ[i]);
            node.push_back(i);
        }
        break;
    case 5:
        for(int i=0; i<sLen*sLen; i++){
            node.push_back(i+(sLen-1)*sLen*sLen);
            fXYZ.push_back(nXYZ[i+(sLen-1)*sLen*sLen]);
        }
        break;
    }

    // Creates the 2D face element and its the determinant

    Face face(fXYZ,shape);
    return face.detJ2D;
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

darray Face::selfB(shapeStruct &shape,darray F){

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

        double wdetJ = shape.weight[i]*detJ2D[i];
        darray Be = math::prod(wdetJ,N,F,0);
        math::add(1,1,Be,B);
    }
    return B;
}