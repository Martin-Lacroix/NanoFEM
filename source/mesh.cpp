#include "..\include\mesh.h"
using namespace std;

// Builds the elements, shape functions and quadratures

Mesh::Mesh(meshStruct mesh,otherStruct param){
    
    meshData = mesh;
    otherData = param;
    int fLen = mesh.fNode.size();
    int eLen = mesh.eNode.size();

    eQuad = math::legendre(3,param.order);
    fQuad = math::legendre(2,param.order);
    eShape = shape(3);
    fShape = shape(2);

    // Stores the mesh elements into a vector

    for(int i=0; i<eLen; i++){

        vector<darray> eXYZ;
        int nLen = mesh.eNode[i].length();
        for(int j=0; j<nLen; j++){eXYZ.push_back(mesh.nXYZ[mesh.eNode[i][j]]);}
        Elem elem(eXYZ,eShape);
        eList.push_back(elem);
    }

    // Stores the mesh faces into a vector

    for(int i=0; i<fLen; i++){

        vector<darray> fXYZ;
        int nLen = mesh.fNode[i].length();
        for(int j=0; j<nLen; j++){fXYZ.push_back(mesh.nXYZ[mesh.fNode[i][j]]);}
        Face face(fXYZ,fShape);
        fList.push_back(face);
    }
}

// Structure of linear shape functions

shapeStruct Mesh::shape(int dim){

    int nLen;
    int gLen;
    shapeStruct shape;
    if(dim==3){gLen = eQuad.gRST.size(); nLen = 8;}
    if(dim==2){gLen = fQuad.gRST.size(); nLen = 4;}
    vector<vector<int>> n(nLen,vector<int>(dim));
    shape.gLen = gLen;
    shape.nLen = nLen;

    // Memory allocation

    shape.drN.setlength(nLen,gLen);
    shape.dsN.setlength(nLen,gLen);
    shape.dtN.setlength(nLen,gLen);
    shape.N.setlength(nLen,gLen);

    // 3D shape functions constructor

    if(dim==3){

        n[0]={-1,-1,-1}; n[1]={1,-1,-1}; n[2]={1,1,-1}; n[3]={-1,1,-1};
        n[4]={-1,-1,1}; n[5]={1,-1,1}; n[6]={1,1,1}; n[7]={-1,1,1};

        // Shape functions at Gauss points

        for(int i=0; i<gLen; i++){

            double r = eQuad.gRST[i][0];
            double s = eQuad.gRST[i][1];
            double t = eQuad.gRST[i][2];

            for(int j=0; j<nLen; j++){

                // Shape functions at Gauss points

                shape.N(j,i) = (1+n[j][0]*r)*(1+n[j][1]*s)*(1+n[j][2]*t)/8;
                shape.drN(j,i) = n[j][0]*(1+n[j][1]*s)*(1+n[j][2]*t)/8;
                shape.dsN(j,i) = n[j][1]*(1+n[j][0]*r)*(1+n[j][2]*t)/8;
                shape.dtN(j,i) = n[j][2]*(1+n[j][0]*r)*(1+n[j][1]*s)/8;
            }
        }
    }

    // 2D shape functions constructor

    if(dim==2){

        n[0]={-1,-1}; n[1]={1,-1};
        n[2]={1,1}; n[3]={-1,1};

        // Shape functions at Gauss points

        for(int i=0; i<gLen; i++){

            double r = fQuad.gRST[i][0];
            double s = fQuad.gRST[i][1];

            for(int j=0; j<nLen; j++){

                // Shape functions at Gauss points

                shape.N(j,i) = (1+n[j][0]*r)*(1+n[j][1]*s)/4;
                shape.drN(j,i) = n[j][0]*(1+n[j][1]*s)/4;
                shape.dsN(j,i) = n[j][1]*(1+n[j][0]*r)/4;
            }
        }
    }
    return shape;
}

// Computes the local stiffness matrix K

sparse Mesh::localK(){

    sparse K;
    double eLen = eList.size();
    double nLen = meshData.nXYZ.size();
    alglib::sparsecreate(3*nLen,3*nLen,K);

    for(int i=0; i<eLen; i++){

        int nbr = eList[i].nLen;
        iarray eNode = meshData.eNode[i];
        matrix K1 = eList[i].selfK(eQuad,otherData.D);

        // Inserts the submatrices into the global matrice

        for(int j=0; j<nbr; j++){
            for(int k=0; k<nbr; k++){
                for(int m=0; m<3; m++){
                    for(int n=0; n<3; n++){

                        double val = K1[j+m*nbr][k+n*nbr];
                        alglib::sparseadd(K,eNode[j]+m*nLen,eNode[k]+n*nLen,val);
                    }
                }
            }
        }
    }
    return K;
}

// Computes the non-local stiffness matrix K

sparse Mesh::nonLocalK(){

    sparse K;
    double eLen = eList.size();
    double nLen = meshData.nXYZ.size();
    alglib::sparsecreate(3*nLen,3*nLen,K);

    for(int i=0; i<eLen; i++){

        matrix K1 = elemK(i);
        int nbr = eList[i].nLen;
        iarray eNode = meshData.eNode[i];

        // Inserts the submatrices into the global matrice

        for(int k=0; k<3*nLen; k++){
            for(int j=0; j<nbr; j++){
                for(int m=0; m<3; m++){

                    double val = K1[j+m*nbr][k];
                    alglib::sparseadd(K,eNode[j]+m*nLen,k,val);
                }
            }
        }
    }
    return K;
}

// Evaluates the non-local S matrix at a point xyz

matrix Mesh::totalS(darray xyz){

    matrix S;
    double eLen = eList.size();
    double nLen = meshData.nXYZ.size();
    S.setlength(6,3*nLen);
    math::zero(S);

    for(int i=0; i<eLen; i++){

        int nbr = eList[i].nLen;
        matrix S1 = eList[i].selfS(eQuad,xyz);

        for(int j=0; j<6; j++){
            for(int k=0; k<nbr; k++){
                for(int n=0; n<3; n++){
                    S(j,meshData.eNode[i][k]+n*nLen) += S1(j,k+n*nbr);
                }
            }
        }
    }
    return S;
}

// Computes the elemental non-local stiffness matrix

matrix Mesh::elemK(int idx){

    matrix B;
    matrix K;
    Elem& elem = eList[idx];
    double nLen = elem.nLen;
    double gLen = eQuad.gRST.size();
    K.setlength(3*nLen,3*meshData.nXYZ.size());
    B.setlength(3*nLen,6);
    math::zero(B);
    math::zero(K);

    // Performs the numerical integration

    for(int i=0; i<gLen; i++){
        for(int j=0; j<nLen; j++){
            
            B(j,0) = B(j+nLen,3) = B(j+2*nLen,4) = elem.dxN(j,i);
            B(j,3) = B(j+nLen,1) = B(j+2*nLen,5) = elem.dyN(j,i);
            B(j,4) = B(j+nLen,5) = B(j+2*nLen,2) = elem.dzN(j,i);
        }
        matrix S = totalS(elem.gXYZ[i]);
        matrix K1 = math::prod(eQuad.weight[i],B,otherData.D);
        matrix K2 = math::prod(elem.detJ[i],K1,S);
        math::add(1,1,K2,K);
    }
    return K;
}

// Builds the vector of Neumann boundary conditions

darray Mesh::neumann(){

    darray B;
    double fLen = fList.size();
    double nLen = meshData.nXYZ.size();
    B.setlength(3*nLen);
    math::zero(B);

    for(int i=0; i<fLen; i++){

        int nbr = fList[i].nLen;
        matrix N = fList[i].selfN(fShape,fQuad);
        darray B1 = math::prod(1,N,otherData.neumann[i]);

        // Inserts the subvectors into the global vector

        for(int j=0; j<nbr; j++){
            for(int k=0; k<3; k++){
                
                double val = B1(j+k*nbr);
                B(meshData.fNode[i][j]+k*nLen) += val;
            }
        }
    }
    return B;
}

// Prepares the vector for Dirichlet BC

void Mesh::dirichlet(darray &B){

    int nLen = meshData.nXYZ.size();
    vector<vector<int>> dirichlet = otherData.dirichlet;

    for(int i=0; i<3; i++){
        for(int j=0; j<dirichlet[i].size(); j++){
            B(dirichlet[i][j]+i*nLen) = 0;
        }
    }
}

// Prepares the matrix for Dirichlet BC

void Mesh::dirichlet(sparse &K){

    int nLen = meshData.nXYZ.size();
    vector<vector<int>> dirichlet = otherData.dirichlet;

    for(int i=0; i<3; i++){
        for(int j=0; j<dirichlet[i].size(); j++){

            // Zeroes the row except the element in the diagonal

            int idx = dirichlet[i][j]+i*nLen;
            alglib::sparseset(K,idx,idx,1);

            for(int n=0; n<3*nLen; n++){
                if(n!=idx){alglib::sparseset(K,idx,n,0);}
                if(n!=idx){alglib::sparseset(K,n,idx,0);}
            }
        }
    }
}