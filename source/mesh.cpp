#include "..\include\mesh.h"
using namespace std;

// Builds the elements, shape functions and quadratures

Mesh::Mesh(meshStruct mesh1,paramStruct param1){
    
    mesh = mesh1;
    param = param1;
    int fLen = mesh.fNode.size();
    int eLen = mesh.eNode.size();

    quad3D = math::legendre(3,param.order);
    quad2D = math::legendre(2,param.order);
    shape3D = shape(3);
    shape2D = shape(2);

    // Stores the mesh elements into a vector

    for(int i=0; i<eLen; i++){

        vector<darray> eXYZ;
        int nLen = mesh.eNode[i].length();
        for(int j=0; j<nLen; j++){eXYZ.push_back(mesh.nXYZ[mesh.eNode[i][j]]);}
        Elem elem(eXYZ,shape3D);
        eList.push_back(elem);
    }

    // Stores the mesh faces into a vector

    for(int i=0; i<fLen; i++){

        vector<darray> fXYZ;
        int nLen = mesh.fNode[i].length();
        for(int j=0; j<nLen; j++){fXYZ.push_back(mesh.nXYZ[mesh.fNode[i][j]]);}
        Face face(fXYZ,shape2D);
        fList.push_back(face);
    }
}

// Structure of linear shape functions

shapeStruct Mesh::shape(int dim){

    int nLen;
    int gLen;
    shapeStruct shape;
    if(dim==3){gLen = quad3D.gRST.size(); nLen = 8;}
    if(dim==2){gLen = quad2D.gRST.size(); nLen = 4;}
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

            double r = quad3D.gRST[i][0];
            double s = quad3D.gRST[i][1];
            double t = quad3D.gRST[i][2];

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

            double r = quad2D.gRST[i][0];
            double s = quad2D.gRST[i][1];

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
    double nLen = mesh.nXYZ.size();
    alglib::sparsecreate(3*nLen,3*nLen,K);

    for(int i=0; i<eLen; i++){

        int nbr = eList[i].nLen;
        iarray eNode = mesh.eNode[i];
        matrix K1 = eList[i].selfK(quad3D,param.D);

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
    double nLen = mesh.nXYZ.size();
    alglib::sparsecreate(3*nLen,3*nLen,K);

    for(int i=0; i<eLen; i++){

        matrix K1 = elemK(i);
        int nbr = eList[i].nLen;
        iarray eNode = mesh.eNode[i];

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
    double nLen = mesh.nXYZ.size();
    S.setlength(6,3*nLen);
    math::zero(S);

    for(int i=0; i<eLen; i++){

        int nbr = eList[i].nLen;
        matrix S1 = eList[i].selfS(quad3D,xyz);

        for(int j=0; j<6; j++){
            for(int k=0; k<nbr; k++){
                for(int n=0; n<3; n++){
                    S(j,mesh.eNode[i][k]+n*nLen) += S1(j,k+n*nbr);
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
    double gLen = quad3D.gRST.size();
    K.setlength(3*nLen,3*mesh.nXYZ.size());
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
        matrix K1 = math::prod(quad3D.weight[i],B,param.D);
        matrix K2 = math::prod(elem.detJ[i],K1,S);
        math::add(1,1,K2,K);
    }
    return K;
}

// Builds the vector of Neumann boundary conditions

darray Mesh::neumann(){

    darray B;
    double fLen = fList.size();
    double nLen = mesh.nXYZ.size();
    B.setlength(3*nLen);
    math::zero(B);

    for(int i=0; i<fLen; i++){

        int nbr = fList[i].nLen;
        matrix N = fList[i].selfN(shape2D,quad2D);
        darray B1 = math::prod(1,N,param.neumann[i]);

        // Inserts the subvectors into the global vector

        for(int j=0; j<nbr; j++){
            for(int k=0; k<3; k++){
                
                double val = B1(j+k*nbr);
                B(mesh.fNode[i][j]+k*nLen) += val;
            }
        }
    }
    return B;
}

// Prepares the vector for Dirichlet BC

void Mesh::dirichlet(darray &B){

    int nLen = mesh.nXYZ.size();
    vector<vector<int>> dirichlet = param.dirichlet;

    for(int i=0; i<3; i++){
        for(int j=0; j<dirichlet[i].size(); j++){
            B(dirichlet[i][j]+i*nLen) = 0;
        }
    }
}

// Cleans and prepares the matrix for Dirichlet BC

void Mesh::dirichlet(sparse &K){

    double val;
    alglib::ae_int_t i=0;
    alglib::ae_int_t j=0;
    alglib::ae_int_t end=0;
    alglib::ae_int_t start=0;

    int nLen = mesh.nXYZ.size();
    vector<vector<int>> row(3*nLen);
    vector<vector<int>> col(3*nLen);

    while(alglib::sparseenumerate(K,start,end,i,j,val)){

        if(abs(val)<1e-9){alglib::sparseset(K,i,j,0);}
        else{row[i].push_back(j);col[j].push_back(i);}
    }

    for(int n=0; n<3; n++){
        int dLen = param.dirichlet[n].size();

        for(int m=0; m<dLen; m++){
            int idx = param.dirichlet[n][m]+n*nLen;

            for(int k=0; k<row[idx].size(); k++){

                if(idx==row[idx][k]){alglib::sparseset(K,idx,row[idx][k],1);}
                else{alglib::sparseset(K,idx,row[idx][k],0);}
            }
            for(int k=0; k<col[idx].size(); k++){
                
                if(col[idx][k]==idx){alglib::sparseset(K,col[idx][k],idx,1);}
                else{alglib::sparseset(K,col[idx][k],idx,0);}
            }
        }
    }
}