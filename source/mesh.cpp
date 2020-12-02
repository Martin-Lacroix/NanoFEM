#include "..\include\mesh.h"
using namespace std;

// Builds the elements, shape functions and quadratures

Mesh::Mesh(meshStruct &input){
    
    mesh = input;
    int fLen = mesh.fNode.size();
    int eLen = mesh.eNode.size();

    quad3D = math::legendre(3,mesh.order);
    quad2D = math::legendre(2,mesh.order);
    shape3D = shape(3);
    shape2D = shape(2);

    // Stores the mesh elements into a vector

    for(int i=0; i<eLen; i++){

        vector<dvector> eXYZ;
        int nLen = mesh.eNode[i].size();
        for(int j=0; j<nLen; j++){eXYZ.push_back(mesh.nXYZ[mesh.eNode[i][j]]);}
        Elem elem(eXYZ,shape3D);
        eList.push_back(elem);
    }

    // Stores the mesh faces into a vector

    for(int i=0; i<fLen; i++){

        vector<dvector> fXYZ;
        int nLen = mesh.fNode[i].size();
        for(int j=0; j<nLen; j++){fXYZ.push_back(mesh.nXYZ[mesh.fNode[i][j]]);}
        Face face(fXYZ,shape2D);
        fList.push_back(face);
    }
}

// Structure of linear shape functions

shapeStruct Mesh::shape(int dim){

    int nLen,gLen;
    shapeStruct shape;
    if(dim==3){gLen = quad3D.gRST.size(); nLen = 8;}
    if(dim==2){gLen = quad2D.gRST.size(); nLen = 4;}
    vector<ivector> n(nLen,ivector(dim));
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
        ivector eNode = mesh.eNode[i];
        matrix K1 = eList[i].selfK(quad3D,mesh.D);

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
        ivector eNode = mesh.eNode[i];

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

matrix Mesh::totalS(dvector xyz){

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

    matrix B,K;
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
        matrix K1 = math::prod(quad3D.weight[i],B,mesh.D);
        matrix K2 = math::prod(elem.detJ[i],K1,S);
        math::add(1,1,K2,K);
    }
    return K;
}

darray Mesh::neumann(){

    darray B;
    double fLen = mesh.neuFace.size();
    double nLen = mesh.nXYZ.size();
    B.setlength(3*nLen);
    math::zero(B);

    for(int i=0; i<fLen; i++){

        int idx = mesh.neuFace[i];
        int nbr = fList[idx].nLen;
        matrix N = fList[idx].selfN(shape2D,quad2D);
        darray B1 = math::prod(1,N,mesh.neuValue[i]);

        // Inserts the subvectors into the global vector

        for(int j=0; j<nbr; j++){
            for(int k=0; k<3; k++){
                
                double val = B1(j+k*nbr);
                B(mesh.fNode[idx][j]+k*nLen) += val;
            }
        }
    }
    return B;
}

// Prepares the vector for Dirichlet BC

void Mesh::dirichlet(sparse&K,darray &B){

    double val;
    alglib::ae_int_t i=0,j=0;
    alglib::ae_int_t I=0,J=0;
    int nLen = mesh.nXYZ.size();
    vector<ivector> row(3*nLen);
    vector<ivector> col(3*nLen);

    // Stores the non-zero index per row and column

    while(alglib::sparseenumerate(K,I,J,i,j,val)){

        if(abs(val)<1e-9){alglib::sparseset(K,i,j,0);}
        else{row[i].push_back(j);col[j].push_back(i);}
    }

    // Edits B and cand cancels column/rows corresponding to the BC

    for(int n=0; n<3; n++){
        int dLen = mesh.dirNode[n].size();

        for(int m=0; m<dLen; m++){
            int idx = mesh.dirNode[n][m]+n*nLen;

            for(int k=0; k<col[idx].size(); k++){
                
                B[col[idx][k]] -= alglib::sparseget(K,col[idx][k],idx)*mesh.dirValue[n][m];
                if(col[idx][k]==idx){alglib::sparseset(K,col[idx][k],idx,1);}
                else{alglib::sparseset(K,col[idx][k],idx,0);}
            }
            for(int k=0; k<row[idx].size(); k++){

                if(idx==row[idx][k]){alglib::sparseset(K,idx,row[idx][k],1);}
                else{alglib::sparseset(K,idx,row[idx][k],0);}
            }
        }
    }

    // Places the BC value in the boundary vector

    for(int n=0; n<3; n++){
        for(int m=0; m<mesh.dirNode[n].size(); m++){
            B(mesh.dirNode[n][m]+n*nLen) = mesh.dirValue[n][m];
        }
    }
}