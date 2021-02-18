#include "..\include\mesh.h"
using namespace std;

// -------------------------------------------------------------|
// Builds the elements list, shape functions and quadratures    |
// -------------------------------------------------------------|

Mesh::Mesh(meshStruct &input){
    
    mesh = input;
    nLen = mesh.nXYZ.size();
    eLen = mesh.eNode.size();
    fLen = mesh.neuFace.size();

    // Stores the quadrature rules and shape functions

    quad3D = math::legendre(3,mesh.ordQuad);
    quad2D = math::legendre(2,mesh.ordQuad);
    shape3D = shape(3,mesh.ordElem);
    shape2D = shape(2,mesh.ordElem);
}

// -------------------------------------------------------------|
// Returns the structure of 2D and 3D linear shape functions    |
// -------------------------------------------------------------|

shapeStruct Mesh::shape(int dim,int order){

    dvector val;
    shapeStruct shape;
    dvector node(order+1);
    double step = 2.0/order;

    // Length of the node list and Gauss points

    if(dim==3){
        shape.gLen = quad3D.gRST.size();
        shape.nLen = 0.1+pow(order+1,3);
    }
    if(dim==2){
        shape.gLen = quad2D.gRST.size();
        shape.nLen = 0.1+pow(order+1,2);
    }

    // Memory allocation and 1D node list

    shape.N.setlength(shape.nLen,shape.gLen);
    shape.drN.setlength(shape.nLen,shape.gLen);
    shape.dsN.setlength(shape.nLen,shape.gLen);
    shape.dtN.setlength(shape.nLen,shape.gLen);

    for(int i=0; i<=order; i++){node[i] = i*step-1;}

    // Sets the coordinates of the Gauss points

    for(int i=0; i<shape.gLen; i++){

        if(dim==2){val = {quad2D.gRST[i][0],quad2D.gRST[i][1]};}
        if(dim==3){val = {quad3D.gRST[i][0],quad3D.gRST[i][1],quad3D.gRST[i][2]};}

        // Computes the Lagrange shape functions at Gauss points

        dvector drN = math::lagrange(0,node,val);
        dvector dsN = math::lagrange(1,node,val);
        dvector dtN = math::lagrange(2,node,val);
        dvector N = math::lagrange(-1,node,val);

        for(int j=0; j<shape.nLen; j++){

            // Stores the shape functions evaluated at the Gauss points
            
            shape.N(j,i) = N[j];
            shape.drN(j,i) = drN[j];
            shape.dsN(j,i) = dsN[j];
            shape.dtN(j,i) = dtN[j];
        }
    }
    return shape;
}

// --------------------------------------------------------|
// Computes the total stiffness matrix K for global FEM    |
// --------------------------------------------------------|

sparse Mesh::totalK(){

    sparse K;
    alglib::sparsecreate(3*nLen,3*nLen,K);

    for(int i=0; i<eLen; i++){

        // Coordinates of the nodes of the element
        
        int sLen = mesh.eNode[i].size();
        vector<dvector> eXYZ(sLen,dvector(3));
        for(int j=0; j<sLen; j++){eXYZ[j] = mesh.nXYZ[mesh.eNode[i][j]];}

        // Computes the elemental K matrices

        Elem elem(eXYZ,shape3D);
        ivector eNode = mesh.eNode[i];
        matrix D = math::stiffness(mesh.Ev[i][0],mesh.Ev[i][1]);
        matrix K1 = elem.selfK(quad3D,D);

        // Inserts the elemental matrix into the global K matrix

        for(int j=0; j<sLen; j++){
            for(int k=0; k<sLen; k++){
                for(int m=0; m<3; m++){
                    for(int n=0; n<3; n++){

                        int row = eNode[j]+m*nLen;
                        int col = eNode[k]+n*nLen;

                        // Adds the element only if located in the upper triangular

                        if(row<=col){

                            double val = K1[j+m*sLen][k+n*sLen];
                            alglib::sparseadd(K,eNode[j]+m*nLen,eNode[k]+n*nLen,val);
                        }
                    }
                }
            }
        }
    }
    return K;
}

// ---------------------------------------------------|
// Computes the total mass matrix M for global FEM    |
// ---------------------------------------------------|

sparse Mesh::totalM(){

    sparse M;
    alglib::sparsecreate(3*nLen,3*nLen,M);

    for(int i=0; i<eLen; i++){

        // Coordinates of the nodes of the element
        
        int sLen = mesh.eNode[i].size();
        vector<dvector> eXYZ(sLen,dvector(3));
        for(int j=0; j<sLen; j++){eXYZ[j] = mesh.nXYZ[mesh.eNode[i][j]];}

        // Computes the elemental M matrices

        double rho = 1;

        Elem elem(eXYZ,shape3D);
        ivector eNode = mesh.eNode[i];
        matrix M1 = elem.selfM(shape3D,quad3D,rho);

        // Inserts the elemental matrix into the global M matrix

        for(int j=0; j<sLen; j++){
            for(int k=0; k<sLen; k++){
                for(int m=0; m<3; m++){
                    for(int n=0; n<3; n++){

                        int row = eNode[j]+m*nLen;
                        int col = eNode[k]+n*nLen;

                        // Adds the element only if located in the upper triangular

                        if(row<=col){

                            double val = M1[j+m*sLen][k+n*sLen];
                            alglib::sparseadd(M,eNode[j]+m*nLen,eNode[k]+n*nLen,val);
                        }
                    }
                }
            }
        }
    }
    return M;
}

// -------------------------------------------------------------|
// Computes the vector of applied stresses B from Neumann BC    |
// -------------------------------------------------------------|

darray Mesh::neumann(){

    darray B;
    B.setlength(3*nLen);
    math::zero(B);

    // Computes the applied surface tractions on the 2D elements

    for(int i=0; i<fLen; i++){

        // Coordinates of the nodes of the faces

        int sLen = mesh.neuFace[i].size();
        vector<dvector> eXYZ(sLen,dvector(3));
        for(int j=0; j<sLen; j++){eXYZ[j] = mesh.nXYZ[mesh.neuFace[i][j]];}

        // Computes the elemental B vectors

        Face face(eXYZ,shape2D);
        matrix M = face.selfM(shape2D,quad2D);
        darray B1 = math::prod(1,M,mesh.neuVal[i]);

        // Inserts the elemental vectors into the global B vector

        for(int j=0; j<sLen; j++){
            for(int k=0; k<3; k++){
                B(mesh.neuFace[i][j]+k*nLen) += B1(j+k*sLen);
            }
        }
    }
    return B;
}

// ----------------------------------------------------------------|
// Edits the matrix K and vector B for fixed nodal displacement    |
// ----------------------------------------------------------------|

void Mesh::dirichlet(sparse &K,darray &B){

    vector<ivector> row = math::sparsemap(K);

    // Gets the number of nodes in the dirichlet list

    for(int n=0; n<3; n++){
        int dLen = mesh.dirNode[n].size();

        // Gets the index of the node in the solution vector

        for(int i=0; i<dLen; i++){
            int idx = mesh.dirNode[n][i]+n*nLen;

            // Substract the coefficients from K to B and cancels the row-colums

            for(int j:row[idx]){

                B(j) -= math::symget(K,j,idx)*mesh.dirVal[n][i];
                math::symset(K,j,idx,0);
            }
            row[idx].clear();
        }
    }

    // Places the dirichlet value the B vector and 1 in the diagonal of K

    for(int n=0; n<3; n++){
        for(int i=0; i<mesh.dirNode[n].size(); i++){
            
            int idx = mesh.dirNode[n][i]+n*nLen;
            alglib::sparseset(K,idx,idx,1);
            B(idx) = mesh.dirVal[n][i];
        }
    }
}

// ------------------------------------------------------------------|
// Edits the matrix K and vector B for coupled nodal displacement    |
// ------------------------------------------------------------------|

void Mesh::coupling(sparse &K,darray &B){

    vector<ivector> row = math::sparsemap(K);

    for(int n=0; n<3; n++){

        // Gets the number of nodes in the periodic list

        for(int i=0; i<mesh.coupNode[n].size(); i++){
            int dLen = mesh.coupNode[n][i].size();

            // Gets the index of the curent and final coupled nodes

            for(int j=0; j<dLen-1; j++){

                int idx1 = mesh.coupNode[n][i][j]+n*nLen;
                int idx2 = mesh.coupNode[n][i].back()+n*nLen;
                double fix = math::symget(K,idx2,idx2)+math::symget(K,idx1,idx1)+2*math::symget(K,idx1,idx2);

                // Adds the current row to the final row of K if non-zero

                for(int k=0; k<row[idx1].size(); k++){

                    int n = row[idx1][k];
                    double val = math::symget(K,idx1,n);

                    // Updates the mapStruct only if an element is added to K

                    if(val!=0){
                        if(math::symget(K,idx2,n)==0){

                            row[n].push_back(idx2);
                            row[idx2].push_back(n);
                        }
                        math::symadd(K,idx2,n,val);
                    }
                }

                // Zeroes the column and row in K corresoponding to the coupled node

                for(int k:row[idx1]){math::symset(K,idx1,k,0);}
                row[idx1].clear();
    
                // Corrects the diagonal elements in K and edits B

                alglib::sparseset(K,idx2,idx2,fix);
                alglib::sparseset(K,idx1,idx1,1);
                B(idx2) += B(idx1);
                B(idx1) = 0;
            }
        }
    }
}

// -------------------------------------------------------------------|
// Edits the matrix K and vector B for a change of variable u = Δu    |
// -------------------------------------------------------------------|

void Mesh::delta(sparse &K,darray &B){

    vector<ivector> row = math::sparsemap(K);

    // Gets the number of nodes in the periodic list

    for(int n=0; n<3; n++){
        for(int i=0; i<mesh.deltaNode[n].size(); i++){

            // Gets the index of the curent and final coupled nodes

            int idx1 = mesh.deltaNode[n][i].first+n*nLen;
            int idx2 = mesh.deltaNode[n][i].second+n*nLen;
            double fix = math::symget(K,idx2,idx2)+math::symget(K,idx1,idx1)+2*math::symget(K,idx1,idx2);

            // Adds the current row to the final row of K if non-zero

            for(int k=0; k<row[idx1].size(); k++){

                int n = row[idx1][k];
                double val = math::symget(K,idx1,n);

                // Updates the mapStruct only if an element is added to K

                if(val!=0){
                    if(math::symget(K,idx2,n)==0){

                        row[n].push_back(idx2);
                        row[idx2].push_back(n);
                    }
                    math::symadd(K,idx2,n,val);
                }
            }

            // Corrects the diagonal elements in K and edits B

            alglib::sparseset(K,idx2,idx2,fix);
            B(idx2) += B(idx1);
        }
    }
}

// -----------------------------------------------------|
// Completes the solution for coupled or delta nodes    |
// -----------------------------------------------------|

void Mesh::complete(darray &u){

    for(int n=0; n<3; n++){

        // Copy the value of the last node into the coupled ones

        for(int i=0; i<mesh.coupNode[n].size(); i++){
            for(int j=0; j<mesh.coupNode[n][i].size(); j++){

                int idx1 = mesh.coupNode[n][i][j]+n*nLen;
                int idx2 = mesh.coupNode[n][i].back()+n*nLen;
                u[idx1] = u[idx2];
            }
        }

        // Comes back from the difference Δu to the actual displacement

        for(int i=0; i<mesh.deltaNode[n].size(); i++){

            int idx1 = mesh.deltaNode[n][i].first+n*nLen;
            int idx2 = mesh.deltaNode[n][i].second+n*nLen;
            u[idx1] += u[idx2];
        }
    }
    
}

// ---------------------------------------------------------|
// Updates the coordinates of the node with the solution    |
// ---------------------------------------------------------|

void Mesh::update(darray &u){

    // Adds the displacement to the corresponding node coordinate

    for(int i=0; i<nLen; i++){
        for(int j=0; j<3; j++){
            mesh.nXYZ[i][j] += u[i+j*nLen];
        }
    }
}

// ----------------------------------------------------------|
// Computes the averaged Von Mises stress in each element    |
// ----------------------------------------------------------|

vector<darray> Mesh::stress(darray &u){

    vector<darray> sigma(eLen);
    for(int i=0; i<eLen; i++){

        // Coordinates of the nodes of the element

        int sLen = mesh.eNode[i].size();
        vector<dvector> eXYZ(sLen,dvector(3));
        for(int j=0; j<sLen; j++){eXYZ[j] = mesh.nXYZ[mesh.eNode[i][j]];}

        // Stores the displacement of the curent element nodes

        darray u1;
        u1.setlength(3*sLen);

        for(int j=0; j<sLen; j++){
            for(int k=0; k<3; k++){
                u1[j+k*sLen] = u(mesh.eNode[i][j]+k*nLen);
            }
        }

        // Computes the averaged Von Mises stress

        Elem elem(eXYZ,shape3D);
        matrix D = math::stiffness(mesh.Ev[i][0],mesh.Ev[i][1]);
        sigma[i] = elem.stress(quad3D,D,u1);
    }
    return sigma;
}