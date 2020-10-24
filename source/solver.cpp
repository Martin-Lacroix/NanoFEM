#include "..\include\mesh.h"
#include "solvers.h"
using namespace std;

meshStruct testMesh(){

    matrix D;
    meshStruct param;
    vector<darray> nXYZ;
    vector<iarray> eNode;
    D.setlength(6,6);

    // Node coordinates

    nXYZ.push_back("[0,0,0]");
    nXYZ.push_back("[0.5,0,0]");
    nXYZ.push_back("[1,0,0]");
    nXYZ.push_back("[1,0.5,0]");
    nXYZ.push_back("[0.5,0.5,0]");
    nXYZ.push_back("[0,0.5,0]");
    nXYZ.push_back("[0,1,0]");
    nXYZ.push_back("[0.5,1,0]");
    nXYZ.push_back("[1,1,0]");
    nXYZ.push_back("[0,0,0.5]");
    nXYZ.push_back("[0.5,0,0.5]");
    nXYZ.push_back("[1,0,0.5]");
    nXYZ.push_back("[1,0.5,0.5]");
    nXYZ.push_back("[0.5,0.5,0.5]");
    nXYZ.push_back("[0,0.5,0.5]");
    nXYZ.push_back("[0,1,0.5]");
    nXYZ.push_back("[0.5,1,0.5]");
    nXYZ.push_back("[1,1,0.5]");

    // Element nodes

    eNode.push_back("[0,1,4,5,9,10,13,14]");
    eNode.push_back("[4,3,12,13,1,2,11,10]");
    eNode.push_back("[6,7,16,15,5,4,13,14]");
    eNode.push_back("[8,7,4,3,17,16,13,12]");

    // Stiffness tensor

    double E = 1;
    double v = 0.3;
    double mu = E/(2*(1+v));
    double lam = E*v/((1+v)*(1-2*v));

    for(int i=0; i<6; i++){
        for(int j=0; j<6; j++){

            if((i<3) and (j<3) and (i!=j)){D(i,j)= lam;}
            else if((i<3) and (j<3) and (i==j)){D(i,j) = 2*mu+lam;}
            else if((i>=3) and (j>=3) and (i==j)){D(i,j)= mu;}
            else{D(i,j) = 0;}
        }
    }

    // Place the result in a struct

    param.D = D;
    param.nXYZ = nXYZ;
    param.eNode = eNode;
    return param;
}

// Creates the test boundary conditions (to be deleted)

bcStruct testBC(){

    bcStruct param;
    vector<darray> fValBC;
    vector<darray> nValBC;
    vector<iarray> nIdxBC;
    vector<iarray> fNodeBC;

    // Value and face-node indices of Neumann BC

    fNodeBC.push_back("[2,11,12,3]");
    fNodeBC.push_back("[8,3,12,17]");

    fValBC.push_back("[0.1,0,0]");
    fValBC.push_back("[0.1,0,0]");

    // Value and node indices of Dirichlet BC

    nIdxBC.push_back("[0,5,6,9,14,15]");
    nIdxBC.push_back("[0,1,2,9,10,11]");
    nIdxBC.push_back("[0,1,2,3,4,5,6,7,8]");

    nValBC.push_back("[0,0,0,0,0,0]");
    nValBC.push_back("[0,0,0,0,0,0]");
    nValBC.push_back("[0,0,0,0,0,0,0,0,0]");

    // Place the result in a struct

    param.fNode = fNodeBC;
    param.fVal = fValBC;
    param.nVal = nValBC;
    param.nIdx = nIdxBC;
    return param;
}

// Solves the sparse linear system Ku = B

darray solve(Mesh mesh){
   
    sparse K = mesh.localK();
    darray B = mesh.neumann();
    int nLen = B.length();

    // Applies boundary conditions

    mesh.dirichlet(B);
    mesh.dirichlet(K);

    darray u; u.setlength(nLen);
    alglib::sparsesolverreport rep;
    alglib::sparsesolve(K,nLen,B,u,rep);

    return u;
}

// Main code

int main(){

    bcStruct bcParam = testBC();
    meshStruct meshParam = testMesh();
    Mesh mesh(meshParam,bcParam);
    darray u = solve(mesh);



    for(int i=0; i<u.length(); i++){cout << u[i] << "\n";}
    cout << "\nDone\n" << endl;

}
