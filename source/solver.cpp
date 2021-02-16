#include "..\include\read.h"
#include "solvers.h"
#include <direct.h>
#include <fstream>
#include <chrono>
using namespace std;

// ---------------------------------------------|
// Solves the sparse symmetric linear system    |
// ---------------------------------------------|

darray solve(Mesh mesh){

    sparse K;
    auto stop = chrono::high_resolution_clock::now();
    auto start = chrono::high_resolution_clock::now();
    auto time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nBuilds the matrix --- ";

    // Builds the full K matrix of the system

    K = mesh.localK();

    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << time.count()/1e6 << " sec\n";
    cout << "Boundary conditions --- ";

    // Applies boundary conditions to K and B

    darray B = mesh.neumann();
    int nLen = B.length();
/*
    ofstream Kfile("output/K.txt");
    ofstream Bfile("output/B.txt");

    for(int i=0; i<B.length(); i++){
        for(int j=0; j<B.length(); j++){
            
            double val = alglib::sparseget(K,i,j);
            Kfile << val << " ";
        }
        Kfile << "\n";
        Bfile << B[i] << "\n";
    }
*/
    mesh.delta(K,B);
    mesh.coupling(K,B);
    mesh.dirichlet(K,B);
    sparseconverttocrs(K);

    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << time.count()/1e6 << " sec\n";
    cout << "Solves the system --- ";
    
    // Solves the symmetric linear system with Alglib

    darray u;
    u.setlength(nLen);
    alglib::lincgreport rep;
    alglib::lincgstate state;
    alglib::lincgcreate(nLen,state);
    alglib::lincgsolvesparse(state,K,1,B);
    alglib::lincgresults(state,u,rep);
    mesh.complete(u);

    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << time.count()/1e6 << " sec\n";
    return u;
}

// ------------------------------------------|
// Main code of the finite element solver    |
// ------------------------------------------|

int main(){

    alglib::setglobalthreading(alglib::parallel);
    auto stop = chrono::high_resolution_clock::now();
    auto start = chrono::high_resolution_clock::now();
    auto time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nReads the files --- ";

    // Reads the input files from Nascam

    string inputPath = "input.txt";
    string meshPath = "input/single Cu.xyz";
    meshStruct mesh = read(inputPath,meshPath);

    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << time.count()/1e6 << " sec";

    // Creates the mesh class and solves with conjugate gradient

    Mesh Mesh(mesh);
    darray u = solve(Mesh);
    start = chrono::high_resolution_clock::now();
    cout << "Stress extraction --- ";
    
    // Computes Von Mises stresses and update the nodes

    int nLen = mesh.nXYZ.size();
    int eLen = mesh.eNode.size();
    vector<darray> sigma = Mesh.stress(u);
    Mesh.update(u);

    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << time.count()/1e6 << " sec\n";
    cout << "Writes the results --- ";

    // Writes the results in a text file

    mkdir("output");
    ofstream stress("output/stress.txt");
    ofstream elements("output/elements.txt");
    ofstream coordinates("output/coordinates.txt");
    ofstream displacement("output/displacement.txt");

    // Writes the displacement field in a text file

    for(int i=0; i<nLen; i++){
        for(int j=0; j<2; j++){displacement << u[i+j*nLen] << ",";}
        displacement << u[i+2*nLen] << "\n";
    }

    // Writes the node coordinates (x,y,z) in a text file

    for(dvector nXYZ:Mesh.mesh.nXYZ){
        for(int j=0; j<2; j++){coordinates << nXYZ[j] << ",";}
        coordinates << nXYZ.back() << "\n";
    }

    // Writes the element nodes in a text file

    for(ivector eNode:Mesh.mesh.eNode){
        for(int j=0; j<eNode.size()-1; j++){elements << eNode[j] << ",";}
        elements << eNode.back() << "\n";
    }

    // Writes the elemental stress (σxx,σyy,σxy,σzx,σyz) in a text file

    for(int i=0; i<eLen; i++){
        for(int j=0; j<5; j++){stress << sigma[i][j] << ",";}
        stress << sigma[i][5] << "\n";
    }

    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << time.count()/1e6 << " sec\n\n";
/*
    cout << "\n";
    for(int i=0; i<nLen; i++){
        cout << "Node " << i << " -- ux = " << u[i] << "\n";
    }
    cout << "\n";
    for(int i=0; i<nLen; i++){
        cout << "Node " << i << " -- uy = " << u[i+nLen] << "\n";
    }
    cout << "\n";
    for(int i=0; i<nLen; i++){
        cout << "Node " << i << " -- uz = " << u[i+2*nLen] << "\n";
    }
    cout << "\n";
*/
}