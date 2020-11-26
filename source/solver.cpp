#include "..\include\read.h"
#include "solvers.h"
#include <direct.h>
#include <fstream>
#include <chrono>
using namespace std;

// Solves the sparse symmetric linear system Ku = B

darray solve(Mesh mesh,double p){

    sparse K;
    auto stop = chrono::high_resolution_clock::now();
    auto start = chrono::high_resolution_clock::now();
    auto time = chrono::duration_cast<std::chrono::microseconds>(stop-start);

    // Builds the full system matrix

    if(p==0){
        K = mesh.localK();
    }
    else if(p==1){
        K = mesh.nonLocalK();
    }
    else{
        K = mesh.nonLocalK();
        sparse KL = mesh.localK();
        math::add(1-p,p,KL,K);
    }

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nBuilds the matrix --- " << time.count()/1e6 << " sec";
    start = chrono::high_resolution_clock::now();

    // Applies boundary conditions

    mesh.dirichlet(K);
    darray B = mesh.neumann();
    sparseconverttocrs(K);
    int nLen = B.length();
    mesh.dirichlet(B);

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nBoundary conditions --- " << time.count()/1e6 << " sec";
    start = chrono::high_resolution_clock::now();
    
    // Solves the symmetric linear system

    darray u;
    u.setlength(nLen);
    alglib::lincgreport rep;
    alglib::lincgstate state;
    alglib::lincgcreate(nLen,state);
    alglib::lincgsolvesparse(state,K,1,B);
    alglib::lincgresults(state,u,rep);

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nSolves the system --- " << time.count()/1e6 << " sec";
    start = chrono::high_resolution_clock::now();
    return u;
}

// Main code of the finite element solver

int main(){

    alglib::setglobalthreading(alglib::parallel);
    auto stop = chrono::high_resolution_clock::now();
    auto start = chrono::high_resolution_clock::now();
    auto time = chrono::duration_cast<std::chrono::microseconds>(stop-start);

    // Reads the input files

    meshStruct mesh;
    paramStruct param;
    readAll(mesh,param);

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nReads the files --- " << time.count()/1e6 << " sec";
    start = chrono::high_resolution_clock::now();

    // Creates the mesh object

    Mesh Mesh(mesh,param);

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nCreates the mesh --- " << time.count()/1e6 << " sec";

    // Solves the symmetric linear system

    darray u = solve(Mesh,0);

    // Writes the results in a file

    mkdir("output");
    start = chrono::high_resolution_clock::now();
    ofstream coordinates("output/coordinates.txt");
    ofstream displacement("output/displacement.txt");

    for(int i=0; i<u.length(); i++){displacement << u[i] << "\n";}
    for(int i=0; i<Mesh.mesh.nXYZ.size(); i++){

        coordinates << Mesh.mesh.nXYZ[i][0] << ",";
        coordinates << Mesh.mesh.nXYZ[i][1] << ",";
        coordinates << Mesh.mesh.nXYZ[i][2] << "\n";
    }

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nWrites the results --- " << time.count()/1e6 << " sec\n\n";
}