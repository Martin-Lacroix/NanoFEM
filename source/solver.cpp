#include "..\include\read.h"
#include "solvers.h"
#include <direct.h>
#include <fstream>
#include <chrono>
using namespace std;

// Solves the sparse symmetric linear system

darray solve(Mesh mesh,double p){

    sparse K;
    auto stop = chrono::high_resolution_clock::now();
    auto start = chrono::high_resolution_clock::now();
    auto time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "Builds the matrix --- ";

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

    // Gets computation time

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << time.count()/1e6 << " sec\n";
    cout << "Boundary conditions --- ";

    // Applies boundary conditions

    darray B = mesh.neumann();
    int nLen = B.length();
    mesh.periodic(K,B);

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
    mesh.dirichlet(K,B);
    sparseconverttocrs(K);

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << time.count()/1e6 << " sec\n";
    cout << "Solves the system --- ";
    
    // Solves the symmetric linear system

    darray u;
    u.setlength(nLen);
    alglib::lincgreport rep;
    alglib::lincgstate state;
    alglib::lincgcreate(nLen,state);
    alglib::lincgsolvesparse(state,K,1,B);
    alglib::lincgresults(state,u,rep);
    mesh.complete(u);

    // Gets computation time

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << time.count()/1e6 << " sec\n";
    return u;
}

// Main code of the finite element solver

int main(){

    alglib::setglobalthreading(alglib::parallel);
    auto stop = chrono::high_resolution_clock::now();
    auto start = chrono::high_resolution_clock::now();
    auto time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nReads the files --- ";

    // Reads the input files

    string inputPath = "input.txt";
    string meshPath = "input/test.xyz";
    meshStruct mesh = read(inputPath,meshPath);

    // Gets computation time

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << time.count()/1e6 << " sec\n";
    cout << "Creates the mesh --- ";

    // Creates the mesh object

    Mesh Mesh(mesh);

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << time.count()/1e6 << " sec\n";

    // Solves the linear system with conjugate gradient

    darray u = solve(Mesh,0);

    cout << "\n\n";
    for(int i=0; i<u.length()/3; i++){
        cout << "Node " << i << " -- ux = " << u[i] << "\n";
    }
    cout << "\n";
    for(int i=0; i<u.length()/3; i++){
        cout << "Node " << i << " -- uy = " << u[i+u.length()/3] << "\n";
    }
    cout << "\n";
    for(int i=0; i<u.length()/3; i++){
        cout << "Node " << i << " -- uz = " << u[i+2*u.length()/3] << "\n";
    }


    // Writes the results in a text file

    mkdir("output");
    start = chrono::high_resolution_clock::now();
    ofstream coordinates("output/coordinates.txt");
    ofstream displacement("output/displacement.txt");
    cout << "Writes the results --- ";

    for(int i=0; i<u.length(); i++){displacement << u[i] << "\n";}
    for(int i=0; i<Mesh.mesh.nXYZ.size(); i++){

        coordinates << Mesh.mesh.nXYZ[i][0] << ",";
        coordinates << Mesh.mesh.nXYZ[i][1] << ",";
        coordinates << Mesh.mesh.nXYZ[i][2] << "\n";
    }

    // Gets computation time

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << time.count()/1e6 << " sec\n\n";
}