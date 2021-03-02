#include "..\include\read.h"
#include "solvers.h"
#include <direct.h>
#include <fstream>
#include <chrono>
using namespace std;

// ------------------------------------------------|
// Writes the simulation results in a text file    |
// ------------------------------------------------|

void write(Mesh &mesh,timeStruct &time,vector<darray> &uList,vector<darray> &sigma){

    mkdir("output");
    ofstream elements("output/elements.txt");
    ofstream coordinates("output/coordinates.txt");
    ofstream displacement("output/displacement.txt");

    // Writes the displacement field in a text file

    if(time.u0.length()>0){
        for(int i=0; i<uList.size(); i+=time.nSave){

            for(int j=0; j<uList[i].length()-1; j++){displacement << uList[i][j] << ",";}
            displacement << uList[i][uList[i].length()-1] << "\n";
        }
    }
    else{
        for(darray u:uList){
            for(int i=0; i<u.length()-1; i++){displacement << u[i] << ",";}
            displacement << u[u.length()-1] << "\n";
        }
    }

    // Writes the node coordinates as (x,y,z)

    for(array3d nXYZ:mesh.mesh.nXYZ){
        for(int i=0; i<nXYZ.size()-1; i++){coordinates << nXYZ[i] << ",";}
        coordinates << nXYZ.back() << "\n";
    }

    // Writes the element nodes as (elem,node)

    for(ivector eNode:mesh.mesh.eNode){
        for(int i=0; i<eNode.size()-1; i++){elements << eNode[i] << ",";}
        elements << eNode.back() << "\n";
    }

    // Writes the averaged elemental stress in a text file

    if(sigma.size()>0){
        ofstream stress("output/stress.txt");

        for(darray s:sigma){
            for(int i=0; i<s.length()-1; i++){stress << s[i] << ",";}
            stress << s[s.length()-1] << "\n";
        }
    }
}

// --------------------------------------------------|
// Solves the linear system with wave propagation    |
// --------------------------------------------------|

vector<darray> solveWave(Mesh &mesh,timeStruct &wave){

    auto stop = chrono::high_resolution_clock::now();
    auto start = chrono::high_resolution_clock::now();
    auto clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nBuilds the matrix --- ";

    // Builds the full K and M matrices of the system

    sparse M1;
    int nLen = 3*mesh.nLen;
    sparse M = mesh.totalM();
    sparse K = mesh.totalK();
    darray B = mesh.neumann();
    double dt = pow(wave.dt,2);

    alglib::sparseconverttocrs(K);
    alglib::sparsecopytocrs(M,M1);

    // Initializes the solution vector

    vector<darray> u;
    darray u1 = wave.u0;
    darray x,y,u2 = u1;
    u.push_back(wave.u0);

    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << clock.count()/1e6 << " sec\n";
    cout << "Time iterations --- ";

    // Applies boundary conditions to M1 and B

    for(int i=0; i<wave.nSteps; i++){

        alglib::sparsecopytohashbuf(M1,M);

        // Builds the right-hand-side of the equation

        math::add(2,-1,u2,u1);
        alglib::sparsesmv(M1,1,u1,x);
        alglib::sparsesmv(K,1,u2,y);
        math::add(1,-1,B,y);
        math::add(1,dt,x,y);
        u1 = u2;

        // Sets the boundary conditions

        mesh.delta(M,y);
        mesh.coupling(M,y);
        mesh.dirichlet(M,y);
        alglib::sparseconverttocrs(M);

        // Solves the symmetric linear system with Alglib

        alglib::lincgreport rep;
        alglib::lincgstate state;
        alglib::lincgcreate(nLen,state);
        alglib::lincgsolvesparse(state,M,1,y);
        alglib::lincgresults(state,u2,rep);
        mesh.complete(u2);
        u.push_back(u2);
    }

    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << clock.count()/1e6 << " sec\n";
    return u;
}

// -------------------------------------------------------|
// Solves the linear system in quasistatic equilibrium    |
// -------------------------------------------------------|

darray solveStatic(Mesh &mesh){

    auto stop = chrono::high_resolution_clock::now();
    auto start = chrono::high_resolution_clock::now();
    auto clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nBuilds the matrix --- ";

    // Builds the full K matrix of the system

    int nLen = 3*mesh.nLen;
    sparse K = mesh.totalK();
    darray B = mesh.neumann();

/*
    sparse K2 = mesh.totalK2();
    darray B2;
    B2.setlength(3*nLen);

    for(int i=0; i<nLen; i++){

        B2(i) = B(i);
        B2(i+nLen) = 0;
    }

    K = K2;
    B = B2;
    nLen = 3*nLen;
    */
    
    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << clock.count()/1e6 << " sec\n";
    cout << "Boundary conditions --- ";

    // Applies boundary conditions to K and B

    mesh.delta(K,B);
    mesh.coupling(K,B);
    mesh.dirichlet(K,B);

/*
    ofstream Kfile("output/K.txt");
    ofstream Bfile("output/B.txt");

    for(int i=0; i<nLen; i++){
        for(int j=0; j<nLen; j++){
            
            double val = alglib::sparseget(K,i,j);
            Kfile << val << " ";
        }
        Kfile << "\n";
        Bfile << B[i] << "\n";
    }
*/



    alglib::sparseconverttocrs(K);

    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << clock.count()/1e6 << " sec\n";
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
    clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << clock.count()/1e6 << " sec\n";
    return u;
}

// ------------------------------------------|
// Main code of the finite element solver    |
// ------------------------------------------|

int main(){

    alglib::setglobalthreading(alglib::parallel);
    auto stop = chrono::high_resolution_clock::now();
    auto start = chrono::high_resolution_clock::now();
    auto clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nReads the files --- ";

    // Reads the input files from Nascam

    meshStruct mesh;
    timeStruct time;
    string path[2] = {"input.txt","input/test.xyz"};
    read(path,mesh,time);

    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << clock.count()/1e6 << " sec";

    // Creates the mesh class and solves with conjugate gradient

    Mesh Mesh(mesh);
    vector<darray> uList;
    vector<darray> sigma;

    // Computes the quasitatic solution

    if(time.u0.length()==0){

        darray u = solveStatic(Mesh);
        uList.push_back(u);

        // Computes Von Mises stresses and updates the nodes

        start = chrono::high_resolution_clock::now();
        cout << "Stress extraction --- ";
        sigma = Mesh.stress(u);
        Mesh.update(u);

        // Prints the computation time of the operation

        stop = chrono::high_resolution_clock::now();
        clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
        start = chrono::high_resolution_clock::now();
        cout << clock.count()/1e6 << " sec\n";
    }

    // Computes the wave propagation

    else{uList = solveWave(Mesh,time);}

    // Writes the results in a text file

    start = chrono::high_resolution_clock::now();
    cout << "Writes the results --- ";
    write(Mesh,time,uList,sigma);

    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << clock.count()/1e6 << " sec\n\n";

    cout << "\n";
    for(int i=0; i<Mesh.nLen; i++){
        cout << "Node " << i << " -- ux = " << uList.back()[i] << "\n";
    }
    cout << "\n";
    for(int i=0; i<Mesh.nLen; i++){
        cout << "Node " << i << " -- uy = " << uList.back()[i+Mesh.nLen] << "\n";
    }
    cout << "\n";
    for(int i=0; i<Mesh.nLen; i++){
        cout << "Node " << i << " -- uz = " << uList.back()[i+2*Mesh.nLen] << "\n";
    }
    cout << "\n";
    

}