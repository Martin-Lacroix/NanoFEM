#include "..\include\parser.h"
#include "..\include\writer.h"
#include "solvers.h"
#include <direct.h>
#include <fstream>
#include <chrono>
using namespace std;

// ------------------------------------------------|
// Writes the simulation results in a text file    |
// ------------------------------------------------|
/*
void write(Mesh &mesh,darray &disp,vector<darray> &sigma){

    mkdir("output");
    ofstream uXYZ("output/disp.txt");
    ofstream elem("output/elem.txt");
    ofstream node("output/node.txt");
    ofstream stress("output/stress.txt");

    // Writes the displacement field in a text file

    for(int i=0; i<mesh.nLen; i++){
        for(int j=0; j<2; j++){uXYZ << disp[i+j*mesh.nLen] << ",";}
        uXYZ << disp[i+2*mesh.nLen] << "\n";
    }

    // Writes the node coordinates as (x,y,z)

    for(array3d nXYZ:mesh.data.nXYZ){
        for(int i=0; i<nXYZ.size()-1; i++){node << nXYZ[i] << ",";}
        node << nXYZ.back() << "\n";
    }

    // Writes the element nodes as (elem,node)

    for(ivector eNode:mesh.data.eNode){
        for(int i=0; i<eNode.size()-1; i++){elem << eNode[i] << ",";}
        elem << eNode.back() << "\n";
    }

    // Writes the averaged elemental stress in a text file

    for(darray s:sigma){
        for(int i=0; i<s.length()-1; i++){stress << s[i] << ",";}
        stress << s[s.length()-1] << "\n";
    }
}
*/
// -------------------------------------------------------|
// Solves the linear system in quasistatic equilibrium    |
// -------------------------------------------------------|

darray solve(Mesh &mesh){

    auto stop = chrono::high_resolution_clock::now();
    auto start = chrono::high_resolution_clock::now();
    auto clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nBuilds the matrix --- ";

    // Builds the system matrix and vector

    int nLen = 3*mesh.nLen;
    int size = 9*mesh.eLen*pow(mesh.data.order+1,6)/4;

    sparse K;
    darray B;
    B.setlength(nLen);
    alglib::sparsecreate(nLen,nLen,size,K);
    math::zero(B);

    mesh.totalKB(K,B);
    mesh.neumann(B);

    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << clock.count()/1e6 << " sec\n";
    cout << "Boundary conditions --- ";

    // Applies boundary conditions to K and B

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


    mesh.delta(K,B);
    mesh.coupling(K,B);
    mesh.dirichlet(K,B);
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

    dataStruct data;
    string path = "input.txt";
    read(path,data);

    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << clock.count()/1e6 << " sec";

    // Creates the mesh class and solves with conjugate gradient

    Mesh mesh(move(data));
    darray disp = solve(mesh);
    start = chrono::high_resolution_clock::now();
    cout << "Stress extraction --- ";

    // Computes Von Mises stresses and updates the nodes

    vector<darray> sigma = mesh.stress(disp);
    mesh.update(disp);

    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    start = chrono::high_resolution_clock::now();
    cout << clock.count()/1e6 << " sec\n";

    // Writes the results in a text file

    start = chrono::high_resolution_clock::now();
    cout << "Writes the results --- ";
    write(mesh,disp,sigma);
    writeJmol(mesh,disp,sigma);

    // Prints the computation time of the operation

    stop = chrono::high_resolution_clock::now();
    clock = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << clock.count()/1e6 << " sec\n\n";

    /*
    cout << "\n";
    for(int i=0; i<mesh.nLen; i++){
        cout << "Node " << i << " -- ux = " << disp[i] << "\n";
    }
    cout << "\n";
    for(int i=0; i<mesh.nLen; i++){
        cout << "Node " << i << " -- uy = " << disp[i+mesh.nLen] << "\n";
    }
    cout << "\n";
    for(int i=0; i<mesh.nLen; i++){
        cout << "Node " << i << " -- uz = " << disp[i+2*mesh.nLen] << "\n";
    }
    cout << "\n";
    */
}