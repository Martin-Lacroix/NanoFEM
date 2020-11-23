#include "..\include\mesh.h"
#include <algorithm>
#include "solvers.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
using namespace std;

// generates a stiffness tensor

matrix stiffness(double E,double v){

    matrix D;
    D.setlength(6,6);
    math::zero(D);

    double mu = E/(2*(1+v));
    double lam = E*v/((1+v)*(1-2*v));
    D(0,1) = D(0,2) = D(1,2)= lam;
    D(1,0) = D(2,0) = D(2,1)= lam;

    for(int i=0; i<3; i++){

        D(i,i) = 2*mu+lam;
        D(i+3,i+3) = mu;
    }
    return D;
}

// Reads the nascam output file

void read(meshStruct &mesh, bcStruct &bc){

    // ------------------------------
    // Needs to read the input file

    mesh.order = 3;
    mesh.D = stiffness(1,0.3);

    // ------------------------------
    // Needs to read the input file

    darray x0BC("[0,0,0]");
    darray x1BC("[0.1,0,0]");
    darray y0BC("[0,0,0]");
    darray y1BC("[0,0,0]");
    darray z0BC("[0,0,0]");
    darray z1BC("[0,0,0]");

    // ------------------------------
    // Needs to read the input file

    vector<bool> x0Lock {1,0,0};
    vector<bool> x1Lock {0,0,0};
    vector<bool> y0Lock {0,1,0};
    vector<bool> y1Lock {0,0,0};
    vector<bool> z0Lock {0,0,1};
    vector<bool> z1Lock {0,0,0};

    // ------------------------------

    int nbr;
    string input;
    ifstream fileMesh;
    bc.dirichlet.resize(3);
    fileMesh.open("C:/Users/ORBBE/Desktop/TFE 2/Input/coating.xyz");
    fileMesh >> nbr;

    // Reads the size of the domain

    for(int i=0; i<2; i++){getline(fileMesh,input,';');}
    replace(input.begin(),input.end(),':',',');
    input = "["+input.substr(8)+"]";
    //darray domain(input.c_str());
    darray domain = "[1,1,0.5]";

    // Reads the origin of the domain

    getline(fileMesh,input,';');
    replace(input.begin(),input.end(),':',',');
    input = "["+input.substr(10)+"]";
    //darray zero(input.c_str());
    darray zero = "[0,0,0]";

    // Reads the size of an element

    getline(fileMesh,input,';');
    replace(input.begin(),input.end(),':',',');
    input = "["+input.substr(8)+"]";
    //darray elem(input.c_str());
    darray elem = "[0.5,0.5,0.5]";

    // Builds the mesh elements and nodes

    fileMesh.close();
    int nx = 0.1+domain(0)/elem(0);
    int ny = 0.1+domain(1)/elem(1);
    int nz = 0.1+domain(2)/elem(2);

    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            for(int k=0; k<=nz; k++){

                // Stores the nodes coordinates

                int idx = i*(ny+1)*(nz+1)+j*(nz+1)+k;
                double xyz[3] = {zero(0)+i*elem(0),zero(1)+j*elem(1),zero(2)+k*elem(2)};
                darray nXYZ; nXYZ.setcontent(3,xyz);
                mesh.nXYZ.push_back(nXYZ);

                // Creates the nodes BC

                if(i==0){
                    for(int n=0; n<3; n++){
                        if(x0Lock[n]){bc.dirichlet[0].push_back(idx);}
                    }
                }
                if(i==nx){
                    for(int n=0; n<3; n++){
                        if(x1Lock[n]){bc.dirichlet[0].push_back(idx);}
                    }
                }
                if(j==0){
                    for(int n=0; n<3; n++){
                        if(y0Lock[n]){bc.dirichlet[1].push_back(idx);}
                    }
                }
                if(j==ny){
                    for(int n=0; n<3; n++){
                        if(y1Lock[n]){bc.dirichlet[1].push_back(idx);}
                    }
                }
                if(k==0){
                    for(int n=0; n<3; n++){
                        if(z0Lock[n]){bc.dirichlet[2].push_back(idx);}
                    }
                }
                if(k==nz){
                    for(int n=0; n<3; n++){
                        if(z1Lock[n]){bc.dirichlet[2].push_back(idx);}
                    }
                }
            }
        }
    }

    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){

                // Associates the element nodes

                int idx = i*(nz+1)*(ny+1)+j*(nz+1)+k;
                int a[4] = {idx,(ny+1)*(nz+1)+idx,(ny+2)*(nz+1)+idx,nz+idx+1};
                int b[4] = {idx+1,(ny+1)*(nz+1)+idx+1,(ny+2)*(nz+1)+idx+1,nz+idx+2};
                alglib::ae_int_t node[8] = {a[0],a[1],a[2],a[3],b[0],b[1],b[2],b[3]};
                iarray eNode; eNode.setcontent(8,node);
                mesh.eNode.push_back(eNode);

                // Creates the element faces and BC

                if(i==0){

                    alglib::ae_int_t node[4] = {a[0],a[3],b[3],a[0]};
                    iarray fNode; fNode.setcontent(4,node);
                    mesh.fNode.push_back(fNode);
                    bc.neumann.push_back(x0BC);
                }
                if(i==nx-1){

                    alglib::ae_int_t node[4] = {a[1],b[1],b[2],a[2]};
                    iarray fNode; fNode.setcontent(4,node);
                    mesh.fNode.push_back(fNode);
                    bc.neumann.push_back(x1BC);
                }
                if(j==0){

                    alglib::ae_int_t node[4] = {a[0],b[0],b[1],a[1]};
                    iarray fNode; fNode.setcontent(4,node);
                    mesh.fNode.push_back(fNode);
                    bc.neumann.push_back(y0BC);
                }
                if(j==ny-1){

                    alglib::ae_int_t node[4] = {a[3],a[2],b[2],b[3]};
                    iarray fNode; fNode.setcontent(4,node);
                    mesh.fNode.push_back(fNode);
                    bc.neumann.push_back(y1BC);
                }
                if(k==0){
                    
                    alglib::ae_int_t node[4] = {a[0],a[1],a[2],a[3]};
                    iarray fNode; fNode.setcontent(4,node);
                    mesh.fNode.push_back(fNode);
                    bc.neumann.push_back(z0BC);
                }
                if(k==nz-1){
                    
                    alglib::ae_int_t node[4] = {b[3],b[2],b[1],b[0]};
                    iarray fNode; fNode.setcontent(4,node);
                    mesh.fNode.push_back(fNode);
                    bc.neumann.push_back(z1BC);
                }
            }
        }
    }
}

// Solves the sparse linear system Ku = B

darray solve(Mesh mesh,double p){

    sparse K;
    auto start = chrono::high_resolution_clock::now();
    auto stop = chrono::high_resolution_clock::now();
    auto time = chrono::duration_cast<std::chrono::microseconds>(stop-start);

    // Builds the system matrix

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
    cout << "\nBuilding K --- " << time.count()/1e6 << " sec";
    start = chrono::high_resolution_clock::now();

    // Applying boundary conditions

    mesh.dirichlet(K); math::clean(K);
    darray B = mesh.neumann();
    sparseconverttocrs(K);
    int nLen = B.length();
    mesh.dirichlet(B);

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nApplying BC --- " << time.count()/1e6 << " sec";
    start = chrono::high_resolution_clock::now();

    // Solves the linear system

    darray u;
    u.setlength(nLen);
    alglib::linlsqrreport rep;
    alglib::linlsqrstate state;
    alglib::linlsqrcreate(nLen,nLen,state);
    linlsqrsolvesparse(state,K,B);
    linlsqrresults(state,u,rep);

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nSolving System --- " << time.count()/1e6 << " sec\n\n";
    start = chrono::high_resolution_clock::now();
    return u;
}

// Main code

int main(){

    alglib::setglobalthreading(alglib::parallel);
    auto start = chrono::high_resolution_clock::now();
    auto stop = chrono::high_resolution_clock::now();
    auto time = chrono::duration_cast<std::chrono::microseconds>(stop-start);

    // Reads the input files

    bcStruct bcParam;
    meshStruct meshParam;
    read(meshParam,bcParam);

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nReading File --- " << time.count()/1e6 << " sec";
    start = chrono::high_resolution_clock::now();

    // Creates the mesh object

    Mesh mesh(meshParam,bcParam);

    stop = chrono::high_resolution_clock::now();
    time = chrono::duration_cast<std::chrono::microseconds>(stop-start);
    cout << "\nCreating Mesh --- " << time.count()/1e6 << " sec";

    // Solves the linear system

    darray u = solve(mesh,0);
    ofstream displacment("build/displacment.txt");
    for (int i=0; i<u.length(); i++){displacment << u[i] << "\n";}
}
