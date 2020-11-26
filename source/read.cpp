#include "..\include\read.h"
#include <fstream>
using namespace std;

// Generates the stiffness tensor

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

// Reads the parameter input file

readStruct readParam(paramStruct &param){

    string input;
    ifstream file;
    readStruct bc;
    file.open("input.txt");

    // Reads the general parameters

    getline(file,input,';');
    param.order = stoi(input);

    getline(file,input,';');
    double E = stod(input);

    getline(file,input,' ');
    double v = stod(input);

    param.D = stiffness(E,v);

    // Reads the Dirichlet boundary conditions

    getline(file,input,'\n');
    for(int i=0; i<5; i++){

        getline(file,input,';');
        bc.fixed.push_back(input.c_str());
    }

    getline(file,input,' ');
    bc.fixed.push_back(input.c_str());

    // Reads the Neumann boundary conditions

    getline(file,input,'\n');
    for(int i=0; i<5; i++){

        getline(file,input,';');
        bc.force.push_back(input.c_str());
    }

    getline(file,input,' ');
    bc.force.push_back(input.c_str());
    file.close();
    return bc;
}

// Reads all the Nascam input files

void readAll(meshStruct &mesh, paramStruct &param){

    readStruct bc = readParam(param);

    int nbr;
    string input;
    ifstream file;
    file.open("input/coating.xyz");
    param.dirichlet.resize(3);
    file >> nbr;

    // Reads the size of the domain

    for(int i=0; i<2; i++){getline(file,input,';');}
    replace(input.begin(),input.end(),':',',');
    input = "["+input.substr(8)+"]";
    darray domain(input.c_str());

    // Reads the origin of the domain

    getline(file,input,';');
    replace(input.begin(),input.end(),':',',');
    input = "["+input.substr(10)+"]";
    darray zero(input.c_str());

    // Reads the size of an element

    getline(file,input,';');
    replace(input.begin(),input.end(),':',',');
    input = "["+input.substr(8)+"]";
    darray elem(input.c_str());

    // Builds the mesh elements and nodes

    file.close();
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

                // Creates the nodes Dirichlet BC

                if(i==0){
                    for(int n=0; n<3; n++){
                        if(bc.fixed[0][n]){param.dirichlet[0].push_back(idx);}
                    }
                }
                if(i==nx){
                    for(int n=0; n<3; n++){
                        if(bc.fixed[1][n]){param.dirichlet[0].push_back(idx);}
                    }
                }
                if(j==0){
                    for(int n=0; n<3; n++){
                        if(bc.fixed[2][n]){param.dirichlet[1].push_back(idx);}
                    }
                }
                if(j==ny){
                    for(int n=0; n<3; n++){
                        if(bc.fixed[3][n]){param.dirichlet[1].push_back(idx);}
                    }
                }
                if(k==0){
                    for(int n=0; n<3; n++){
                        if(bc.fixed[4][n]){param.dirichlet[2].push_back(idx);}
                    }
                }
                if(k==nz){
                    for(int n=0; n<3; n++){
                        if(bc.fixed[5][n]){param.dirichlet[2].push_back(idx);}
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

                // Creates the element faces and Neumann BC

                if(i==0){

                    alglib::ae_int_t node[4] = {a[0],a[3],b[3],b[0]};
                    iarray fNode; fNode.setcontent(4,node);
                    param.neumann.push_back(bc.force[0]);
                    mesh.fNode.push_back(fNode);
                }
                if(i==nx-1){

                    alglib::ae_int_t node[4] = {a[1],b[1],b[2],a[2]};
                    iarray fNode; fNode.setcontent(4,node);
                    param.neumann.push_back(bc.force[1]);
                    mesh.fNode.push_back(fNode);
                }
                if(j==0){

                    alglib::ae_int_t node[4] = {a[0],b[0],b[1],a[1]};
                    iarray fNode; fNode.setcontent(4,node);
                    param.neumann.push_back(bc.force[2]);
                    mesh.fNode.push_back(fNode);
                }
                if(j==ny-1){

                    alglib::ae_int_t node[4] = {a[3],a[2],b[2],b[3]};
                    iarray fNode; fNode.setcontent(4,node);
                    param.neumann.push_back(bc.force[3]);
                    mesh.fNode.push_back(fNode);
                }
                if(k==0){
                    
                    alglib::ae_int_t node[4] = {a[0],a[1],a[2],a[3]};
                    iarray fNode; fNode.setcontent(4,node);
                    param.neumann.push_back(bc.force[4]);
                    mesh.fNode.push_back(fNode);
                }
                if(k==nz-1){
                    
                    alglib::ae_int_t node[4] = {b[3],b[2],b[1],b[0]};
                    iarray fNode; fNode.setcontent(4,node);
                    param.neumann.push_back(bc.force[5]);
                    mesh.fNode.push_back(fNode);
                }
            }
        }
    }
}
