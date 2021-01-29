#include "..\include\read.h"
#include <unordered_map>
#include <iterator>
#include <fstream>
#include <sstream>
using namespace std;

// Converts a string to a vector of doubles

dvector tovec(string input){

    replace(input.begin(),input.end(),':',' ');
    replace(input.begin(),input.end(),',',' ');

    double nbr;
    dvector vector;
    stringstream iss(input);
    while (iss>>nbr){vector.push_back(nbr);}
    return vector;
}

// Reads the parameter input file

void readInput(readStruct &read,meshStruct &mesh,string path){

    double E,v;
    string input;
    ifstream file;
    file.open(path);

    // Reads the general parameters

    getline(file,input,';');
    mesh.order = stoi(input);
    getline(file,input,' ');
    read.threshold = stod(input);
    getline(file,input,'\n');

    // Reads the bulk parameters

    getline(file,input,';'); E = stod(input);
    getline(file,input,' '); v = stod(input);
    read.Db = math::stiffness(E,v);
    getline(file,input,'\n');

    // Reads the hole parameters

    getline(file,input,';'); E = stod(input);
    getline(file,input,' '); v = stod(input);
    read.Dh = math::stiffness(E,v);
    getline(file,input,'\n');

    // Reads the boundary conditions

    for(int i=0; i<3; i++){
        getline(file,input,'\n');

        if(input=="axial strain"){

            getline(file,input,'\n');
            spair pair = make_pair("axial strain",input);
            read.boundary.push_back(pair);
        }
        else{

            spair pair = make_pair(input,"none");
            read.boundary.push_back(pair);
        }
    }
}

// Reads the mesh file

void readMesh(readStruct &read,meshStruct &mesh,string path){

    int nbr;
    string input;
    ifstream file;
    vector<ivector> neighbour;
    mesh.perNode.resize(3);
    mesh.dirNode.resize(3);
    mesh.dirVal.resize(3);
    file.open(path);
    file >> nbr;

    // Reads the size of the domain

    for(int i=0; i<2; i++){getline(file,input,';');}
    dvector dom = tovec(input.substr(8));

    // Reads the origin of the domain

    getline(file,input,';');
    dvector zero = tovec(input.substr(10));

    // Reads the size of an element

    getline(file,input,';');
    dvector elem = tovec(input.substr(8));

    // Number of elements in each dimension

    int nx = 0.1+dom[0]/elem[0];
    int ny = 0.1+dom[1]/elem[1];
    int nz = 0.1+dom[2]/elem[2];

    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            for(int k=0; k<=nz; k++){

                // Creates the node coordinates

                int ID,max;
                int idx = i*(ny+1)*(nz+1)+j*(nz+1)+k;
                mesh.nXYZ.push_back({zero[0]+i*elem[0],zero[1]+j*elem[1],zero[2]+k*elem[2]});

                for(int n=0; n<3; n++){
                    switch(n){

                        case 0: ID=i; max=nx; break;
                        case 1: ID=j; max=ny; break;
                        case 2: ID=k; max=nz; break;
                    }

                    // Stores the nodes boundary conditions

                    if(ID==0){
                        if(read.boundary[n].first=="clamped"){

                            mesh.dirNode[n].push_back(idx);
                            mesh.dirVal[n].push_back(0);
                        }
                        else if(read.boundary[n].first=="periodic"){
                            mesh.perNode[n].first.push_back(idx);
                        }
                        else if(read.boundary[n].first=="axial strain"){
                            
                            double val = stod(read.boundary[n].second)*dom[n]/2;
                            mesh.dirNode[n].push_back(idx);
                            mesh.dirVal[n].push_back(-val);
                        }
                    }
                    else if(ID==max){
                        if(read.boundary[n].first=="axial strain"){

                            double val = stod(read.boundary[n].second)*dom[n]/2;
                            mesh.dirNode[n].push_back(idx);
                            mesh.dirVal[n].push_back(val);
                        }
                        else if(read.boundary[n].first=="periodic"){
                            mesh.perNode[n].second.push_back(idx);
                        }
                    }
                }
            }
        }
    }

    // Builds the elements and the faces of the mesh

    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){

                int eLen = mesh.eNode.size();
                int idx = i*(ny+1)*(nz+1)+j*(nz+1)+k;
                neighbour.push_back({-1,-1,-1,-1,-1,-1});

                // Stores the nodes of each element

                ivector node = {idx,(ny+1)*(nz+1)+idx,(ny+2)*(nz+1)+idx,nz+idx+1};
                node.insert(node.end(),{node[0]+1,node[1]+1,node[2]+1,node[3]+1});
                mesh.eNode.push_back(node);

                // Stores the list of neightbour elements

                if(i!=0){neighbour.back()[0] = eLen-ny*nz;}
                if(j!=0){neighbour.back()[2] = eLen-nz;}
                if(k!=0){neighbour.back()[4] = eLen-1;}

                if(i!=nx-1){neighbour.back()[1] = eLen+ny*nz;}
                if(j!=ny-1){neighbour.back()[3] = eLen+nz;}
                if(k!=nz-1){neighbour.back()[5] = eLen+1;}
            }
        }
    }

    // Reads the filling fraction of the elements

    getline(file,input,'\n');
    int eLen = mesh.eNode.size();
    dvector fraction(eLen,0);
    mesh.eSurf.resize(eLen);
    mesh.D.resize(eLen);

    while(getline(file,input,' ')){

        // Coodrinates of the species

        getline(file,input,';');
        dvector coord = tovec(input);
        getline(file,input,';');
        getline(file,input,'\n');

        // Filling fraction of the elements

        int dx = nx*coord[0]/dom[0];
        int dy = ny*coord[1]/dom[1];
        int dz = nz*coord[2]/dom[2];
        int index = dx*ny*nz+dy*nz+dz;
        fraction[index] += stod(input);
    }

    // Creates the element stiffness tensors

    for(int i=0; i<eLen; i++){

        if(fraction[i]<read.threshold){
            mesh.D[i] = read.Dh;
        }

        // Stores the free surfaces of bulk elements

        else{
            mesh.D[i] = read.Db;
            for(int j=0; j<6; j++){

                int idx = neighbour[i][j];
                if(fraction[idx]<read.threshold){mesh.eSurf[i].push_back(j);}
            }
        }
    }
    file.close();
}

// Reads all the Nascam input files

meshStruct read(string inputPath,string meshPath){

    meshStruct mesh;
    readStruct read;
    readInput(read,mesh,inputPath);
    readMesh(read,mesh,meshPath);
    return mesh;
}