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

// Keeps only unique faces in the mesh

void cleanFace(meshStruct &mesh){

    ivector fElem;
    vector<string> liste;
    vector<ivector> fNode;
    int fLen = mesh.fNode.size();
    unordered_map<string,int> map;

   for(int i=0; i<fLen; i++){

        ostringstream key;
        copy(mesh.fNode[i].begin(),mesh.fNode[i].end(),ostream_iterator<int>(key,":"));
        liste.push_back(key.str());
        map[liste[i]] += 1;
   }

   for(int i=0; i<fLen; i++){
        if(map[liste[i]]==1){

            fElem.push_back(mesh.fElem[i]);
            fNode.push_back(mesh.fNode[i]);
        }
   }
    mesh.fElem = fElem;
    mesh.fNode = fNode;
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

    // Reads the boundary Conditions

    for(int i=0; i<3; i++){

        getline(file,input,'\n');
        if(input!="imposed strain"){
            read.dirichlet.push_back(input);
        }
        else{
            getline(file,input,'\n');
            read.dirichlet.push_back(input);
        }
    }
}

// Reads the mesh file

void readMesh(readStruct &read,meshStruct &mesh,string path){

    int nbr;
    string input;
    ifstream file;
    mesh.dirNode.resize(3);
    mesh.dirValue.resize(3);
    svector dirichlet = read.dirichlet;
    file.open(path);
    file >> nbr;

    // Reads the size of the domain

    for(int i=0; i<2; i++){getline(file,input,';');}
    dvector domain = tovec(input.substr(8));

    // Reads the origin of the domain

    getline(file,input,';');
    dvector zero = tovec(input.substr(10));

    // Reads the size of an element

    getline(file,input,';');
    dvector elem = tovec(input.substr(8));

    // Number of elements in each dimensions

    int nx = 0.1+domain[0]/elem[0];
    int ny = 0.1+domain[1]/elem[1];
    int nz = 0.1+domain[2]/elem[2];

    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            for(int k=0; k<=nz; k++){

                // Creates the node coordinates

                int now = 0;
                int idx = i*(ny+1)*(nz+1)+j*(nz+1)+k;
                mesh.nXYZ.push_back({zero[0]+i*elem[0],zero[1]+j*elem[1],zero[2]+k*elem[2]});

                // Creates the nodes Dirichlet BC

                for(int n=0; n<3; n++){

                    switch(n){
                    case 0: now = i; break;
                    case 1: now = j; break;
                    case 2: now = k; break;
                    }

                    if(now==0){
                        if(dirichlet[n]=="clamped bottom"){

                            mesh.dirValue[n].push_back(0);
                            mesh.dirNode[n].push_back(idx);
                        }
                        else if(dirichlet[n]=="free periodic"){
                        }
                        else{
                            mesh.dirValue[n].push_back(-stod(dirichlet[n])*domain[n]/2);
                            mesh.dirNode[n].push_back(idx);
                        }
                    }
                    else if(now==nx){
                        if(dirichlet[n]=="free periodic"){
                        }
                        else if(dirichlet[n]!="clamped bottom"){

                            mesh.dirValue[n].push_back(stod(dirichlet[n])*domain[n]/2);
                            mesh.dirNode[n].push_back(idx);
                        }
                    }
                }
            }
        }
    }

    cout << "\n";
    for(int i=0; i<mesh.dirValue[0].size(); i++){
        cout << "Node " << mesh.dirNode[0][i] << " - DirValue ux = " << mesh.dirValue[0][i] << "\n";
    }
    cout << "\n";
    for(int i=0; i<mesh.dirValue[1].size(); i++){
        cout << "Node " << mesh.dirNode[1][i] << " - DirValue uy = " << mesh.dirValue[2][i] << "\n";
    }
    cout << "\n";
    for(int i=0; i<mesh.dirValue[2].size(); i++){
        cout << "Node " << mesh.dirNode[2][i] << " - DirValue uy = " << mesh.dirValue[2][i] << "\n";
    }

    // Builds the elements and the faces of the mesh

    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            for(int k=0; k<nz; k++){

                vector<ivector> face(6,ivector(4));
                int idx = i*(ny+1)*(nz+1)+j*(nz+1)+k;
                face[0] = {idx,(ny+1)*(nz+1)+idx,(ny+2)*(nz+1)+idx,nz+idx+1};
                face[1] = {face[0][0]+1,face[0][1]+1,face[0][2]+1,face[0][3]+1};

                // Makes the left, right, top and bottom faces

                face[2] = {face[0][1],face[0][2],face[1][2],face[1][1]};
                face[3] = {face[0][0],face[0][3],face[1][3],face[1][0]};
                face[4] = {face[0][0],face[0][1],face[1][1],face[1][0]};
                face[5] = {face[0][3],face[0][2],face[1][2],face[1][3]};

                // Stores the nodes of each element

                ivector elem = face[0];
                elem.insert(elem.end(),face[1].begin(),face[1].end());
                mesh.eNode.push_back(elem);

                // Stores the faces of each element

                for(int n=0; n<6; n++){

                    mesh.fNode.push_back(face[n]);
                    mesh.fElem.push_back(mesh.eNode.size());
                }
            }
        }
    }

    // Reads the filling fraction of the elements

    getline(file,input,'\n');
    int eLen = mesh.eNode.size();
    dvector fraction(eLen,0);
    mesh.D.resize(eLen);

    while(getline(file,input,' ')){

        // Coodrinates of the species

        getline(file,input,';');
        dvector coord = tovec(input);
        getline(file,input,';');
        getline(file,input,'\n');

        // Filling fraction of the elements

        int idx = nx*coord[0]/domain[0];
        int idy = ny*coord[1]/domain[1];
        int idz = nz*coord[2]/domain[2];
        int index = idx*ny*nz+idy*nz+idz;
        fraction[index] += stod(input);
    }

    // Creates the element stiffness tensors

    for(int i=0; i<eLen; i++){

        if(fraction[i]>=read.threshold){mesh.D[i] = read.Db;}
        else{mesh.D[i] = read.Dh;}
    }
    file.close();
}

// Reads all the Nascam input files

meshStruct read(string inputPath,string meshPath){

    meshStruct mesh;
    readStruct read;
    readInput(read,mesh,inputPath);
    readMesh(read,mesh,meshPath);
    cleanFace(mesh);
    return mesh;
}