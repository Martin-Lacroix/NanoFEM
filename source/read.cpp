#include "..\include\read.h"
#include <unordered_map>
#include <iterator>
#include <fstream>
#include <sstream>
using namespace std;

// --------------------------------------------|
// Converts a string to a vector of doubles    |
// --------------------------------------------|

dvector tovec(string input){

    replace(input.begin(),input.end(),':',' ');
    replace(input.begin(),input.end(),',',' ');

    // Stores the values into the vector

    double nbr;
    dvector vector;
    stringstream iss(input);
    while (iss>>nbr){vector.push_back(nbr);}
    return vector;
}

// -------------------------------------------------------------|
// Stores the boundary conditions of the nodes in meshStruct    |
// -------------------------------------------------------------|

void setBC(readStruct &read,meshStruct &mesh,ivector loop){

    ivector dLen = read.dLen;
    ivector add = {dLen[0]*(dLen[1]+1)*(dLen[2]+1),dLen[1]*(dLen[2]+1),dLen[2]};
    int idx = loop[0]*(dLen[1]+1)*(dLen[2]+1)+loop[1]*(dLen[2]+1)+loop[2];

    // Check node position for the 3 indices

    for(int i=0; i<3; i++){
        string bc = read.boundary[i].first;

        // If the curent node is at the bottom

        if(loop[i]==0){

            // If the bottom face is fixed

            if(bc=="clamped"){

                mesh.dirNode[i].push_back(idx);
                mesh.dirVal[i].push_back(0);
            }

            // If the bottom face is fixed and top face is flat

            else if(bc=="axial stress" || bc=="periodic"){

                mesh.dirNode[i].push_back(idx);
                mesh.dirVal[i].push_back(0);

                if(read.row[i][idx+add[i]]<0){

                    mesh.perNode[i][0].push_back(idx+add[i]);
                    read.row[i][idx+add[i]] = 0;
                }
            }

            // If the bottom face is fixed and imposed strain on top face

            else if(bc=="axial strain"){
                
                mesh.dirNode[i].push_back(idx);
                mesh.dirVal[i].push_back(0);
            }

            // Coupled displacement in transversal dimensions if periodic

            for(int j=0; j<2; j++){
                int k = (i+j+1)%3;

                if(bc=="periodic"){

                    // Gets location of the node pair in the list of coupled nodes

                    int loc1 = read.row[k][idx];
                    int loc2 = read.row[k][idx+add[i]];

                    if(loc1>=0 && loc2<0){

                        // If the first node of the pair is already in the list

                        mesh.perNode[k][loc1].push_back(idx+add[i]);
                        read.row[k][idx+add[i]] = loc1;
                    }
                    else if(loc1<0 && loc2>=0){

                        // If the second node of the pair is already in the list

                        mesh.perNode[k][loc2].push_back(idx);
                        read.row[k][idx] = loc2;
                    }
                    else if(loc1<0 && loc2<0){

                        // If no of the nodes of the pair are already in the list

                        mesh.perNode[k].push_back({idx,idx+add[i]});
                        read.row[k][idx+add[i]] = mesh.perNode[k].size()-1;
                        read.row[k][idx] = mesh.perNode[k].size()-1;
                    }
                }
            }
        }
        
        // If the curent node is at the top and imposed strain

        else if(loop[i]==read.dLen[i]){
            if(bc=="axial strain"){

                double val = read.boundary[i].second*read.dSize[i];
                mesh.dirNode[i].push_back(idx);
                mesh.dirVal[i].push_back(val);
            }
        }
    }
}

// -----------------------------------------------|
// Reads the parameter from the input.txt file    |
// -----------------------------------------------|

void readInput(readStruct &read,meshStruct &mesh,string path){

    string input;
    ifstream file;
    file.open(path);

    // Reads the order of the quadrature rule

    getline(file,input,' ');
    mesh.order = stoi(input);
    getline(file,input,'\n');

    // Reads the parameters of empty elements

    getline(file,input,';'); read.E.push_back(stod(input));
    getline(file,input,' '); read.v.push_back(stod(input));
    getline(file,input,'\n');

    // Reads the parameters of bulk elements

    getline(file,input,';'); read.E.push_back(stod(input));
    getline(file,input,' '); read.v.push_back(stod(input));
    getline(file,input,'\n');

    // Reads the boundary conditions

    for(int i=0; i<3; i++){
        getline(file,input,'\n');

        // Imposed strain along the i-th dimension

        if(input=="axial strain"){

            getline(file,input,'\n');
            sdpair pair = make_pair("axial strain",stod(input));
            read.boundary.push_back(pair);
        }

        // Imposed stress at top face along the i-th dimension

        else if(input=="axial stress"){

            getline(file,input,'\n');
            sdpair pair = make_pair("axial stress",stod(input));
            read.boundary.push_back(pair);
        }

        // If clamped bottom face or no boundary conditions
        
        else{
            sdpair pair = make_pair(input,0);
            read.boundary.push_back(pair);
        }
    }
}

// -------------------------------------------------|
// Reads the parameter from the coating.xyz file    |
// -------------------------------------------------|

void readMesh(readStruct &read,meshStruct &mesh,string path){

    string input;
    ifstream file;
    vector<ivector> neighbour;
    mesh.perNode.resize(3);
    mesh.dirNode.resize(3);
    mesh.dirVal.resize(3);
    file.open(path);

    // Reads the size of the domain

    getline(file,input,'\n');
    for(int i=0; i<2; i++){getline(file,input,';');}
    dvector dom = tovec(input.substr(8));

    // Reads the origin of the domain

    getline(file,input,';');
    dvector zero = tovec(input.substr(10));
    read.dSize = {dom[0]-zero[0],dom[1]-zero[1],dom[2]-zero[2]};

    // Reads the size of a 8-node finite element

    getline(file,input,';');
    dvector eSize = tovec(input.substr(8));

    // Computes number of elements in each dimension and prepares the BC

    for(int i=0; i<3; i++){

        string bc = read.boundary[i].first;
        if(bc=="axial stress" || bc=="periodic"){mesh.perNode[i].resize(1);}
        read.dLen.push_back(0.1+read.dSize[i]/eSize[i]);
    }

    // Number of elements and initialization of the BC tracker

    ivector dLen = read.dLen;
    int nLen = (dLen[0]+1)*(dLen[1]+1)*(dLen[2]+1);
    read.row.resize(3,ivector(nLen,-1));

    // Stores the node coordinates and the BC in the meshStruct

    for(int i=0; i<=dLen[0]; i++){
        for(int j=0; j<=dLen[1]; j++){
            for(int k=0; k<=dLen[2]; k++){

                ivector loop = {i,j,k};
                mesh.nXYZ.push_back({zero[0]+i*eSize[0],zero[1]+j*eSize[1],zero[2]+k*eSize[2]});
                setBC(read,mesh,loop);
            }
        }
    }

    // Builds the elements and the faces of the mesh

    for(int i=0; i<dLen[0]; i++){
        for(int j=0; j<dLen[1]; j++){
            for(int k=0; k<dLen[2]; k++){
                
                ivector node;
                int eLen = mesh.eNode.size();
                neighbour.push_back({-1,-1,-1,-1,-1,-1});
                int idx = i*(dLen[1]+1)*(dLen[2]+1)+j*(dLen[2]+1)+k;

                // Stores the nodes of each element in meshStruct

                node = {idx,(dLen[1]+1)*(dLen[2]+1)+idx,(dLen[1]+2)*(dLen[2]+1)+idx,dLen[2]+idx+1};
                node.insert(node.end(),{node[0]+1,node[1]+1,node[2]+1,node[3]+1});
                mesh.eNode.push_back(node);

                // Stores the list of neightbour elements on 3 bottom faces

                if(i!=0){neighbour.back()[0] = eLen-dLen[1]*dLen[2];}
                if(j!=0){neighbour.back()[2] = eLen-dLen[2];}
                if(k!=0){neighbour.back()[4] = eLen-1;}

                // Stores the list of neightbour elements on 3 top faces

                if(i!=dLen[0]-1){neighbour.back()[1] = eLen+dLen[1]*dLen[2];}
                if(j!=dLen[1]-1){neighbour.back()[3] = eLen+dLen[2];}
                if(k!=dLen[2]-1){neighbour.back()[5] = eLen+1;}
            }
        }
    }

    // Initializes the filling fraction of the elements

    getline(file,input,'\n');
    int eLen = mesh.eNode.size();
    dvector tol = {-1e-5,1e-5};
    dvector fraction(eLen,0);
    mesh.eSurf.resize(eLen);
    mesh.D.resize(eLen);

    while(getline(file,input,' ')){

        // Reafds the coodrinates of the chemical species

        unordered_set<int> set;
        getline(file,input,';');
        dvector coord = tovec(input);
        getline(file,input,';');
        getline(file,input,'\n');

        // Slightly moves the species to check if located at boundary between elements

        for(int i=0; i<2; i++){

            int dx = dLen[0]*coord[0]/read.dSize[0]+tol[i];
            if(dx>=dLen[0]){continue;}

            // Moves the species along y axis

            for(int j=0; j<2; j++){

                int dy = dLen[1]*coord[1]/read.dSize[1]+tol[j];
                if(dy>=dLen[1]){continue;}

                // Moves the species along z axis

                for(int k=0; k<2; k++){
                            
                    int dz = dLen[2]*coord[2]/read.dSize[2]+tol[k];
                    if(dz<dLen[2]){set.insert(dx*dLen[1]*dLen[2]+dy*dLen[2]+dz);}
                }
            }
        }

        // Splits the filling fraction between the elements

        for (int i:set){
            fraction[i] += stod(input)/set.size();
        }
    }

    for(int i=0; i<eLen; i++){

        // Computes the mixed stiffness tensor of the elements

        if(fraction[i]>1){fraction[i] = 1;}
        double E = fraction[i]*read.E[1]+(1-fraction[i])*read.E[0];
        double v = fraction[i]*read.v[1]+(1-fraction[i])*read.v[0];
        mesh.D[i] = math::stiffness(E,v);

        for(int j=0; j<3; j++){
            if(read.boundary[j].first=="axial stress"){

                // Neumann BC applied to the top face perpendicular to the j-th dimension
                    
                ivector node = mesh.eNode[i];
                vector<ivector> face(3,ivector(4));
                double val = read.boundary[j].second;

                // Computes the faces of the 8-node element

                face[0] = {node[1],node[2],node[6],node[5]};
                face[1] = {node[2],node[6],node[7],node[3]};
                face[2] = {node[4],node[5],node[6],node[7]};

                if(neighbour[i][2*j+1]==-1){

                    // Stores the applied stress on the face without neighbour
                                    
                    mesh.neuNode.push_back(face[j]);
                    mesh.neuVal.push_back("[0,0,0]");
                    mesh.neuVal.back()[j] = val;
                }
            }
        }
    }
    file.close();
}

// -------------------------------------------------------|
// Reads the Nascam input files return the meshStruct     |
// -------------------------------------------------------|

meshStruct read(string inputPath,string meshPath){

    meshStruct mesh;
    readStruct read;

    // Calls the two mesh reading and building functions

    readInput(read,mesh,inputPath);
    readMesh(read,mesh,meshPath);
    return mesh;
}