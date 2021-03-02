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
    replace(input.begin(),input.end(),';',' ');
    replace(input.begin(),input.end(),',',' ');

    // Stores the values into the vector

    double nbr;
    dvector vector;
    stringstream iss(input);
    while (iss>>nbr){vector.push_back(nbr);}
    return vector;
}

// -----------------------------------------------|
// Reads the parameter from the input.txt file    |
// -----------------------------------------------|

void readInput(readStruct &read,meshStruct &mesh,string path){

    string input;
    ifstream file;
    file.open(path);

    // Reads the order of the quadrature rule

    getline(file,input,';');
    mesh.order = stoi(input);
    getline(file,input,';');
    read.Lc = stod(input);
    getline(file,input,'!');
    read.cropZ = stod(input);
    getline(file,input,'\n');

    // Reads the parameters of empty elements

    getline(file,input,'!');
    read.emptyEv = tovec(input);
    getline(file,input,'\n');

    // Reads the parameters of substrate element
        
    getline(file,input,'!');
    read.Ev.push_back(tovec(input));
    read.EvS.push_back({0,0,0});
    getline(file,input,'\n');

    // Reads the type of applied stress

    getline(file,input,'\n');
    read.type = input;

    // Reads the axis of applied stress

    getline(file,input,';');
    read.load = input;

    for(int i=0; i<2; i++){

        if(input=="All"){read.axis[i] = {0,1,2};}
        else if(input[i+1]=='X'){read.axis[i] = {0};}
        else if(input[i+1]=='Y'){read.axis[i] = {1};}
        else if(input[i+1]=='Z'){read.axis[i] = {2};}
    }

    // Reads the value of the action

    getline(file,input,'!');
    read.Fval = stod(input);
    getline(file,input,'\n');

    // Stores the Young modulus and Poisson ratio of the layers

    while(getline(file,input,'!')){

        dvector param = tovec(input);
        read.Ev.push_back({param[0],param[1]});
        read.EvS.push_back({param[2],param[3]});
        getline(file,input,'\n');
    }
}

// -----------------------------------------------|
// Reads the parameter from the input.xyz file    |
// -----------------------------------------------|

void readMeshSize(readStruct &read,meshStruct &mesh,string path){

    string input;
    ifstream file;
    file.open(path);
    int order = mesh.order;

    // Reads the size of the cubic domain

    getline(file,input,'\n');
    for(int i=0; i<2; i++){getline(file,input,';');}
    read.dSize = tovec(input.substr(8));

    // Reads the origin of the domain

    getline(file,input,';');
    read.zero = tovec(input.substr(10));
    dvector zero = read.zero;

    // Reads the size of a 8-node finite element

    getline(file,input,';');
    read.eSize = tovec(input.substr(8));
    dvector eSize = read.eSize;

    // Truncates the height of the domain to the closest element

    if(read.cropZ>eSize[0]/2){
        
        read.cropZ -= eSize[2]/2;
        read.cropZ += eSize[2]-fmod(read.cropZ,eSize[2]);
        read.dSize[2] -= read.cropZ;
    }

    // Number of elements and initialization of the BC tracker

    for(int i=0; i<3; i++){read.dLen.push_back(0.1+read.dSize[i]/eSize[i]);}
    ivector dLen = read.dLen;

    // Stores the node coordinates and the BC in the mesh

    for(int i=0; i<=dLen[0]*order; i++){
        for(int j=0; j<=dLen[1]*order; j++){
            for(int k=0; k<=dLen[2]*order; k++){

                mesh.nXYZ.push_back({0,0,0});
                mesh.nXYZ.back()[0] = zero[0]+i*eSize[0]/order;
                mesh.nXYZ.back()[1] = zero[1]+j*eSize[1]/order;
                mesh.nXYZ.back()[2] = zero[2]+k*eSize[2]/order;
            }
        }
    }

    // Builds the elements and the faces of the mesh

    for(int i=0; i<dLen[0]; i++){
        for(int j=0; j<dLen[1]; j++){
            for(int k=0; k<dLen[2]; k++){
                
                mesh.eNode.push_back({});
                int eLen = mesh.eNode.size();
                read.neighbour.push_back(ivector(6,-1));

                // Computes the first index and geometrical spacing

                int dy = order*dLen[2]+1;
                int dx = (order*dLen[1]+1)*(order*dLen[2]+1);
                int idx = i*(dLen[1]*order+1)*(dLen[2]*order+1)*order+j*(dLen[2]*order+1)*order+k*order;

                // Stores the nodes of each element in the mesh

                for(int n=0; n<order+1; n++){
                    for(int m=0; m<order+1; m++){
                        for(int p=0; p<order+1; p++){

                            mesh.eNode.back().push_back(idx+p+m*dy+n*dx);
                        }
                    }
                }

                // Stores the list of neightbour elements on 3 bottom faces

                if(i!=0){read.neighbour.back()[0] = eLen-dLen[1]*dLen[2];}
                if(j!=0){read.neighbour.back()[2] = eLen-dLen[2];}
                if(k!=0){read.neighbour.back()[4] = eLen-1;}

                // Stores the list of neightbour elements on 3 top faces

                if(i!=dLen[0]-1){read.neighbour.back()[1] = eLen+dLen[1]*dLen[2];}
                if(j!=dLen[1]-1){read.neighbour.back()[3] = eLen+dLen[2];}
                if(k!=dLen[2]-1){read.neighbour.back()[5] = eLen+1;}
            }
        }
    }
    file.close();
}

// -------------------------------------------------------|
// Localizes the chemical species in the mesh elements    |
// -------------------------------------------------------|

unordered_set<int> locSpecies(readStruct &read,dvector coord){

    unordered_set<int> eList;
    ivector dLen = read.dLen;
    dvector zero = read.zero;
    dvector tol = {-1e-5,1e-5};

    for(int i=0; i<2; i++){

        // Slightly moves the species along y axis

        int dx = dLen[0]*(coord[0]-zero[0])/read.dSize[0]+tol[i];
        if(dx>=dLen[0]){continue;}

        // Slightly moves the species along y axis

        for(int j=0; j<2; j++){

            int dy = dLen[1]*(coord[1]-zero[1])/read.dSize[1]+tol[j];
            if(dy>=dLen[1]){continue;}

            // Slightly moves the species along z axis

            for(int k=0; k<2; k++){
                        
                int dz = dLen[2]*(coord[2]-zero[2])/read.dSize[2]+tol[k];
                if(dz<dLen[2]){eList.insert(dx*dLen[1]*dLen[2]+dy*dLen[2]+dz);}
            }
        }
    }
    return eList;
}

// ---------------------------------------------------------|
// Computes the stiffness tensor from the input.xyz file    |
// ---------------------------------------------------------|

void readSpecies(readStruct &read,meshStruct &mesh,string path){

    string input;
    ifstream file;
    file.open(path);
    getline(file,input,'\n');
    getline(file,input,'\n');

    // Initializes the filling fraction of the elements

    int eLen = mesh.eNode.size();
    vector<dvector> frac(eLen,dvector(read.Ev.size(),0));
    read.empty.resize(eLen);
    mesh.EvR.resize(eLen);
    mesh.EvS.resize(eLen);

    // Reads the coodrinates of the chemical species

    while(getline(file,input,' ')){

        getline(file,input,';');
        dvector coord = tovec(input);
        getline(file,input,';');

        // Gets the layer of the chemical species

        int layer = stoi(input);
        getline(file,input,'\n');

        // Splits the filling fraction between the elements

        unordered_set<int> eList = locSpecies(read,coord);
        for (int i:eList){frac[i][layer] += stod(input)/eList.size();}
    }

    // Stores the mixed mechanical parameters of the elements

    for(int i=0; i<eLen; i++){
        double sum = accumulate(frac[i].begin(),frac[i].end(),0.0);

        if(sum<0.5){
            mesh.EvS[i] = {0,0,0};
            mesh.EvR[i] = {read.emptyEv[0],read.emptyEv[1],0};
            read.empty[i] = 1;
        }
        else{
            int max = max_element(frac[i].begin(),frac[i].end())-frac[i].begin();
            mesh.EvS[i] = {read.EvS[max][0],read.EvS[max][1],0};
            mesh.EvR[i] = {read.Ev[max][0],read.Ev[max][1],0};
            read.empty[i] = 0;
        }
    }
    file.close();
}

// ----------------------------------------------|
// Sets the free surface list of the elements    |
// ----------------------------------------------|

void surface(readStruct &read,meshStruct &mesh){

    int eLen = mesh.eNode.size();
    mesh.eSurf.resize(eLen);

    for(int i=0; i<eLen; i++){
        if(read.empty[i]==1){continue;}

        else{
            for(int j=0; j<read.neighbour[i].size(); j++){
                
                int idx = read.neighbour[i][j];
                if(j==3 && idx==-1){mesh.eSurf[i].push_back(j);}
                else if(idx!=-1 && read.empty[idx]==1){mesh.eSurf[i].push_back(j);}
            }
        }
    }
}

// --------------------------------------------------|
// Stores the displacement constrains on the nodes   |
// --------------------------------------------------|

void dirichlet(readStruct &read,meshStruct &mesh){

    ivector dLen = read.dLen;
    dvector tol = read.eSize;

    // Tolerance to check if a node is at a boundary

    for(int i=0; i<3; i++){

        dLen[i] *= mesh.order;
        tol[i] = read.eSize[i]/(mesh.order+1);
    }

    // Computes some repeatitly needed values

    vector<ivector> row(3,ivector(mesh.nXYZ.size(),-1));
    ivector add1 = {(dLen[1]+1)*(dLen[2]+1),(dLen[2]+1),1};
    ivector add2 = {dLen[0]*(dLen[1]+1)*(dLen[2]+1),dLen[1]*(dLen[2]+1),dLen[2]};

    // Initialization and lock the node at the origin

    for(int i=0; i<3; i++){
        
        mesh.dirNode[i].push_back(0);
        mesh.dirVal[i].push_back(0);
        mesh.coupNode[i].resize(1);
    }

   // Prevents the elongation of the sheared face

    int f = read.deltaZero.first;
    int d = read.deltaZero.second;
    double fLen = read.zero[f]+read.dSize[f];

    // Selects the nodes at the edge of the conatrained face

    for(int i=0; i<mesh.nXYZ.size(); i++){

        if(abs(mesh.nXYZ[i][f]-fLen)<tol[f]){
            if(abs(mesh.nXYZ[i][d]-read.zero[d])<tol[d]){

                // Change of variable u => Δu = 0 for the other nodes of the face

                for(int j=1; j<=dLen[d]; j++){

                    mesh.deltaNode[d].push_back(make_pair(i+j*add1[d],i));
                    mesh.dirNode[d].push_back(i+j*add1[d]);
                    mesh.dirVal[d].push_back(0);
                }
            }
        }
    }
    
    // Locks the displacement in all directions if clamped

    for(int j:read.clamped){
        for(int i=0; i<mesh.nXYZ.size(); i++){

            if(abs(mesh.nXYZ[i][j]-read.zero[j])<tol[j]){
                for(int k=0; k<3; k++){

                    mesh.dirNode[k].push_back(i);
                    mesh.dirVal[k].push_back(0);
                }
            }
        }
    }

    // Imposes the axial displacement on a top face

    for(pair<int,double> j:read.axial){
        for(int i=0; i<mesh.nXYZ.size(); i++){

            if(abs(mesh.nXYZ[i][j.first]-read.zero[j.first])<tol[j.first]){

                mesh.dirNode[j.first].push_back(i+add2[j.first]);
                mesh.dirVal[j.first].push_back(j.second);
            }
        }
    }

    // Stores the Dirichlet boundary conditions

    for(int j:read.flat){
        for(int i=0; i<mesh.nXYZ.size(); i++){
            if(abs(mesh.nXYZ[i][j]-read.zero[j])<tol[j]){

                // Lock the axis if the current node is at the bottom

                mesh.dirNode[j].push_back(i);
                mesh.dirVal[j].push_back(0);

                // Coupled displacement in the axial dimension

                if(row[j][i+add2[j]]<0){

                    mesh.coupNode[j][0].push_back(i+add2[j]);
                    row[j][i+add2[j]] = 0;
                }
            }
        }
    }
                
    // Coupled displacement in 2 transversal dimensions

    for(int i=0; i<mesh.nXYZ.size(); i++){
        for(pair<int,int> pair:read.coupled){

            int j = pair.first;
            int k = pair.second;
            int loc1 = row[k][i];
            int loc2 = row[k][i+add2[j]];

            if(abs(mesh.nXYZ[i][j]-read.zero[j])<tol[j]){

                // Do not add the nodes with change of variable

                if(abs(mesh.nXYZ[i][f]-fLen)<tol[f]){
                    if(abs(mesh.nXYZ[i][d]-read.zero[d])<tol[d]){
                        continue;
                    }
                }

                // Check whether the first or second node is already in the list

                if(loc1>=0 && loc2<0){

                    mesh.coupNode[k][loc1].push_back(i+add2[j]);
                    row[k][i+add2[j]] = loc1;
                }
                else if(loc1<0 && loc2>=0){

                    mesh.coupNode[k][loc2].push_back(i);
                    row[k][i] = loc2;
                }
                else if(loc1<0 && loc2<0){

                    mesh.coupNode[k].push_back({i,i+add2[j]});
                    row[k][i+add2[j]] = mesh.coupNode[k].size()-1;
                    row[k][i] = mesh.coupNode[k].size()-1;
                }
            }
        }
    }
}

// ------------------------------------------|
// Stores the Neumann boundary conditions    |
// ------------------------------------------|

void neumann(readStruct &read,meshStruct &mesh){

    int idx;
    int order = mesh.order;
    int sLen = mesh.order+1;

    // Neumann BC on the face perpendicular to the j-th dimension

    for(int n=0; n<read.axis[0].size(); n++){
        for(int i=0; i<mesh.eNode.size(); i++){
            if(read.neighbour[i][2*read.axis[0][n]+1]==-1){

                // Computes the nodes of the face without neighbour

                ivector face;

                for(int j=0; j<sLen; j++){
                    for(int k=0; k<sLen; k++){

                        if(read.axis[0][n]==0){idx = j*sLen+k+order*sLen*sLen;}
                        else if(read.axis[0][n]==1){idx = j+k*sLen*sLen+order*sLen;}
                        else if(read.axis[0][n]==2){idx = j*sLen*sLen+k*sLen+order;}
                        face.push_back(mesh.eNode[i][idx]);
                    }
                }

                // Stores the applied axial stress on the right face
                                
                mesh.neuFace.push_back(face);
                mesh.neuVal.push_back("[0,0,0]");
                mesh.neuVal.back()[read.axis[1][n]] = read.Fval;
            }
        }
    }
}

// -------------------------------------------------------|
// Reads the Nascam input files to build the mesh data    |
// -------------------------------------------------------|

void read(string path[2],meshStruct &mesh, timeStruct &time){

    readStruct read;
    readInput(read,mesh,path[0]);
    readMeshSize(read,mesh,path[1]);
    readSpecies(read,mesh,path[1]);

    // Sets the boundary conditions parameters

    if(read.type=="Axial stress"){
        
        read.flat = {0,1,2};
        read.coupled = {make_pair(0,1),make_pair(0,2),make_pair(1,0),make_pair(1,2)};
    }
    else if(read.load=="eXY"){
        
        read.flat = {0,2};
        read.clamped = {0};
        read.deltaZero = make_pair(0,1);
        read.coupled = {make_pair(0,2),make_pair(1,0),make_pair(1,1),make_pair(1,2)};
    }
    else if(read.load=="eYX"){
        
        read.flat = {1,2};
        read.clamped = {1};
        read.deltaZero = make_pair(1,0);
        read.coupled = {make_pair(1,2),make_pair(0,0),make_pair(0,1),make_pair(0,2)};
    }
    else if(read.load=="eZY"){
        
        read.flat = {0,2};
        read.clamped = {2};
        read.deltaZero = make_pair(2,1);
        read.coupled = {make_pair(0,1),make_pair(0,2),make_pair(1,0),make_pair(1,2)};
    }
    else if(read.load=="eZX"){
        
        read.flat = {1,2};
        read.clamped = {2};
        read.deltaZero = make_pair(2,0);
        read.coupled = {make_pair(0,1),make_pair(0,2),make_pair(1,0),make_pair(1,2)};
    }

    else if(read.type=="Hydrostatic"){
        
        read.flat = {0,1,2};
        read.coupled = {make_pair(0,1),make_pair(0,2),make_pair(1,0),make_pair(1,2)};
    }

    // Sets the boundary conditions in the mesh

    dirichlet(read,mesh);
    neumann(read,mesh);
    surface(read,mesh);

    // Rescales the mesh from lattice constant to nanometer

    for(int i=0; i<mesh.nXYZ.size(); i++){
        for(double &n:mesh.nXYZ[i]){n *= read.Lc;}
    }

/*
    cout << "\n\nNodes\n";
    for(int i=0; i<mesh.nXYZ.size(); i++){
        for(int j=0; j<mesh.nXYZ[i].size(); j++){
            cout << mesh.nXYZ[i][j] << ", ";
        }
        cout << "\n";
    }
    
    cout << "\n\nElements\n";
    for(int i=0; i<mesh.eNode.size(); i++){
        for(int j=0; j<mesh.eNode[i].size(); j++){
            cout << mesh.eNode[i][j] << ", ";
        }
        cout << "\n";
    }

    cout << "\n\nFaces\n";
    for(int i=0; i<mesh.neuFace.size(); i++){
        for(int j=0; j<mesh.neuFace[i].size(); j++){
            cout << mesh.neuFace[i][j] << ", ";
        }
        cout << "\n";
    }

    cout << "\n\ndirNode\n";
    for(int i=0; i<3; i++){
        cout << "dim " << i << " -- ";
        for(int j=0; j<mesh.dirNode[i].size(); j++){
            cout << mesh.dirNode[i][j] << ",";
        }
        cout << "\n";
    }

    cout << "\n\ndirVal\n";
    for(int i=0; i<3; i++){
        cout << "dim " << i << " -- ";
        for(int j=0; j<mesh.dirVal[i].size(); j++){
            cout << mesh.dirVal[i][j] << ",";
        }
        cout << "\n";
    }

    cout << "\n\ncoupNode\n";
    for(int i=0; i<3; i++){
        for(int j=0; j<mesh.coupNode[i].size(); j++){
            cout << "dim " << i << " -- ";
            for(int k=0; k<mesh.coupNode[i][j].size(); k++){
                cout << mesh.coupNode[i][j][k] << ",";
            }
            cout << "\n";
        }
        cout << "\n";
    }
    */

}