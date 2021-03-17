#include "..\include\parser.h"
#include <unordered_map>
#include <direct.h>
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

void readInput(readStruct &read,dataStruct &data,string path){

    string input;
    ifstream file;
    file.open(path);

    // Reads the order of the quadrature rule

    getline(file,input,';');
    data.order = stoi(input);
    getline(file,input,';');
    read.Lc = stod(input);
    getline(file,input,'!');
    read.cropZ = stod(input);
    getline(file,input,'\n');
    read.cropZ /= read.Lc;

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
    transform(input.begin(),input.end(),input.begin(),::tolower);
    read.type = input;

    // Reads the axis of applied stress

    getline(file,input,';');
    transform(input.begin(),input.end(),input.begin(),::tolower);
    read.load = input;

    for(int i=0; i<2; i++){

        if(input=="all"){read.axis[i] = {0,1,2};}
        else if(input[i+1]=='x'){read.axis[i] = {0};}
        else if(input[i+1]=='y'){read.axis[i] = {1};}
        else if(input[i+1]=='z'){read.axis[i] = {2};}
    }

    // reads the surface location

    getline(file,input,';');
    transform(input.begin(),input.end(),input.begin(),::tolower);

    if(input=="top"){read.free = {1};}
    else if(input=="bottom"){read.free = {0};}
    else if(input=="both"){read.free = {0,1};}

    // Reads the value of the stress

    getline(file,input,'!');
    read.Fval = stod(input);
    getline(file,input,'\n');

    // Stores the Young modulus and Poisson ratio of the layers

    while(getline(file,input,'!')){

        dvector param = tovec(input);
        read.Ev.push_back({param[0],param[1]});
        read.EvS.push_back({param[2],param[3],param[4]});
        getline(file,input,'\n');
    }
}

// -----------------------------------------------|
// Reads the parameter from the input.xyz file    |
// -----------------------------------------------|

void readMeshSize(readStruct &read,dataStruct &data,string path){

    string input;
    ifstream file;
    file.open(path);
    int order = data.order;

    // Reads the size of the cubic domain

    getline(file,input,'\n');
    getline(file,input,'>');
    getline(file,input,';');
    read.dSize = tovec(input);

    // Reads the origin of the domain

    getline(file,input,'>');
    getline(file,input,';');
    read.zero = tovec(input);
    dvector zero = read.zero;

    // Reads the size of a hexahedron finite element

    getline(file,input,'>');
    getline(file,input,';');
    read.eSize = tovec(input);
    dvector eSize = read.eSize;

    // Truncates the size of the domain to the closest element

    for(int i=0; i<3; i++){
        read.dSize[i] -= fmod(read.dSize[i],eSize[i]);
    }

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

                data.nXYZ.push_back({0,0,0});
                data.nXYZ.back()[0] = zero[0]+i*eSize[0]/order;
                data.nXYZ.back()[1] = zero[1]+j*eSize[1]/order;
                data.nXYZ.back()[2] = zero[2]+k*eSize[2]/order;
            }
        }
    }

    // Builds the elements and the faces of the mesh

    for(int i=0; i<dLen[0]; i++){
        for(int j=0; j<dLen[1]; j++){
            for(int k=0; k<dLen[2]; k++){
                
                data.eNode.push_back({});
                int eID = data.eNode.size()-1;
                read.neighbour.push_back(ivector(6,-1));

                // Computes the first index and geometrical spacing

                int dy = order*dLen[2]+1;
                int dx = (order*dLen[1]+1)*(order*dLen[2]+1);
                int idx = i*(dLen[1]*order+1)*(dLen[2]*order+1)*order+j*(dLen[2]*order+1)*order+k*order;

                // Stores the nodes of each element in the mesh

                for(int n=0; n<order+1; n++){
                    for(int m=0; m<order+1; m++){
                        for(int p=0; p<order+1; p++){

                            data.eNode.back().push_back(idx+p+m*dy+n*dx);
                        }
                    }
                }

                // Stores the list of neightbour elements

                if(k!=0){read.neighbour.back()[0] = eID-1;}
                if(k!=dLen[2]-1){read.neighbour.back()[1] = eID+1;}
                if(j!=0){read.neighbour.back()[2] = eID-dLen[2];}
                if(j!=dLen[1]-1){read.neighbour.back()[3] = eID+dLen[2];}
                if(i!=0){read.neighbour.back()[4] = eID-dLen[2]*dLen[1];}
                if(i!=dLen[0]-1){read.neighbour.back()[5] = eID+dLen[2]*dLen[1];}
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

void readSpecies(readStruct &read,dataStruct &data,string path){

    string input;
    ifstream file;
    file.open(path);
    getline(file,input,'\n');
    getline(file,input,'\n');

    // Initializes the filling fraction of the elements

    int eLen = data.eNode.size();
    vector<dvector> frac(eLen,dvector(read.Ev.size(),0));
    read.empty.resize(eLen);
    data.EvR.resize(eLen);
    data.EvS.resize(eLen);

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
            data.EvS[i] = {0,0,0};
            data.EvR[i] = {read.emptyEv[0],read.emptyEv[1],0};
            read.empty[i] = 1;
        }
        else{
            int max = max_element(frac[i].begin(),frac[i].end())-frac[i].begin();
            data.EvS[i] = {read.EvS[max][0],read.EvS[max][1],read.EvS[max][2]};
            data.EvR[i] = {read.Ev[max][0],read.Ev[max][1],0};
            read.empty[i] = 0;
        }
    }
    file.close();
}

// ----------------------------------------------|
// Sets the free surface list of the elements    |
// ----------------------------------------------|

void surface(readStruct &read,dataStruct &data){

    int eLen = data.eNode.size();
    data.eSurf.resize(eLen);

    for(int i=0; i<eLen; i++){
        if(read.empty[i]==1){continue;}

        else{
            for(int j=0; j<read.neighbour[i].size(); j++){
                
                int idx = read.neighbour[i][j];
                if(idx==-1){for(int k:read.free){if(j==k){data.eSurf[i].push_back(k);}}}
                else if(read.empty[idx]==1){data.eSurf[i].push_back(j);}
            }
        }
    }
}

// --------------------------------------------------|
// Stores the displacement constrains on the nodes   |
// --------------------------------------------------|

void dirichlet(readStruct &read,dataStruct &data){

    ivector dLen = read.dLen;
    dvector tol = read.eSize;

    // Tolerance to check if a node is at a boundary

    for(int i=0; i<3; i++){

        dLen[i] *= data.order;
        tol[i] = read.eSize[i]/(data.order+1);
    }

    // RowLoc stores tracks the node added each row of in coupNode

    vector<ivector> rowLoc(3,ivector(data.nXYZ.size(),-1));
    ivector opposite = {dLen[0]*(dLen[1]+1)*(dLen[2]+1),dLen[1]*(dLen[2]+1),dLen[2]};

    // Initialization and lock the node (0,0,0) at the origin

    for(int i=0; i<3; i++){
        
        data.dirNode[i].push_back(0);
        data.dirVal[i].push_back(0);
        data.coupNode[i].resize(1);
    }

    // Uniform (j,k) : sets the displacement along k of the top face j uniform

    for(pair<int,int> pair:read.uniform){

        int j = pair.first;
        int k = pair.second;

        // For each node at the top surface of the dimension j

        for(int i=0; i<data.nXYZ.size(); i++){
            if(abs(data.nXYZ[i][j]-read.zero[j]-read.dSize[j])<tol[j]){

                // Adds the node to coupNode[k] if not already there

                if(rowLoc[k][i]<0){

                    data.coupNode[k][0].push_back(i);
                    rowLoc[k][i] = 0;
                }
            }
        }
    }

    // Axial (j,k) : imposes the displacement k along j of the top face j

    for(pair<int,double> pair:read.axial){

            int j = pair.first;
            double k = pair.second;

        for(int i=0; i<data.nXYZ.size(); i++){

            // For each node at the top surface of the dimension j

            if(abs(data.nXYZ[i][j]-read.zero[j]-read.dSize[j])<tol[j]){

                data.dirNode[j].push_back(i);
                data.dirVal[j].push_back(k);
            }
        }
    }

    // LockBot (j,k) : locks the displacement along k of the bottom face j at 0

    for(pair<int,int> pair:read.lockBot){

        int j = pair.first;
        int k = pair.second;

        // For each node at the bottom surface of the dimension j

        for(int i=0; i<data.nXYZ.size(); i++){
            if(abs(data.nXYZ[i][j]-read.zero[j])<tol[j]){

                data.dirNode[k].push_back(i);
                data.dirVal[k].push_back(0);
            }
        }
    }

    // LockTop (j,k) : locks the displacement along k of the top face j at 0

    for(pair<int,int> pair:read.lockTop){

        int j = pair.first;
        int k = pair.second;

        // For each node at the top surface of the dimension j

        for(int i=0; i<data.nXYZ.size(); i++){
            if(abs(data.nXYZ[i][j]-read.zero[j]-read.dSize[j])<tol[j]){

                data.dirNode[k].push_back(i);
                data.dirVal[k].push_back(0);
            }
        }
    }

    // Coupled (j,k) : coupled u along k of opposite nodes of the top-bottom faces j

    for(int i=0; i<data.nXYZ.size(); i++){
        for(pair<int,int> pair:read.coupled){

            int j = pair.first;
            int k = pair.second;
            int loc1 = rowLoc[k][i];
            int loc2 = rowLoc[k][i+opposite[j]];

            // For each node at the bottom surface of the dimension j

            if(abs(data.nXYZ[i][j]-read.zero[j])<tol[j]){

                // Check whether the first or second node is already in the list

                if(loc1>=0 && loc2<0){

                    data.coupNode[k][loc1].push_back(i+opposite[j]);
                    rowLoc[k][i+opposite[j]] = loc1;
                }
                else if(loc1<0 && loc2>=0){

                    data.coupNode[k][loc2].push_back(i);
                    rowLoc[k][i] = loc2;
                }
                else if(loc1<0 && loc2<0){

                    data.coupNode[k].push_back({i,i+opposite[j]});
                    rowLoc[k][i+opposite[j]] = data.coupNode[k].size()-1;
                    rowLoc[k][i] = data.coupNode[k].size()-1;
                }
            }
        }
    }

    // Selects the nodes at the edge of the conatrained face

    for(int i=0; i<data.nXYZ.size(); i++){
        for(pair<int,int> pair:read.deltaZero){

            int j = pair.first;
            int k = pair.second;

            // For each node at the bottom surface of the dimension j

            if(abs(data.nXYZ[i][j]-read.zero[j])<tol[j]){

                // Does not take the node if at the surface k

                if(abs(data.nXYZ[i][k]-read.zero[k])>tol[j]){
                    if(abs(data.nXYZ[i][k]-read.zero[k]-read.dSize[k])>tol[j]){

                        // Change of variable u => Δu = 0 for the other nodes of the face

                        data.deltaNode[j].push_back(make_pair(i+opposite[j],i));
                        data.dirNode[j].push_back(i+opposite[j]);
                        data.dirVal[j].push_back(0);
                    }
                }
            }
        }
    }
}

// ------------------------------------------|
// Stores the Neumann boundary conditions    |
// ------------------------------------------|

void neumann(readStruct &read,dataStruct &data){

    int idx;
    int order = data.order;
    int sLen = data.order+1;

    // Neumann BC on the face perpendicular to the j-th dimension

    for(int n=0; n<read.axis[0].size(); n++){
        for(int i=0; i<data.eNode.size(); i++){
            if(read.neighbour[i][5-2*read.axis[0][n]]==-1){

                // Computes the nodes of the face without neighbour

                ivector face;

                for(int j=0; j<sLen; j++){
                    for(int k=0; k<sLen; k++){

                        if(read.axis[0][n]==0){idx = j*sLen+k+order*sLen*sLen;}
                        else if(read.axis[0][n]==1){idx = j+k*sLen*sLen+order*sLen;}
                        else if(read.axis[0][n]==2){idx = j*sLen*sLen+k*sLen+order;}
                        face.push_back(data.eNode[i][idx]);
                    }
                }

                // Stores the applied axial stress on the right face
                                
                data.neuFace.push_back(face);
                data.neuVal.push_back("[0,0,0]");
                data.neuVal.back()[read.axis[1][n]] = read.Fval;
            }
        }
    }
}

// -------------------------------------------------------|
// Reads the Nascam input files to build the mesh data    |
// -------------------------------------------------------|

void read(string path[2],dataStruct &data){

    readStruct read;
    readInput(read,data,path[0]);
    readMeshSize(read,data,path[1]);
    readSpecies(read,data,path[1]);

    // Sets the boundary conditions parameters

    if(read.type=="axial stress"){
        
        read.uniform = {make_pair(0,0),make_pair(1,1),make_pair(2,2)};
        read.lockBot = {make_pair(0,0),make_pair(1,1),make_pair(2,2)};
        read.coupled = {make_pair(0,1),make_pair(0,2),make_pair(1,0),make_pair(1,2)};
    }
    else if(read.type=="shear stress" && read.load=="exy"){

        read.uniform = {make_pair(0,1)};
        read.lockTop = {make_pair(0,0),make_pair(0,2)};
        read.lockBot = {make_pair(0,0),make_pair(0,1),make_pair(0,2),make_pair(2,2)};
        read.coupled = {make_pair(1,0),make_pair(1,2)};
        read.deltaZero = {make_pair(1,0)};
    }
    else if(read.type=="shear stress" && read.load=="eyx"){

        read.uniform = {make_pair(1,0)};
        read.lockTop = {make_pair(1,1),make_pair(1,2)};
        read.lockBot = {make_pair(1,0),make_pair(1,1),make_pair(1,2),make_pair(2,2)};
        read.coupled = {make_pair(0,1),make_pair(0,2)};
        read.deltaZero = {make_pair(0,1)};
    }
    else if(read.type=="shear stress" && read.load=="ezy"){
        
        read.uniform = {make_pair(2,1)};
        read.lockTop = {make_pair(2,0),make_pair(2,2)};
        read.lockBot = {make_pair(2,0),make_pair(2,1),make_pair(2,2)};
        read.coupled = {make_pair(0,1),make_pair(0,2),make_pair(1,0),make_pair(1,2)};
        read.deltaZero = {make_pair(0,2),make_pair(1,2)};
    }
    else if(read.type=="shear stress" && read.load=="ezx"){

        read.uniform = {make_pair(2,0)};
        read.lockTop = {make_pair(2,1),make_pair(2,2)};
        read.lockBot = {make_pair(2,0),make_pair(2,1),make_pair(2,2)};
        read.coupled = {make_pair(0,1),make_pair(0,2),make_pair(1,0),make_pair(1,2)};
        read.deltaZero = {make_pair(0,2),make_pair(1,2)};
    }
    else if(read.type=="hydrostatic" && read.load=="all"){
        
        read.uniform = {make_pair(0,0),make_pair(1,1),make_pair(2,2)};
        read.lockBot = {make_pair(0,0),make_pair(1,1),make_pair(2,2)};
        read.coupled = {make_pair(0,1),make_pair(0,2),make_pair(1,0),make_pair(1,2)};
    }
    else if(read.type=="hydrostatic" && read.load=="exx"){
        
        read.uniform = {make_pair(0,0)};
        read.lockTop = {make_pair(1,1),make_pair(2,2)};
        read.lockBot = {make_pair(0,0),make_pair(0,1),make_pair(0,2),make_pair(1,1),make_pair(2,2)};
        read.coupled = {make_pair(1,0),make_pair(1,2)};
    }
    else if(read.type=="hydrostatic" && read.load=="eyy"){
        
        read.uniform = {make_pair(1,1)};
        read.lockTop = {make_pair(0,0),make_pair(2,2)};
        read.lockBot = {make_pair(0,0),make_pair(1,0),make_pair(1,1),make_pair(1,1),make_pair(2,2)};
        read.coupled = {make_pair(0,1),make_pair(0,2)};
    }
    else if(read.type=="hydrostatic" && read.load=="ezz"){
        
        read.uniform = {make_pair(2,2)};
        read.lockTop = {make_pair(0,0),make_pair(1,1)};
        read.lockBot = {make_pair(0,0),make_pair(1,1),make_pair(2,2),make_pair(2,0),make_pair(2,1)};
        read.coupled = {make_pair(0,1),make_pair(0,2),make_pair(1,0),make_pair(1,2)};
    }

    // Sets the boundary conditions in the data

    dirichlet(read,data);
    neumann(read,data);
    surface(read,data);

    // Rescales the data from lattice constant to nanometer

    for(int i=0; i<data.nXYZ.size(); i++){
        for(double &n:data.nXYZ[i]){n *= read.Lc;}
    }

    /*
    cout << "\n\nNodes\n";
    for(int i=0; i<data.nXYZ.size(); i++){
        for(int j=0; j<data.nXYZ[i].size(); j++){
            cout << data.nXYZ[i][j] << ", ";
        }
        cout << "\n";
    }

    cout << "\n\nElements\n";
    for(int i=0; i<data.eNode.size(); i++){
        for(int j=0; j<data.eNode[i].size(); j++){
            cout << data.eNode[i][j] << ", ";
        }
        cout << "\n";
    }

    cout << "\n\nNeighbour\n";
    for(int i=0; i<read.neighbour.size(); i++){
        for(int j=0; j<read.neighbour[i].size(); j++){
            cout << read.neighbour[i][j] << ", ";
        }
        cout << "\n";
    }

    cout << "\n\nEmpty\n";
    for(int i=0; i<read.empty.size(); i++){
        cout << "Elem " << i << " -- " << read.empty[i] << "\n";
    
    }

    cout << "\n\neSurf\n";
    for(int i=0; i<data.eSurf.size(); i++){
        cout << "Elem " << i << " -- ";
        for(int j=0; j<data.eSurf[i].size(); j++){
            cout << data.eSurf[i][j] << ", ";
        }
        cout << "\n";
    }

    cout << "\n\nFaces\n";
    for(int i=0; i<data.neuFace.size(); i++){
        for(int j=0; j<data.neuFace[i].size(); j++){
            cout << data.neuFace[i][j] << ", ";
        }
        cout << "\n";
    }

    cout << "\n\ndirNode\n";
    for(int i=0; i<3; i++){
        cout << "dim " << i << " -- ";
        for(int j=0; j<data.dirNode[i].size(); j++){
            cout << data.dirNode[i][j] << ",";
        }
        cout << "\n";
    }

    cout << "\n\ndirVal\n";
    for(int i=0; i<3; i++){
        cout << "dim " << i << " -- ";
        for(int j=0; j<data.dirVal[i].size(); j++){
            cout << data.dirVal[i][j] << ",";
        }
        cout << "\n";
    }

    cout << "\n\ncoupNode\n";
    for(int i=0; i<3; i++){
        for(int j=0; j<data.coupNode[i].size(); j++){
            cout << "dim " << i << " -- ";
            for(int k=0; k<data.coupNode[i][j].size(); k++){
                cout << data.coupNode[i][j][k] << ",";
            }
            cout << "\n";
        }
        cout << "\n";
    }

    cout << "\n\ndeltaNode\n";
    for(int i=0; i<3; i++){
        cout << "dim " << i << " -- ";
        for(int j=0; j<data.deltaNode[i].size(); j++){
            cout << "(" << data.deltaNode[i][j].first << "," << data.deltaNode[i][j].second << ") ";
        }
        cout << "\n";
    }
    cout << "\n";
    */
}

// ----------------------------------------------------------------------|
// Converts a normed quantity into a chemical species (Jérôme MULLER)    |
// ----------------------------------------------------------------------|

const char* FCT_atm_name(double norm){

    int i_color;
    int Nb_color = 30;
    i_color = (int)ceil(norm*(Nb_color-1));

    // Selects the chemical species

    switch(i_color){

       case 0 : return("Pb");
       case 1 : return("Ir");
       case 2 : return("Os");
       case 3 : return("Re");
       case 4 : return("Pu");
       case 5 : return("Np");
       case 6 : return("U" );
       case 7 : return("Pa");
       case 8 : return("Th");
       case 9 : return("Ta");
       case 10 : return("Am");
       case 11 : return("Cm");
       case 12 : return("Bk");
       case 13 : return("Cf");
       case 14 : return("Es");
       case 15 : return("Fm");
       case 16 : return("Md");
       case 17 : return("No");
       case 18 : return("Lr");
       case 19 : return("Rf");
       case 20 : return("Db");
       case 21 : return("Sg");
       case 22 : return("Bh");
       case 23 : return("Hs");
       case 24 : return("Mt");
       case 25 : return("Fe");
       case 26 : return("P" );
       case 27 : return("Se");
       case 28 : return("Au");
       case 29 : return("S" );
       default : return("XX");
    }
}

// ------------------------------------------------------------|
// Writes the output displacement file compatible with Jmol    |
// ------------------------------------------------------------|

void writeJmol(Mesh &mesh,darray &disp,vector<darray> &sigma){

    mkdir("output");
    int nLen = mesh.nLen;
    int eLen = mesh.eLen;
    int sLen = mesh.shape3D.N.cols();
    double scale = abs(mesh.data.nXYZ[1][2]-mesh.data.nXYZ[0][2]);

    vector<string> header(18);
    vector<string> uName = {"X","Y","Z"};
    vector<string> sName = {"XX","YY","ZZ","XY","YZ","ZX"};

    // Parameters for the colour legend

    int barLength = 300;
    double zMin = mesh.data.nXYZ[0][2];
    double zMax = mesh.data.nXYZ.back()[2];
    double xLoc = mesh.data.nXYZ[0][0]-mesh.data.nXYZ.back()[0]/4.0;
    double yLoc = mesh.data.nXYZ[0][1]-mesh.data.nXYZ.back()[1]/4.0;

    // Header parameters for Jmol

    header[1] = "jmolscript: set autobond off";
    header[2] = "set specular OFF";
    header[3] = "set antialiasDisplay OFF";
    header[4] = "set platformSpeed 2";
    header[5] = "spacefill "+to_string(scale);
    header[6] = "vector off";
    header[7] = "moveto 0 BOTTOM";
    header[8] = "anim off";
    header[9] = "set labeloffset -10 0";
    header[10] = "set labelfront";
    header[11] = "color label white";
    header[12] = "font label 30 serif";
    header[13] = "select all";

    // Loops in the 3 dimensions (ux, uy, uz)

    for(int k=0; k<3; k++){

        ofstream uXYZ("output/displacement-"+uName[k]+".xyz");
        double min = 0;
        double max = 0;

        // Computes the maximum and minimum displacement

        for(int i=0; i<nLen; i++){

            if(disp[i+k*nLen]>max){max = disp[i+k*nLen];}
            else if(disp[i+k*nLen]<min){min = disp[i+k*nLen];}
        }

        // Header parameters for Jmol

        header[0] = "Displacement field u"+uName[k];
        header[14] = "select atomno="+to_string(1);
        header[15] = "label "+to_string(min)+" nm";
        header[16] = "select atomno="+to_string(barLength);
        header[17] = "label "+to_string(max)+" nm";

        // Writes the header and number of nodes

        uXYZ << nLen+barLength << "\n";
        for(string option:header){uXYZ << option << ";";}
        uXYZ << "\n";

        // Writes the legend colour bar in the file

        for(int i=0; i<barLength; i++){

            double val = i/(double)barLength;
            uXYZ << FCT_atm_name(val) << " ";
            uXYZ << xLoc << " " << yLoc << " " << zMin+zMax*i/(double)barLength;
            uXYZ << "\n";
        }

        // Writes the displacement field in the file

        for(int i=0; i<mesh.nLen; i++){

            double val = (disp[i+k*nLen]-min)/(max-min);
            const char *species = FCT_atm_name(val);
            uXYZ << species << " ";

            // Extracts the u(x,y,z) values from the solution vector

            for(int j=0; j<3; j++){uXYZ << mesh.data.nXYZ[i][j] << " ";}
            uXYZ << disp[i+k*nLen] << " ";
            uXYZ << "\n";
        }
    }

    // Loops in the 6 dimensions (sxx, syy, szz, sxy, syz, szx)

    for(int k=0; k<6; k++){

        ofstream sXYZ("output/stress-"+sName[k]+".xyz");
        double min = 0;
        double max = 0;

        // Computes the maximum and minimum stress

        for(int i=0; i<eLen; i++){

            if(sigma[i][k]>max){max = sigma[i][k];}
            else if(sigma[i][k]<min){min = sigma[i][k];}
        }

        // Header parameters for Jmol

        header[0] = "Displacement field u"+uName[k];
        header[14] = "select atomno="+to_string(1);
        header[15] = "label "+to_string(min)+" GPa";
        header[16] = "select atomno="+to_string(barLength);
        header[17] = "label "+to_string(max)+" GPa";

        // Writes the header and number of nodes

        sXYZ << eLen+barLength << "\n";
        for(string option:header){sXYZ << option << ";";}
        sXYZ << "\n";

        // Writes the legend colour bar in the file

        for(int i=0; i<barLength; i++){

            double val = i/(double)barLength;
            sXYZ << FCT_atm_name(val) << " ";
            sXYZ << xLoc << " " << yLoc << " " << zMin+zMax*i/(double)barLength;
            sXYZ << "\n";
        }

        // Writes the stress field in the file

        for(int i=0; i<mesh.eLen; i++){

            array3d xyz = {0,0,0};
            double val = (sigma[i][k]-min)/(max-min);
            const char *species = FCT_atm_name(val);
            sXYZ << species << " ";


            // Coordinates of the center of the element

            for(int j:mesh.data.eNode[i]){
                for(int n=0; n<3; n++){
                    xyz[n] += mesh.data.nXYZ[j][n]/sLen;
                }
            }

            // Extracts the stress values from the solution vector

            for(int j=0; j<3; j++){sXYZ << xyz[j] << " ";}
            sXYZ << sigma[i][k] << " ";
            sXYZ << "\n";
        }
    }
}