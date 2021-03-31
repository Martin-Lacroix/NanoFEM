#include "..\include\parser.h"
#include <unordered_map>
#include <direct.h>
#include <fstream>
using namespace std;

// -------------------------------------------|
// Input = [E,v] to Lamé parameters [λ,μ]     |
// -------------------------------------------|

array3d toLame(dvector input){

    array3d LmX;
    LmX[1] = input[0]/(2*(1+input[1]));
    LmX[0] = input[0]*input[1]/((1+input[1])*(1-2*input[1]));
    LmX[2] = 0;
    return LmX;
}

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

string readInput(readStruct &read,dataStruct &data,string path){

    string input;
    ifstream file;
    string coating;
    file.open(path);

    // Gets the coating input file

    getline(file,input,'\n');
    transform(input.begin(),input.end(),input.begin(),::tolower);

    if(input=="manual selection"){

        getline(file,input,'\n');
        coating = input;
    }
    else if(input=="all layers"){
        coating = "input/coating.xyz";
    }

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
    read.emptyLmR = toLame(tovec(input));
    getline(file,input,'\n');

    // Reads the parameters of substrate element
        
    getline(file,input,'!');
    read.LmS.push_back({0,0,0});
    read.LmR.push_back(toLame(tovec(input)));
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
    else if(input=="None"){read.free = {};}

    // Reads the value of the stress

    getline(file,input,'!');
    read.Fval = stod(input);
    getline(file,input,'\n');

    // Stores the Young modulus and Poisson ratio of the layers

    while(getline(file,input,'!')){

        dvector param = tovec(input);
        read.LmR.push_back(toLame({param[0],param[1]}));
        read.LmS.push_back({param[2],param[3],param[4]});
        getline(file,input,'\n');
    }

    return coating;
}

// -----------------------------------------------|
// Reads the parameter from the input.xyz file    |
// -----------------------------------------------|

void readMeshSize(readStruct &read,dataStruct &data,string path){

    string input;
    ifstream file;
    file.open(path);
    vector<string> vec;
    int order = data.order;

    // Stores the input parameters

    getline(file,input,'\n');
    getline(file,input,'\n');
    istringstream iss(input);
    while(getline(iss,input,';')){vec.push_back(input);}

    // Reads the size of the elements

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("MCsize") != string::npos){

            istringstream iss(vec[i]);
            getline(iss,input,'>');
            getline(iss,input);
            read.eSize = tovec(input);
        }
    }

    // Reads the size of the cubic domain

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("SBsize") != string::npos){

            istringstream iss(vec[i]);
            getline(iss,input,'>');
            getline(iss,input);
            read.dSize = tovec(input);
        }
    }

    // Reads the origin of the domain

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("SBcorner") != string::npos){

            istringstream iss(vec[i]);
            getline(iss,input,'>');
            getline(iss,input);
            read.zero = tovec(input);
        }
    }

    dvector zero = read.zero;
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

    // Stores the node coordinates of the nodes in the mesh

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
        if(dx>=dLen[0] || dx<0){continue;}

        // Slightly moves the species along y axis

        for(int j=0; j<2; j++){

            int dy = dLen[1]*(coord[1]-zero[1])/read.dSize[1]+tol[j];
            if(dy>=dLen[1] || dy<0){continue;}

            // Slightly moves the species along z axis

            for(int k=0; k<2; k++){
                        
                int dz = dLen[2]*(coord[2]-zero[2])/read.dSize[2]+tol[k];
                if(dz<dLen[2] && dz>=0){eList.insert(dx*dLen[1]*dLen[2]+dy*dLen[2]+dz);}
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
    vector<dvector> frac(eLen,dvector(read.LmR.size(),0));
    read.empty.resize(eLen);
    data.LmR.resize(eLen);
    data.LmS.resize(eLen);

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
            data.LmR[i] = read.emptyLmR;
            data.LmS[i] = {0,0,0};
            read.empty[i] = 1;
        }
        else{
            int max = max_element(frac[i].begin(),frac[i].end())-frac[i].begin();
            data.LmR[i] = read.LmR[max];
            data.LmS[i] = read.LmS[max];
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
    ivector dLen = read.dLen;
    data.eSurf.resize(eLen);

    // Parameters for reaching the opposite element in the mesh

    for(int i=0; i<3; i++){dLen[i] -= 1;}
    ivector opposite = {dLen[0]*(dLen[1]+1)*(dLen[2]+1),dLen[1]*(dLen[2]+1),dLen[2]};

    // Sets the free surface of each element

    for(int i=0; i<eLen; i++){
        if(read.empty[i]==1){continue;}

        else{
            for(int j=0; j<read.neighbour[i].size(); j++){

                int idx = read.neighbour[i][j];
                if(idx==-1){
                    
                    // Check if the top and bottom Z surfaces are free

                    if(j==1 || j==0){
                        if(find(read.free.begin(),read.free.end(),j)!=read.free.end()){
                            data.eSurf[i].push_back(j);
                        }
                    }

                    // Checks if the periodicity leads to an empty element
                    
                    else if(j==2 && read.empty[i+opposite[1]]==1){data.eSurf[i].push_back(j);}
                    else if(j==3 && read.empty[i-opposite[1]]==1){data.eSurf[i].push_back(j);}
                    else if(j==4 && read.empty[i+opposite[0]]==1){data.eSurf[i].push_back(j);}
                    else if(j==5 && read.empty[i-opposite[0]]==1){data.eSurf[i].push_back(j);}
                }

                // Checks if neighbour element is empty

                else if(read.empty[idx]==1){data.eSurf[i].push_back(j);}
            }
        }
    }
}

// -------------------------------------------------------|
// Reads the Nascam input files to build the mesh data    |
// -------------------------------------------------------|

void read(string path,dataStruct &data){

    readStruct read;
    string coating = readInput(read,data,path);
    readMeshSize(read,data,coating);
    readSpecies(read,data,coating);

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