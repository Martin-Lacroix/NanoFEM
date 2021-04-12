#include "..\include\parser.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>
using namespace std;

// -------------------------------------------|
// Force prints some logs in the terminal     |
// -------------------------------------------|

template <class T>
void logs(string text,T vec,double Lc){

    cout << text << " = [";
    cout << vec[0]*Lc << ", " << vec[1]*Lc << ", " << vec[2]*Lc << "]";
    cout << endl;
}

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
    vector<string> vec;

    // Stores the input parameters

    while(getline(file,input,'\n')){

        transform(input.begin(),input.end(),input.begin(),::tolower);
        vec.push_back(input);
    }

    // Large or small deformation and coating file

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("!simulation type") != string::npos){

            istringstream iss(vec[i]);

            getline(iss,input,';');
            if(input=="large strain"){read.defo = "large";}
            else if(input=="small strain"){read.defo = "small";}

            getline(iss,input,'!');
            input.erase(input.find_last_not_of(" ")+1);
            if(input=="manual selection"){coating = vec[i+1];}
            else if(input=="whole structure"){coating = "input/coating.xyz";}
            break;
        }
    }

    // Reads the order of the quadrature rule

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("!general parameters") != string::npos){

            istringstream iss(vec[i]);
            getline(iss,input,';');
            data.order = stoi(input);
            getline(iss,input,';');
            read.Lc = stod(input);

            if(read.defo=="small"){

                getline(iss,input,'!');
                read.cropZ = stod(input)/read.Lc;
            }
            else if(read.defo=="large"){

                getline(iss,input,';');
                read.cropZ = stod(input)/read.Lc;
                getline(iss,input,'!');
                data.step = stod(input);
            }
            break;
        }
    }

    // Reads the parameters of empty elements

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("!hole parameters") != string::npos){

            istringstream iss(vec[i]);
            getline(iss,input,'!');
            read.emptyLmR = toLame(tovec(input));
            break;
        }
    }

    // Reads the parameters of substrate element

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("!substrate parameters") != string::npos){
            
            istringstream iss(vec[i]);
            getline(iss,input,'!');
            read.LmS.push_back({0,0,0});
            read.LmR.push_back(toLame(tovec(input)));
            break;
        }
    }

    // Reads the type of applied stress

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("!load type") != string::npos){
            
            istringstream iss(vec[i]);
            getline(iss,input,'!');
            input.erase(input.find_last_not_of(" ")+1);
            read.type = input;
            break;
        }
    }

    // Reads the axis, face and value of applied stress

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("!load parameters") != string::npos){
            
            istringstream iss(vec[i]);
            getline(iss,input,';');
            read.load = input;

            for(int i=0; i<2; i++){

                if(input=="all"){read.axis[i] = {0,1,2};}
                else if(input[i+1]=='x'){read.axis[i] = {0};}
                else if(input[i+1]=='y'){read.axis[i] = {1};}
                else if(input[i+1]=='z'){read.axis[i] = {2};}
            }
        
            getline(iss,input,';');
            if(input=="top"){read.free = {1};}
            else if(input=="bottom"){read.free = {0};}
            else if(input=="both"){read.free = {0,1};}
            else if(input=="None"){read.free = {};}

            getline(iss,input,'!');
            read.Fval = stod(input);
            break;
        }
    }

    // Finds the index of the first layer

    int idx;
    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("! layer parameters-layer 1") != string::npos){

            idx = i;
            break;
        }
    }

    // Stores the Young modulus and Poisson ratio of the layers

    for(int i=idx; i<vec.size(); i++){
        if(vec[i].find("! layer parameters-layer") != string::npos){

            istringstream iss(vec[i]);
            getline(iss,input,'!');
            dvector param = tovec(input);
            read.LmR.push_back(toLame({param[0],param[1]}));
            if(read.defo=="small"){read.LmS.push_back({param[2],param[3],param[4]});}
        }
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
            break;
        }
    }

    // Reads the size of the cubic domain

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("SBsize") != string::npos){

            istringstream iss(vec[i]);
            getline(iss,input,'>');
            getline(iss,input);
            read.dSize = tovec(input);
            break;
        }
    }

    // Reads the origin of the domain

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("SBcorner") != string::npos){

            istringstream iss(vec[i]);
            getline(iss,input,'>');
            getline(iss,input);
            read.zero = tovec(input);
            break;
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

    // Prints some logs of the simulation

    cout << "\n\n----------------------\n";
    cout << "File logs";
    cout << "\n----------------------\n\n";
    logs("Element per dimension",read.dLen,1);
    logs("Element size",read.eSize,read.Lc);
    logs("Domain size",read.dSize,read.Lc);
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

string read(string path,dataStruct &data){

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
    
    return read.defo;
}