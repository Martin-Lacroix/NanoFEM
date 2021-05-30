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
        input.erase(remove_if(input.begin(),input.end(),::isspace),input.end());
        vec.push_back(input);
    }

    // Material model and coating file

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("!simulationtype") != string::npos){

            istringstream iss(vec[i]);
            getline(iss,input,';');
            read.model = input;

            getline(iss,input,'!');
            if(input=="manualselection"){coating = vec[i+1];}
            else if(input=="wholestructure"){coating = "input/coating.xyz";}
            break;
        }
    }

    // Reads the general parameters

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("!generalparameters") != string::npos){

            istringstream iss(vec[i]);
            getline(iss,input,'!');
            dvector inp = tovec(input);

            read.Lc = abs(inp[0]);
            read.cropZ = abs(inp[1]);

            if(read.model=="saintvenant-kirchhoff"){

                data.step = abs(inp[2]+1e-5);
                data.tol = abs(inp[3]);
            }
            break;
        }
    }

    // Reads the parameters of Lagrange hexahedrons

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("!finiteelements") != string::npos){

            istringstream iss(vec[i]);
            getline(iss,input,'!');
            dvector inp = tovec(input);

            data.order = abs(inp[0]+1e-5);
            for(int j=1; j<4; j++){read.eSize.push_back(abs(inp[j]));}
            break;
        }
    }

    // Reads the parameters of empty elements

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("!holeparameters") != string::npos){

            istringstream iss(vec[i]);
            getline(iss,input,'!');
            dvector inp = tovec(input);

            read.emptyLmR = toLame({inp[0],inp[1]});
            read.frac = inp[2];
            break;
        }
    }

    // Reads the parameters of substrate element

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("!substrateparameters") != string::npos){
            
            istringstream iss(vec[i]);
            getline(iss,input,'!');
            read.LmS.push_back({0,0,0});
            read.LmR.push_back(toLame(tovec(input)));
            break;
        }
    }

    // Reads the boundary consitions

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("!boundaryconditions") != string::npos){
            
            istringstream iss(vec[i]);

            // Reads the axis of the load

            getline(iss,input,';');
            read.axis = input;

            // Position of the free surface

            getline(iss,input,';');
            if(input=="top"){read.free = {1};}
            else if(input=="bottom"){read.free = {0};}
            else if(input=="both"){read.free = {0,1};}

            // Substrate deformation

            getline(iss,input,';');
            read.transverse = input;

            // Value of the applied stress

            getline(iss,input,'!');
            read.Fval = stod(input);
            break;
        }
    }

    // Finds the index of the first layer

    int idx;
    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("!layerparameters-layer1") != string::npos){

            idx = i;
            break;
        }
    }

    // Stores the Young modulus and Poisson ratio of the layers

    for(int i=idx; i<vec.size(); i++){
        if(vec[i].find("!layerparameters-layer") != string::npos){

            istringstream iss(vec[i]);
            getline(iss,input,'!');
            dvector inp = tovec(input);
            read.LmR.push_back(toLame({inp[0],inp[1]}));

            // Adds the surface parameter in SET/SST

            if(read.model=="linearset/sstmodel"){
                read.LmS.push_back({inp[2],inp[3],inp[4]});
            }
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

    // Reads the initial size of the hexahedrons

    for(int i=0; i<vec.size(); i++){
        if(vec[i].find("MCsize") != string::npos){

            istringstream iss(vec[i]);
            getline(iss,input,'>');
            getline(iss,input);

            read.scale = 1;
            dvector eSize = tovec(input);

            // Makes the size at least equal to Nascam resolution

            for(int j=0; j<3; j++){
                
                if(read.eSize[j]<eSize[j]){read.eSize[j] = eSize[j];}
                read.scale *= abs(eSize[j]/read.eSize[j]);
            }
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
        read.dSize[i] += read.eSize[i]/5;
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

    // Computes the node offset between opposed faces

    for(int i=0; i<3; i++){dLen[i] *= data.order;}
    read.opposite = {dLen[0]*(dLen[1]+1)*(dLen[2]+1),dLen[1]*(dLen[2]+1),dLen[2]};
    file.close();

    // Prints some logs of the simulation

    cout << "\n\n----------------------\n";
    cout << "File logs";
    cout << "\n----------------------\n\n";
    logs("Element per dimension",read.dLen,1);
    logs("Element size (nm)",read.eSize,read.Lc);
    logs("Domain size (nm)",read.dSize,read.Lc);
}

// -------------------------------------------------------|
// Localizes the chemical species in the mesh elements    |
// -------------------------------------------------------|

unordered_set<int> locSpecies(readStruct &read,dvector coord){

    double d;
    unordered_set<int> eList;
    ivector dLen = read.dLen;
    dvector zero = read.zero;
    dvector tol = {-1e-5,1e-5};

    for(int i=0; i<2; i++){

        // Slightly moves the species along x axis

        d = dLen[0]*(coord[0]-zero[0])/read.dSize[0]+tol[i];
        if(d<0){d -= 1;} int dx = static_cast<int>(d);
        if(dx>=dLen[0] || dx<0){continue;}

        // Slightly moves the species along y axis

        for(int j=0; j<2; j++){

            d = dLen[1]*(coord[1]-zero[1])/read.dSize[1]+tol[j];
            if(d<0){d -= 1;} int dy = static_cast<int>(d);
            if(dy>=dLen[1] || dy<0){continue;}

            // Slightly moves the species along z axis

            for(int k=0; k<2; k++){

                d = dLen[2]*(coord[2]-zero[2])/read.dSize[2]+tol[k];
                if(d<0){d -= 1;} int dz = static_cast<int>(d);
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
        for (int i:eList){frac[i][layer] += read.scale*stod(input)/eList.size();}
    }

    // Stores the mixed mechanical parameters of the elements

    for(int i=0; i<eLen; i++){
        double sum = accumulate(frac[i].begin(),frac[i].end(),0.0);

        if(sum<read.frac){
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

// --------------------------------------------------|
// Stores the displacement constrains on the nodes   |
// --------------------------------------------------|

void dirichlet(readStruct &read,dataStruct &data){

    // Tolerance to check if a node is at a boundary

    dvector tol = read.eSize;
    vector<ivector> rowLoc(3,ivector(data.nXYZ.size(),-1));
    for(int i=0; i<3; i++){tol[i] = read.eSize[i]/(data.order+1);}

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

    // LockBot (j,k) : locks the displacement along k of the bottom face j

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

    // LockTop (j,k) : locks the displacement along k of the bottom face j

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

            // For each node at the bottom surface of the dimension j

            if(abs(data.nXYZ[i][j]-read.zero[j])<tol[j]){

                int loc1 = rowLoc[k][i];
                int loc2 = rowLoc[k][i+read.opposite[j]];

                // Check whether the first or second node is already in the list

                if(loc1>=0 && loc2<0){

                    data.coupNode[k][loc1].push_back(i+read.opposite[j]);
                    rowLoc[k][i+read.opposite[j]] = loc1;
                }
                else if(loc1<0 && loc2>=0){

                    data.coupNode[k][loc2].push_back(i);
                    rowLoc[k][i] = loc2;
                }
                else if(loc1<0 && loc2<0){

                    data.coupNode[k].push_back({i,i+read.opposite[j]});
                    rowLoc[k][i+read.opposite[j]] = data.coupNode[k].size()-1;
                    rowLoc[k][i] = data.coupNode[k].size()-1;
                }
            }
        }
    }

    for(int i=0; i<3; i++){
        for(int j=0; j<data.coupNode[i].size(); j++){
            reverse(data.coupNode[i][j].begin(),data.coupNode[i][j].end());
        }
    }

    // Delta (j,k) : set (ut,ub) => (ut-ub,ub) of opposite faces along the j-axis

    for(int i=0; i<data.nXYZ.size(); i++){
        for(int j:read.delta){

            // Change of variable u => Δu for the other nodes of the face

            if(abs(data.nXYZ[i][j]-read.zero[j])<tol[j]){
                data.deltaNode[j].push_back(make_pair(i+read.opposite[j],i));
            }
        }
    }
}

// ------------------------------------------|
// Stores the Neumann boundary conditions    |
// ------------------------------------------|

void neumann(readStruct &read,dataStruct &data){

    int order = data.order;
    int sLen = data.order+1;

    // Axial stress along the z-axis

    if(read.axis=="z-axis"){
        for(int i=0; i<data.eNode.size(); i++){
            if(read.neighbour[i][1]==-1){

                // Computes the nodes of the face without neighbour

                ivector face;
                for(int j=0; j<sLen; j++){
                    for(int k=0; k<sLen; k++){

                        int idx = j*sLen*sLen+k*sLen+order;
                        face.push_back(data.eNode[i][idx]);
                    }
                }

                // Stores the applied axial stress on the right face
                                
                data.neuFace.push_back(face);
                data.neuVal.push_back("[0,0,0]");
                data.neuVal.back()[2] = read.Fval;
            }
        }
    }

    // Axial stress along the y-axis

    if(read.axis=="y-axis"){
        for(int i=0; i<data.eNode.size(); i++){
            if(read.neighbour[i][3]==-1){

                // Computes the nodes of the face without neighbour

                ivector face1;
                ivector face2;

                for(int j=0; j<sLen; j++){
                    for(int k=0; k<sLen; k++){

                        int idx = j+k*sLen*sLen+order*sLen;
                        face1.push_back(data.eNode[i][idx]);
                        face2.push_back(data.eNode[i][idx]-read.opposite[1]);
                    }
                }

                // Stores the applied axial stress on the right face
                                
                data.neuFace.push_back(face1);
                data.neuVal.push_back("[0,0,0]");
                data.neuVal.back()[1] = read.Fval;

                // Stores the applied axial stress on the left face
                                
                data.neuFace.push_back(face2);
                data.neuVal.push_back("[0,0,0]");
                data.neuVal.back()[1] = -read.Fval;
            }
        }
    }

    // Axial stress along the x-axis

    if(read.axis=="x-axis"){
        for(int i=0; i<data.eNode.size(); i++){
            if(read.neighbour[i][5]==-1){

                // Computes the nodes of the face without neighbour

                ivector face1;
                ivector face2;

                for(int j=0; j<sLen; j++){
                    for(int k=0; k<sLen; k++){

                        int idx = j*sLen+k+order*sLen*sLen;
                        face1.push_back(data.eNode[i][idx]);
                        face2.push_back(data.eNode[i][idx]-read.opposite[0]);
                    }
                }

                // Stores the applied axial stress on the top face
                                
                data.neuFace.push_back(face1);
                data.neuVal.push_back("[0,0,0]");
                data.neuVal.back()[0] = read.Fval;

                // Stores the applied axial stress on the bottom face
                                
                data.neuFace.push_back(face2);
                data.neuVal.push_back("[0,0,0]");
                data.neuVal.back()[0] = -read.Fval;
            }
        }
    }
}

// -------------------------------------------------------|
// Reads the Nascam input files to build the mesh data    |
// -------------------------------------------------------|

readStruct reader(string path,dataStruct &data){

    readStruct read;
    string coating = readInput(read,data,path);
    readMeshSize(read,data,coating);
    readSpecies(read,data,coating);

    // Sets the boundary conditions parameters

    if(read.transverse=="periodicstrain"){

        read.delta = {0,1};
        read.lockBot = {make_pair(2,2)};
        read.uniform = {make_pair(0,0),make_pair(1,1),make_pair(2,2)};
        read.coupled = {make_pair(0,1),make_pair(0,2),make_pair(1,0),make_pair(1,2)};
    }
    if(read.transverse=="periodic(nostrain)"){

        read.lockBot = {make_pair(2,2)};
        read.uniform = {make_pair(2,2)};
        read.coupled = {make_pair(0,1),make_pair(0,2),
        make_pair(1,0),make_pair(1,2),make_pair(0,0),make_pair(1,1)};
    }
    if(read.transverse=="uniformstrain"){

        read.lockBot = {make_pair(0,0),make_pair(1,1),make_pair(2,2)};
        read.uniform = {make_pair(0,0),make_pair(1,1),make_pair(2,2)};
        read.coupled = {make_pair(0,1),make_pair(0,2),make_pair(1,0),make_pair(1,2)};
    }
    if(read.transverse=="uniform(nostrain)"){

        read.lockBot = {make_pair(0,0),make_pair(1,1),make_pair(2,2)};
        read.uniform = {make_pair(0,0),make_pair(1,1),make_pair(2,2)};
        read.lockTop = {make_pair(0,0),make_pair(1,1)};
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

    cout << "\n\nForce\n";
    for(int i=0; i<data.neuVal.size(); i++){
        for(int j=0; j<data.neuVal[i].length(); j++){
            cout << data.neuVal[i][j] << ", ";
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


    
    return read;
}