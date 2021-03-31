#include "..\include\storebc.h"
using namespace std;

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