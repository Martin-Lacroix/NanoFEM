#include "..\include\writer.h"
using namespace std;

// ------------------------------------------------|
// Writes the simulation results in a text file    |
// ------------------------------------------------|

void write(Mesh &mesh,darray &disp,dvector &VM){

    ofstream uXYZ("output/disp.txt");
    ofstream elem("output/elem.txt");
    ofstream node("output/node.txt");
    ofstream stress("output/stress.txt");

    // Writes the displacement field in a text file

    for(int i=0; i<mesh.nLen; i++){
        for(int j=0; j<2; j++){uXYZ << disp[i+j*mesh.nLen] << ",";}
        uXYZ << disp[i+2*mesh.nLen] << "\n";
    }

    // Writes the node coordinates as (x,y,z)

    for(array3d nXYZ:mesh.data.nXYZ){
        for(int i=0; i<nXYZ.size()-1; i++){node << nXYZ[i] << ",";}
        node << nXYZ.back() << "\n";
    }

    // Writes the element nodes as (elem,node)

    for(ivector eNode:mesh.data.eNode){
        for(int i=0; i<eNode.size()-1; i++){elem << eNode[i] << ",";}
        elem << eNode.back() << "\n";
    }

    // Writes the averaged elemental Von Mises stress

    for(int i=0; i<mesh.eLen; i++){
        stress << VM[i] << "\n";
    }
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

void writeJmol(Mesh &mesh,darray &disp,dvector &VM){

    int nLen = mesh.nLen;
    int eLen = mesh.eLen;
    int sLen = mesh.shape3D.N.cols();
    double scale = abs(mesh.data.nXYZ[1][2]-mesh.data.nXYZ[0][2]);

    vector<string> header(18);
    vector<string> uName = {"X","Y","Z"};
    vector<string> sName = {"XX","YY","ZZ","XY","YZ","ZX"};

    // Parameters for the colour legend

    int barLength = 3000;
    double zMin = mesh.data.nXYZ[0][2];
    double zMax = mesh.data.nXYZ.back()[2];
    double xLoc = mesh.data.nXYZ[0][0]-mesh.data.nXYZ.back()[0]/4.0;
    double yLoc = mesh.data.nXYZ[0][1]-mesh.data.nXYZ.back()[1]/4.0;

    // Header parameters for Jmol

    header[1] = "jmolscript: set autobond off";
    header[2] = "set specular OFF";
    header[3] = "set antialiasDisplay OFF";
    header[4] = "set platformSpeed 2";
    header[5] = "spacefill "+to_string(scale/2);
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

    // Writes the Von Mises stresses

    ofstream sXYZ("output/stress.xyz");
    double min = 0;
    double max = 0;

    // Computes the maximum and minimum stress

    for(int i=0; i<eLen; i++){

        if(VM[i]>max){max = VM[i];}
        else if(VM[i]<min){min = VM[i];}
    }

    // Header parameters for Jmol

    header[0] = "Von Mises stress field";
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
        double val = (VM[i]-min)/(max-min);
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
        sXYZ << VM[i] << " ";
        sXYZ << "\n";
    }
}

// ------------------------------------------------|
// Writes the simulation results in a text file    |
// ------------------------------------------------|

void graph(Mesh &mesh,darray &disp,ivector opposite){

    double stress = 0;
    int nLen = mesh.nLen;
    int test = cbrt(nLen)+0.1;
    array3d strain;

    // Computes the total average applied stress

    for(int i=0; i<mesh.fLen; i++){
        stress += math::norm(mesh.data.neuVal[i])/mesh.fLen;
    }

    // Computes the strain along (x,y,z)

    for(int i=0; i<3; i++){

        double u = disp[i*nLen+opposite[i]]-disp[i*nLen];
        double L = mesh.data.nXYZ.back()[i]-mesh.data.nXYZ[0][i];
        strain[i] = abs(u/L)*100;
    }

    // Writes the stress-strain relation

    ofstream graph("output/stress-strain.csv",ios::app);
    graph << stress << ";" << strain[0] << ";" << strain[1] << ";" << strain[2] << "\n";
    graph.close();
}