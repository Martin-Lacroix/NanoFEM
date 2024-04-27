#include "mesh.h"
#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
using namespace std;

// -------------------------------------------|
// Sets the starting time before a process    |
// -------------------------------------------|

double start(string text){

    // Gets the current clock time

    cout << text << " ... " << flush;
    auto time = chrono::system_clock::now();
    auto now = chrono::duration<double>(time.time_since_epoch());
    return now.count();
}

// --------------------------------------------|
// Prints the computation time of a process    |
// --------------------------------------------|

void end(double start){

    auto time = chrono::system_clock::now();
    auto now = chrono::duration<double>(time.time_since_epoch());
    double laps = now.count()-start;

    // Gets the elapsed computation time

    int min = laps/60;
    double sec = laps-min*60;

    // Prints the computation time
    
    if(min > 0) cout << "Done in " << min << " min " << setprecision(2) << sec << " sec";
    else cout << "Done in " << setprecision(2) << sec << " sec";
    cout << endl;
}

// ------------------------------------------------|
// Writes the simulation results in a text file    |
// ------------------------------------------------|

void write(Mesh &mesh,darray &disp,dvector &VM){

    mkdir("output",S_IRWXU);
    ofstream uXYZ("output/disp.txt");
    ofstream elem("output/elem.txt");
    ofstream node("output/node.txt");
    ofstream stress("output/stress.txt");

    // Writes the displacement field in a text file

    for(int i = 0; i < mesh.nLen; i++){
        for(int j = 0; j < 2; j++) uXYZ << disp[i+j*mesh.nLen] << ",";
        uXYZ << disp[i+2*mesh.nLen] << "\n";
    }

    // Writes the node coordinates as (x,y,z)

    for(array3d nXYZ : mesh.data.nXYZ){
        for(int i = 0; i < nXYZ.size()-1; i++) node << nXYZ[i] << ",";
        node << nXYZ.back() << "\n";
    }

    // Writes the element nodes as (elem,node)

    for(ivector eNode : mesh.data.eNode){
        for(int i = 0; i < eNode.size()-1; i++) elem << eNode[i] << ",";
        elem << eNode.back() << "\n";
    }

    // Writes the averaged elemental Von Mises stress

    for(int i = 0; i < mesh.eLen; i++){
        stress << VM[i] << "\n";
    }
}

// ----------------------------------------------|
// Solves the linear system for SST/SET model    |
// ----------------------------------------------|

darray solveS(Mesh &mesh){

    double time;
    int nLen = 3*mesh.nLen;
    int size = 9*mesh.eLen*pow(mesh.data.order+1,6)/4;

    // Initializes the solver parameters

    sparse K(nLen,nLen);
    darray B(nLen);
    B.setZero();

    // Builds the system matrix and vector

    time = start("Builds the matrix");
    mesh.totalKB(K,B);
    mesh.neumann(B);
    end(time);

    // Applies boundary conditions to K and B

    time = start("Boundary conditions");

    mesh.delta(K,B);
    mesh.coupling(K,B);
    mesh.dirichlet(K,B);
    end(time);
    
    // Solves the symmetric linear system with Alglib

    K.makeCompressed();
    time = start("Solves the system");
    Eigen::ConjugateGradient<sparse,Eigen::Upper> cg;
    cg.compute(K);
    darray u = cg.solve(B);
    mesh.complete(u);
    end(time);

    // Prints if the solver success or not

    cout << "Output of the solver = " << cg.info();
    cout << "\n" << endl;
    return u;
}

// ----------------------------------------------|
// Solves the non-linear system for SVK model    |
// ----------------------------------------------|

darray solveL(Mesh &mesh){

    double time;
    int max = 30;
    int nLen = 3*mesh.nLen;
    int step = mesh.data.step;
    int size = 9*mesh.eLen*pow(mesh.data.order+1,6)/4;
    vector<darray3d> neuVal = mesh.data.neuVal;

    // Initializes the solver parameters

    darray u(nLen);
    darray du(nLen);
    u.setZero();

    // Performs the load iterations

    for(int i = 0; i < step; i++){
        double norm = 0;
        
        // Updates the surface traction vector

        for(int k = 0; k < 3; k++){
            for(int j = 0; j < neuVal.size(); j++){
                mesh.data.neuVal[j](k) = (i+1)*neuVal[j][k]/step;
            }
        }

        for(int j = 0; j < max; j++){
            
            cout << "[" << i+1 << "/" << step << "] - Newton step ";
            time = start(to_string(j+1));

            // Performs a Newton-Raphson iteration

            sparse K(nLen,nLen);
            darray B(nLen);
            B.setZero();

            // Builds the system matrix and vector

            mesh.totalKT(K,B,u);
            mesh.neumann(B);

            // Applies boundary conditions to K and B

            mesh.delta(K,B);
            mesh.coupling(K,B);
            mesh.dirichlet(K,B);
            
            // Solves the symmetric linear system with Alglib

            K.makeCompressed();
            Eigen::ConjugateGradient<sparse,Eigen::Upper> cg;
            cg.compute(K);
            du = cg.solve(B);

            // Update solution and success state

            mesh.complete(du);
            u = u+du;
            end(time);

            // Check the convergence criteria

            double prec = norm;
            norm = u.norm();
            double rez = abs(norm-prec)/norm;
            cout << "Output of the solver = " << cg.info() << endl;
            cout << "Relative correction = " << rez << endl;
            if(rez < mesh.data.tol) break;
        }

        cout << endl;
    }
    return u;
}

// ------------------------------------------|
// Main code of the finite element solver    |
// ------------------------------------------|

int main(){

    dvector VM;
    double time;
    darray disp;
    dataStruct data;
    string model = "SVK";

    // Reads the input files from Nascam

    time = start("\nReads input data");
    // data = [your input data here] <====== [!]
    end(time);

    // Creates the mesh and solves with conjugate gradient

    Mesh mesh(move(data));
    if(model == "SST/SET") disp = solveS(mesh);
    if(model == "SVK") disp = solveL(mesh);

    // Computes Von Mises stresses and updates the nodes

    time = start("Stress extraction");
    if(model == "SST/SET") VM = mesh.stress(disp,0);
    if(model == "SVK") VM = mesh.stress(disp,1);
    mesh.update(disp);
    end(time);

    // Writes the results in a text file

    time = start("Writes the results");
    write(mesh,disp,VM);
    end(time);
    
    cout << endl;
    return 0;
}