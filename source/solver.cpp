#include "..\include\parser.h"
#include "solvers.h"
#include <iomanip>
#include <chrono>
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
    
    if(min>0){cout << "Done in " << min << " min " << setprecision(2) << sec << " sec";}
    else{cout << "Done in " << setprecision(2) << sec << " sec";}
    cout << endl;
}

// ---------------------------------------------------|
// Solves the linear system for small deformations    |
// ---------------------------------------------------|

darray solveS(Mesh &mesh,ivector opposite){

    double time;
    int nLen = 3*mesh.nLen;
    int size = 9*mesh.eLen*pow(mesh.data.order+1,6)/4;

    // Initializes the solver parameters

    sparse K;
    darray B,u;
    u.setlength(nLen);
    B.setlength(nLen);
    math::zero(B);

    alglib::lincgreport rep;
    alglib::lincgstate state;
    alglib::sparsecreate(nLen,nLen,size,K);

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

    time = start("Solves the system");
    alglib::sparseconverttocrs(K);
    alglib::lincgcreate(nLen,state);
    alglib::lincgsolvesparse(state,K,1,B);
    alglib::lincgresults(state,u,rep);
    mesh.complete(u);
    end(time);

    graph(mesh,u,opposite);

    return u;
}

// -------------------------------------------------------|
// Solves the non-linear system for large deformations    |
// -------------------------------------------------------|

darray solveL(Mesh &mesh,ivector opposite){

    double time;
    int max = 30;
    int nLen = 3*mesh.nLen;
    int step = mesh.data.step;
    int size = 9*mesh.eLen*pow(mesh.data.order+1,6)/4;
    vector<darray> neuVal = mesh.data.neuVal;

    // Initializes the solver parameters

    darray u,du;
    u.setlength(nLen);
    du.setlength(nLen);
    math::zero(u);

    // Performs the load iterations

    for(int i=0; i<step; i++){
        double norm1,norm2=1;

        // Updates the surface traction vector

        for(int k=0; k<3; k++){
            for(int j=0; j<neuVal.size(); j++){
                mesh.data.neuVal[j](k) = (i+1)*neuVal[j][k]/step;
            }
        }

        for(int j=0; j<max; j++){
            
            cout << "[" << i+1 << "/" << step << "] - Newton step ";
            time = start(to_string(j+1));

            // Performs a Newton-Raphson iteration

            sparse K;
            darray B;
            B.setlength(nLen);
            math::zero(B);

            alglib::lincgreport rep;
            alglib::lincgstate state;
            alglib::sparsecreate(nLen,nLen,size,K);

            // Builds the system matrix and vector

            mesh.totalKT(K,B,u);
            mesh.neumann(B);

            // Applies boundary conditions to K and B

            mesh.delta(K,B);
            mesh.coupling(K,B);
            mesh.dirichlet(K,B);
            
            // Solves the symmetric linear system with Alglib

            alglib::sparseconverttocrs(K);
            alglib::lincgcreate(nLen,state);
            alglib::lincgsolvesparse(state,K,1,B);
            alglib::lincgresults(state,du,rep);

            mesh.complete(du);
            math::add(1,1,du,u);
            end(time);

            // Check the convergence criteria

            norm1 = norm2;
            norm2 = math::norm(u);
            double rez = abs(norm2-norm1)/norm1;
            cout << "Relative Correction = " << rez << endl;
            if(rez<mesh.data.tol){break;}
        }

        graph(mesh,u,opposite);
        cout << endl;
    }
    return u;
}

// ------------------------------------------|
// Main code of the finite element solver    |
// ------------------------------------------|

int main(){

    double time;
    darray disp;
    string type;
    dataStruct data;
    ivector opposite;
    string path = "input.txt";
    alglib::setglobalthreading(alglib::parallel);

    // Reads the input files from Nascam

    {
        time = start("\nReads the files");
        readStruct param = read(path,data);
        opposite = param.opposite;
        type = param.deformation;
        end(time);
    }

    // Writes some outputs and logs

    cout << "\n----------------------\n";
    cout << "FEM algorithm";
    cout << "\n----------------------\n";
    cout << endl;

    mkdir("output");
    ofstream graph("output/stress-strain.csv");
    graph << "Absolute stress (GPa);Absolute strain X (%);Absolute strain Y (%);Absolute strain Z (%)\n";
    graph << 0 << ";" << 0 << ";" << 0 << ";" << 0 << "\n";
    graph.close();

    // Creates the mesh and solves with conjugate gradient

    Mesh mesh(move(data));
    if(type=="small"){disp = solveS(mesh,opposite);}
    if(type=="large"){disp = solveL(mesh,opposite);}

    // Computes Von Mises stresses and updates the nodes

    time = start("Stress extraction");
    dvector VM = mesh.stress(disp);
    mesh.update(disp);
    end(time);

    // Writes the results in a text file

    time = start("Writes the results");
    writeJmol(mesh,disp,VM);
    write(mesh,disp,VM);
    end(time);
    
    cout << endl;

    
/*
    cout << "\n";
    for(int i=0; i<mesh.nLen; i++){
        cout << setprecision(7) << "Node " << i << " -- ux = " << disp[i] << "\n";
    }
    cout << "\n";
    for(int i=0; i<mesh.nLen; i++){
        cout << setprecision(7) << "Node " << i << " -- uy = " << disp[i+mesh.nLen] << "\n";
    }
    cout << "\n";
    for(int i=0; i<mesh.nLen; i++){
        cout << setprecision(7) << "Node " << i << " -- uz = " << disp[i+2*mesh.nLen] << "\n";
    }
    cout << "\n";
    */
    
   
    return 0;
}