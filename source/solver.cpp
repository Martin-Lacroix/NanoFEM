#include "..\include\parser.h"
#include "..\include\writer.h"
#include "solvers.h"
#include <direct.h>
#include <chrono>
using namespace std;

// -------------------------------------------|
// Sets the starting time before a process    |
// -------------------------------------------|

double start(string text){

    // Gets the current clock time

    cout << text << " ... " << std::flush;
    auto time = std::chrono::system_clock::now();
    auto now = std::chrono::duration<double>(time.time_since_epoch());
    return now.count();
}

// --------------------------------------------|
// Prints the computation time of a process    |
// --------------------------------------------|

void end(double start){

    auto time = std::chrono::system_clock::now();
    auto now = std::chrono::duration<double>(time.time_since_epoch());
    double laps = now.count()-start;

    // Gets the elapsed computation time

    int min = laps/60;
    double sec = laps-min*60;

    // Prints the computation time
    
    if(min>0){cout << "done in " << min << " min " << setprecision(2) << sec << " sec";}
    else{cout << "done in " << setprecision(2) << sec << " sec";}
    cout << endl;
}

// -------------------------------------------------------|
// Solves the linear system in quasistatic equilibrium    |
// -------------------------------------------------------|

darray solve(Mesh &mesh){

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
    alglib::sparseconverttocrs(K);
    end(time);
    
    // Solves the symmetric linear system with Alglib

    time = start("Solves the system");
    alglib::lincgcreate(nLen,state);
    alglib::lincgsolvesparse(state,K,1,B);
    alglib::lincgresults(state,u,rep);
    mesh.complete(u);
    end(time);

    return u;
}

// ------------------------------------------|
// Main code of the finite element solver    |
// ------------------------------------------|

int main(){

    double time;
    dataStruct data;
    string path = "input.txt";
    alglib::setglobalthreading(alglib::parallel);

    // Reads the input files from Nascam

    time = start("\nReads the files");
    read(path,data);
    end(time);

    cout << "\n----------------------\n";
    cout << "FEM algorithm";
    cout << "\n----------------------\n";
    cout << endl;

    // Creates the mesh class and solves with conjugate gradient

    Mesh mesh(move(data));
    darray disp = solve(mesh);

    // Computes Von Mises stresses and updates the nodes

    time = start("Stress extraction");
    vector<darray> sigma = mesh.stress(disp);
    mesh.update(disp);
    end(time);

    // Writes the results in a text file

    time = start("Writes the results");
    writeJmol(mesh,disp,sigma);
    write(mesh,disp,sigma);
    end(time);
    
    cout << endl;
    return 0;
}