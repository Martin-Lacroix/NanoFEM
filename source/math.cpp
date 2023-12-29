#include "math.h"
using namespace std;

namespace math{

    // -----------------------------------------------------|
    // Stiffness tensor D for isotropic linear elasticity   |
    // -----------------------------------------------------|

    matrix stiffness(array3d LmX){

        matrix D(6,6);
        D.setZero();

        // Fills the elements of the stiffness tensor

        D(0,1) = D(0,2) = D(1,2) = LmX[0];
        D(1,0) = D(2,0) = D(2,1) = LmX[0];
        D(3,3) = D(4,4) = D(5,5) = LmX[1];
        D(0,0) = D(1,1) = D(2,2) = 2*LmX[1]+LmX[0];

        return D;
    }

    // ------------------------------------------------|
    // Generates the Gauss-Legendre quadrature rule    |
    // ------------------------------------------------|

    quadStruct legendre(int dim,int order){

        quadStruct quad;
        matrix M(order+1,order+1);
        M.setZero();

        // Generates the three term recurrence coefficients

        for(int i = 0; i < order; i++){

            M(i,i+1) = (i+1)/sqrt(4*(i+1)*(i+1)-1);
            M(i+1,i) = (i+1)/sqrt(4*(i+1)*(i+1)-1);
        }

        Eigen::SelfAdjointEigenSolver<matrix> solver(M);
        darray vec = solver.eigenvectors().col(0);
        darray val = solver.eigenvalues();

        // Quadrature rule for a linear element in 2D/3D

        if(dim == 2) return math::makeQuad2D(vec,val);
        else if(dim == 3) return math::makeQuad3D(vec,val);
        else throw std::runtime_error("Incorrect dimension");
    }

    // -----------------------------------------------|
    // Create and fill the quadrature class in 2D     |
    // -----------------------------------------------|

    quadStruct makeQuad2D(darray &vec, darray &val){

        quadStruct quad;
        int nLen = vec.size();

        for(int i = 0; i < nLen; i++){
            for(int j = 0; j < nLen; j++){

                double W = 4*vec(i)*vec(i)*vec(j)*vec(j);
                quad.gRST.push_back({val(i),val(j)});
                quad.weight.push_back(W);
            }
        }

        quad.gLen = quad.weight.size();
        return quad;
    }

    // -----------------------------------------------|
    // Create and fill the quadrature class in 3D     |
    // -----------------------------------------------|

    quadStruct makeQuad3D(darray &vec, darray &val){

        quadStruct quad;
        int nLen = vec.size();

        for(int i = 0; i < nLen; i++){
            for(int j = 0; j < nLen; j++){
                for(int k = 0; k < nLen; k++){
                    
                    double W = 8*vec(i)*vec(i)*vec(j)*vec(j)*vec(k)*vec(k);
                    quad.gRST.push_back({val(i),val(j),val(k)});
                    quad.weight.push_back(W);
                }
            }
        }

        quad.gLen = quad.weight.size();
        return quad;
    }

    // -----------------------------------------------------------------------|
    // Cleans and stores the non-zero indices of a sparse symmetric matrix    |
    // -----------------------------------------------------------------------|

    vector<ivector> sparsemap(sparse &M){

        int nLen = M.rows();
        vector<ivector> row(3*nLen);
        
        // Stores the non-zero indice locations per row and column

        for (int i = 0; i < M.outerSize(); i++){
            for (sparse::InnerIterator it(M,i); it; ++it){

                int j = it.row();
                row[j].push_back(i);
                if(j != i) row[i].push_back(j);
            }
        }
        return row;
    }

    // ---------------------------------------------------------------|
    // Sets the coordinate (row,col) for a symmetric sparse matrix    |
    // ---------------------------------------------------------------|

    void symset(sparse &M,int row,int col,double val){

        if(row > col) M.coeffRef(col,row) = val;
        else M.coeffRef(row,col) = val;
    }

    // ------------------------------------------------------------------|
    // Adds to the coordinate (row,col) for a symmetric sparse matrix    |
    // ------------------------------------------------------------------|

    void symadd(sparse &M,int row,int col,double val){

        if(row > col) M.coeffRef(col,row) += val;
        else M.coeffRef(row,col) += val;
    }

    // ------------------------------------------------------------------------|
    // gets the value at coordinate (row,col) for a symmetric sparse matrix    |
    // ------------------------------------------------------------------------|

    double get(sparse &M,int row,int col){

        double val;
        if(row > col) val = M.coeff(col,row);
        else val = M.coeff(row,col);
        return val;
    }

    // ---------------------------------------------------------------------|
    // Evaluates the shape functions for Lagrange isoparametric elements    |
    // ---------------------------------------------------------------------|

    dvector lagrange(int dvar,dvector node,dvector val){

        int dim = val.size();
        int nLen = pow(node.size(),dim);
        int sLen = node.size();
        dvector N(nLen,1);

        // Var is the index of the variable we want to derivate

        for(int n = 0; n < dim; n++){
            dvector N1;

            // Computes the derivative of Lagrange polynomial for this variable

            if(n == dvar){
                N1.resize(sLen,0);

                for(int j = 0; j < sLen; j++){
                    for(int i = 0; i < sLen; i++){

                        // Computes the first term of the summ of products

                        if(i != j){
                            double La = 1/(node[j]-node[i]);

                            for(int k = 0; k < sLen; k++){
                                if((k != j) && (k != i)) La *= (val[n]-node[k])/(node[j]-node[k]);
                            }
                            N1[j] += La;
                        }
                    }
                }
            }
   
            // Computes the Lagrange polynomial for other variables

            else{
                N1.resize(sLen,1);

                for(int i = 0; i < sLen; i++){
                    for(int j = 0; j < sLen; j++){
                        if(i != j) N1[j] *= (val[n]-node[i])/(node[j]-node[i]);
                    }
                }
            }

            // Product of the Lagrange polynomials in the 2D case

            if(dim == 2){
                for(int i = 0; i < sLen; i++){
                    for(int j = 0; j < sLen; j++){
                            
                        ivector loop = {i,j};
                        N[i*sLen+j] *= N1[loop[n]];
                    }
                }
            }

            // Product of the Lagrange polynomials in the 3D case

            if(dim == 3){
                for(int i = 0; i < sLen; i++){
                    for(int j = 0; j < sLen; j++){
                        for(int k = 0; k < sLen; k++){
                            
                            ivector loop = {i,j,k};
                            N[i*sLen*sLen+j*sLen+k] *= N1[loop[n]];
                        }
                    }
                }
            }
        }
        return N;
    }

    // --------------------------------------------------------------|
    // Projection matrix for the surface gradient in a 3D element    |
    // --------------------------------------------------------------|

    matrix projection(darray3d norm){

        matrix3d P;
        matrix T(6,6);

        // Computes the surface gredient ∇s = (I-n⊗n)∇

        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                
                if(j == k) P(j,k) = 1-norm(j)*norm(k);
                else P(j,k) = -norm(j)*norm(k);
            }
        }

        // Computes the projection matrix with the normal

        T(0,0) = P(0,0)*P(0,0);
        T(0,1) = P(0,1)*P(0,1);
        T(0,2) = P(2,0)*P(2,0);
        T(0,3) = 2*P(0,0)*P(0,1);
        T(0,4) = 2*P(2,0)*P(0,1);
        T(0,5) = 2*P(0,0)*P(2,0);
        
        T(1,0) = P(0,1)*P(0,1);
        T(1,1) = P(1,1)*P(1,1);
        T(1,2) = P(1,2)*P(1,2);
        T(1,3) = 2*P(0,1)*P(1,1);
        T(1,4) = 2*P(1,2)*P(1,1);
        T(1,5) = 2*P(0,1)*P(1,2);
        
        T(2,0) = P(2,0)*P(2,0);
        T(2,1) = P(1,2)*P(1,2);
        T(2,2) = P(2,2)*P(2,2);
        T(2,3) = 2*P(1,2)*P(2,0);
        T(2,4) = 2*P(1,2)*P(2,2);
        T(2,5) = 2*P(2,2)*P(2,0);
        
        T(3,0) = P(0,0)*P(0,1);
        T(3,1) = P(0,1)*P(1,1);
        T(3,2) = P(2,0)*P(1,2);
        T(3,3) = P(0,0)*P(1,1)+P(0,1)*P(0,1);
        T(3,4) = P(2,0)*P(1,1)+P(1,2)*P(0,1);
        T(3,5) = P(2,0)*P(0,1)+P(0,0)*P(1,2);
        
        T(4,0) = P(0,1)*P(2,0);
        T(4,1) = P(1,2)*P(1,1);
        T(4,2) = P(1,2)*P(2,2);
        T(4,3) = P(1,1)*P(2,0)+P(0,1)*P(1,2);
        T(4,4) = P(1,1)*P(2,2)+P(1,2)*P(1,2);
        T(4,5) = P(1,2)*P(2,0)+P(0,1)*P(2,2);
        
        T(5,0) = P(0,0)*P(2,0);
        T(5,1) = P(0,1)*P(1,2);
        T(5,2) = P(2,2)*P(2,0);
        T(5,3) = P(0,1)*P(2,0)+P(0,0)*P(1,2);
        T(5,4) = P(2,0)*P(1,2)+P(2,2)*P(0,1);
        T(5,5) = P(2,2)*P(0,0)+P(2,0)*P(2,0);

        return T;
    }
}
