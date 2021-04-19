#include "..\include\math.h"
using namespace std;

namespace math{

    // -----------------------------------------------------|
    // Stiffness tensor D for isotropic linear elasticity   |
    // -----------------------------------------------------|

    matrix stiffness(array3d LmX){

        matrix D;
        D.setlength(6,6);
        math::zero(D);

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

        matrix M;
        quadStruct quad;
        int nbr = order+1;
        darray Ja; Ja.setlength(nbr);
        darray Jb; Jb.setlength(nbr-1);

        // Generates the three term recurrence coefficients
        
        for(int i=1; i<nbr; i++){Jb(i-1) = i/sqrt(4*i*i-1);}
        for(int i=0; i<nbr; i++){Ja(i) = 0;}
        alglib::smatrixtdevd(Ja,Jb,nbr,3,M);

        // Quadrature rule for 8-node brick linear element in 3D

        if(dim==3){     
            for(int i=0; i<nbr; i++){
                for(int j=0; j<nbr; j++){
                    for(int k=0; k<nbr; k++){

                        quad.gRST.push_back({Ja(i),Ja(j),Ja(k)});
                        quad.weight.push_back(8*M(0,i)*M(0,i)*M(0,j)*M(0,j)*M(0,k)*M(0,k));
                    }
                }
            }
        }

        // Quadrature rule for 4-node quadrangle linear element in 2D

        if(dim==2){
            for(int i=0; i<nbr; i++){
                for(int j=0; j<nbr; j++){

                    quad.gRST.push_back({Ja(i),Ja(j)});
                    quad.weight.push_back(4*M(0,i)*M(0,i)*M(0,j)*M(0,j));
                }
            }
        }

        quad.gLen = quad.weight.size();
        return quad;
    }

    // --------------------------------------------------|
    // Performs the vector cross product V3 = V1 × V2    |
    // --------------------------------------------------|

    array3d cross(array3d &V1,array3d &V2){

        array3d V3;
        V3[0] = V1[1]*V2[2]-V1[2]*V2[1];
        V3[1] = V1[0]*V2[2]-V1[2]*V2[0];
        V3[2] = V1[0]*V2[1]-V1[1]*V2[0];
        return V3;
    }

    // ------------------------------------------------------|
    // Computes the vector operation V3 = V1 - (V1·V2) V2    |
    // ------------------------------------------------------|

    array3d dotsub(array3d &V1,array3d &V2){

        array3d V3;
        double k = 0;

        for(int i=0; i<3; i++){k += V1[i]*V2[i];}
        for(int i=0; i<3; i++){V3[i] = V2[i]-k*V1[i];}
        return V3;
    }

    // ------------------------------------------------------|
    // Standard matrix-matrix addition M2 = k1 M1 + k2 M2    |
    // ------------------------------------------------------|

    void add(double k1,double k2,matrix &M1,matrix &M2){

        int m = M1.rows();
        int n = M1.cols();
        alglib::rmatrixgencopy(m,n,k1,M1,0,0,k2,M2,0,0);
    }

    // ------------------------------------------------------|
    // Standard vector-vector addition V2 = k1 V1 + k2 V2    |
    // ------------------------------------------------------|

    void add(double k1,double k2,darray &V1,darray &V2){

        for(int i=0; i<V1.length(); i++){
            V2(i) = k1*V1(i) + k2*V2(i);
        }
    }

    // --------------------------------------------------------------------------|
    // Matrix-matrix product M3 = k M1 M2, tj = 1 means that Mj is transposed    |
    // --------------------------------------------------------------------------|

    matrix prod(double k,matrix &M1,matrix &M2,int t1,int t2){

        int m = M1.rows();
        int n = M2.cols();
        int w = M1.cols();

        // Transpose the matrices if required

        if(t1==1){w = M1.rows(); m = M1.cols();}
        if(t2==1){n = M2.rows();}

        // Performs the product with Alglib

        matrix M3; M3.setlength(m,n);
        alglib::rmatrixgemm(m,n,w,k,M1,0,0,t1,M2,0,0,t2,0,M3,0,0);
        return M3;
    }

    // --------------------------------------------------------------------------|
    // Matrix-vector product V2 = k M1 V1, t1 = 1 means that M1 is transposed    |
    // --------------------------------------------------------------------------|

    darray prod(double k,matrix &M1,darray &V1,int t1){

        int m = M1.rows();
        int n = M1.cols();

        // Transpose the matrices if required

        if(t1==1){
            n = M1.rows();
            m = M1.cols();
        }

        // Performs the product with Alglib

        darray V2; V2.setlength(m);
        alglib::rmatrixgemv(m,n,k,M1,0,0,t1,V1,0,0,V2,0);
        return V2;
    }

    // --------------------------------------------------|
    // Fills an initialized dense matrix M with zeros    |
    // --------------------------------------------------|

    void zero(matrix &M){

        for(int i=0; i<M.rows(); i++){
            for(int j=0; j<M.cols(); j++){
                M(i,j) = 0;
            }
        }
    }

    // -------------------------------------------------|
    // Fills an initialized dense array V with zeros    |
    // -------------------------------------------------|

    void zero(darray &V){
        for(int i=0; i<V.length(); i++){V(i) = 0;}
    }

    // ------------------------------------------------|
    // Computes the inverse of a 3×3 general matrix    |
    // ------------------------------------------------|

    matrix invert(matrix &M,double det){

        matrix invM;
        invM.setlength(3,3);

        for(int i=0; i<3; i++){
            for(int j=0; j<3; j++){

                invM(i,j) = M[(i+1)%3][(j+1)%3]*M[(i+2)%3][(j+2)%3];
                invM(i,j) -= M[(i+2)%3][(j+1)%3]*M[(i+1)%3][(j+2)%3];
                invM(i,j) /= det;
            }
        }
        return invM;
    }

    // ---------------------------------------------------|
    // Computes the L2 vector norm of a general vector    |
    // ---------------------------------------------------|

    double norm(darray &V){

        double norm = 0;
        for(int i=0; i<V.length(); i++){norm += V[i]*V[i];}
        norm = sqrt(norm);
        return norm;
    }

    // ------------------------------------------------------|
    // Computes the L2 vector norm of a general 3D vector    |
    // ------------------------------------------------------|

    double norm(array3d &V){

        double norm = 0;
        for(int i=0; i<3; i++){norm += V[i]*V[i];}
        norm = sqrt(norm);
        return norm;
    }

    // -----------------------------------------------------------------------|
    // Cleans and stores the non-zero indices of a sparse symmetric matrix    |
    // -----------------------------------------------------------------------|

    vector<ivector> sparsemap(sparse &M){

        double val;
        alglib::ae_int_t i=0,j=0;
        alglib::ae_int_t I=0,J=0;
        int nLen = alglib::sparsegetnrows(M);
        vector<ivector> row(3*nLen);
        
        // Stores the non-zero indice locations per row and column

        while(alglib::sparseenumerate(M,I,J,i,j,val)){

            row[i].push_back(j);
            if(i!=j){row[j].push_back(i);}
        }
        return row;
    }

    // ---------------------------------------------------------------|
    // Sets the coordinate (row,col) for a symmetric sparse matrix    |
    // ---------------------------------------------------------------|

    void symset(sparse &M,int row,int col,double val){

        if(row>col){alglib::sparseset(M,col,row,val);}
        else{alglib::sparseset(M,row,col,val);}
    }

    // ------------------------------------------------------------------|
    // Adds to the coordinate (row,col) for a symmetric sparse matrix    |
    // ------------------------------------------------------------------|

    void symadd(sparse &M,int row,int col,double val){

        if(row>col){alglib::sparseadd(M,col,row,val);}
        else{alglib::sparseadd(M,row,col,val);}
    }

    // ------------------------------------------------------------------------|
    // gets the value at coordinate (row,col) for a symmetric sparse matrix    |
    // ------------------------------------------------------------------------|

    double get(sparse &M,int row,int col){

        double val;
        if(row>col){val = alglib::sparseget(M,col,row);}
        else{val = alglib::sparseget(M,row,col);}
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

        for(int n=0; n<dim; n++){
            dvector N1;

            // Computes the derivative of Lagrange polynomial for this variable

            if(n==dvar){
                N1.resize(sLen,0);

                for(int j=0; j<sLen; j++){
                    for(int i=0; i<sLen; i++){

                        // Computes the first term of the summ of products

                        if(i!=j){
                            double La = 1/(node[j]-node[i]);

                            for(int k=0; k<sLen; k++){
                                if(k!=j && k!=i){La *= (val[n]-node[k])/(node[j]-node[k]);}
                            }
                            N1[j] += La;
                        }
                    }
                }
            }
   
            // Computes the Lagrange polynomial for other variables

            else{
                N1.resize(sLen,1);

                for(int i=0; i<sLen; i++){
                    for(int j=0; j<sLen; j++){
                        if(i!=j){N1[j] *= (val[n]-node[i])/(node[j]-node[i]);}
                    }
                }
            }

            // Product of the Lagrange polynomials in the 2D case

            if(dim==2){
                for(int i=0; i<sLen; i++){
                    for(int j=0; j<sLen; j++){
                            
                        ivector loop = {i,j};
                        N[i*sLen+j] *= N1[loop[n]];
                    }
                }
            }

            // Product of the Lagrange polynomials in the 3D case

            if(dim==3){
                for(int i=0; i<sLen; i++){
                    for(int j=0; j<sLen; j++){
                        for(int k=0; k<sLen; k++){
                            
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

    matrix projection(array3d norm){

        matrix T,P;
        P.setlength(3,3);
        T.setlength(6,6);

        // Computes the surface gredient ∇s = (I-n⊗n)∇

        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                
                if(j==k){P(j,k) = 1-norm[j]*norm[k];}
                else{P(j,k) = -norm[j]*norm[k];}
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