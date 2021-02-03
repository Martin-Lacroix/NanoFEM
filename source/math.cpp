#include "..\include\math.h"
using namespace std;

namespace math{

    // Non-local kernel density function

    double kernel(dvector x,dvector y){

        double d = 0.05;
        double norm = 1/(8*d*d*d);
        double dist = abs(x[0]-y[0])+abs(x[1]-y[1])+abs(x[2]-y[2]);
        double k = norm*exp(-dist/d);
        return k;
    }

    // Generates the stiffness tensor

    matrix stiffness(double E,double v){

        matrix D;
        D.setlength(6,6);
        math::zero(D);

        double mu = E/(2*(1+v));
        double lam = E*v/((1+v)*(1-2*v));
        D(0,1) = D(0,2) = D(1,2) = lam;
        D(1,0) = D(2,0) = D(2,1) = lam;

        for(int i=0; i<3; i++){

            D(i,i) = 2*mu+lam;
            D(i+3,i+3) = mu;
        }
        return D;
    }

    // Legendre quadrature rule generator

    quadStruct legendre(int dim,int order){

        matrix M;
        quadStruct quad;
        int nbr = order+1;
        darray Ja; Ja.setlength(nbr);
        darray Jb; Jb.setlength(nbr-1);

        // Generates the three term coefficients
        
        for(int i=1; i<nbr; i++){Jb(i-1) = i/sqrt(4*i*i-1);}
        for(int i=0; i<nbr; i++){Ja(i) = 0;}
        alglib::smatrixtdevd(Ja,Jb,nbr,3,M);

        // Cube element quadrature rule

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

        // Square face quadrature rule

        if(dim==2){
            for(int i=0; i<nbr; i++){
                for(int j=0; j<nbr; j++){

                    quad.gRST.push_back({Ja(i),Ja(j)});
                    quad.weight.push_back(4*M(0,i)*M(0,i)*M(0,j)*M(0,j));
                }
            }
        }
        return quad;
    }

    // Converts a 3D square into 2D square

    vector<dvector> to2D(vector<dvector> &nXYZ){

        dvector n{0,0,0};
        darray v1; v1.setlength(3);
        darray v2; v2.setlength(3);
        darray v3; v3.setlength(3);
        double n21=0, dot2=0;
        double n31=0, dot3=0;

        for(int i=0; i<3; i++){

            v1[i] = nXYZ[1][i]-nXYZ[0][i];
            v2[i] = nXYZ[2][i]-nXYZ[0][i];
            v3[i] = nXYZ[3][i]-nXYZ[0][i];
            
            n[0] += v1[i]*v1[i];
            n[1] += v2[i]*v2[i];
            n[2] += v3[i]*v3[i];
        }

        darray v21 = cross(v2,v1);
        darray v31 = cross(v3,v1);

        for(int i=0; i<3; i++){
            
            n21 += v21[i]*v21[i];
            n31 += v31[i]*v31[i];
            dot2 += v2[i]*v1[i];
            dot3 += v3[i]*v1[i];
        }

        double l2 = copysign(sqrt(abs(n[1]-n21/n[0])),dot2);
        double l3 = copysign(sqrt(abs(n[2]-n31/n[0])),dot3);
        vector<dvector> nXY = {{0,0},{sqrt(n[0]),0},{l2,sqrt(n21/n[0])},{l3,sqrt(n31/n[0])}};
        return nXY;
    }

    // Standard vector cross product V3 = V1 × V2

    darray cross(darray &V1,darray &V2){

        darray V3; V3.setlength(3);
        V3[0] = V1[1]*V2[2]-V1[2]*V2[1];
        V3[1] = V1[0]*V2[2]-V1[2]*V2[0];
        V3[2] = V1[0]*V2[1]-V1[1]*V2[0];
        return V3;
    }

    // Standard matrix-matrix addition M2 = k1 M1 + k2 M2

    void add(double k1,double k2,matrix &M1,matrix &M2){

        int m = M1.rows();
        int n = M1.cols();
        alglib::rmatrixgencopy(m,n,k1,M1,0,0,k2,M2,0,0);
    }

    // Matrix-matrix product M3 = k op(M1) op(M2), tj = matrix Mj is transposed

    matrix prod(double k,matrix &M1,matrix &M2,int t1,int t2){

        int m = M1.rows();
        int n = M2.cols();
        int w = M1.cols();

        if(t1==1){w = M1.rows(); m = M1.cols();}
        if(t2==1){n = M2.rows();}

        matrix M3; M3.setlength(m,n);
        alglib::rmatrixgemm(m,n,w,k,M1,0,0,t1,M2,0,0,t2,0,M3,0,0);
        return M3;
    }

    // Matrix-vector product V2 = k op(M1) V1, t1 = matrix M1 is transposed

    darray prod(double k,matrix &M1,darray &V1,int t1){

        int m = M1.rows();
        int n = M1.cols();

        if(t1==1){
            n = M1.rows();
            m = M1.cols();
        }

        darray V2; V2.setlength(m);
        alglib::rmatrixgemv(m,n,k,M1,0,0,t1,V1,0,0,V2,0);
        return V2;
    }

    // Fills a dense matrix or an array with zeros

    void zero(matrix &M){

        for(int i=0; i<M.rows(); i++){
            for(int j=0; j<M.cols(); j++){
                M(i,j) = 0;
            }
        }
    }

    void zero(darray &V){
        for(int i=0; i<V.length(); i++){V(i) = 0;}
    }

    // Sparse matrix addition M2 = k1 M1 + k2 M2

    void add(double k1,double k2,sparse &M1,sparse &M2){
    
        double val;
        alglib::ae_int_t i=0,j=0;
        alglib::ae_int_t I=0,J=0;

        while(alglib::sparseenumerate(M2,I,J,i,j,val)){
            alglib::sparserewriteexisting(M2,i,j,k2*val);
        }
        while(alglib::sparseenumerate(M1,I,J,i,j,val)){
            alglib::sparseadd(M2,i,j,k1*val);
        }
    }

    // Cleans and stores the non-zero indices of a sparse matrix

    matStruct mapclean(sparse &K){

        double val;
        matStruct mat;
        alglib::ae_int_t i=0,j=0;
        alglib::ae_int_t I=0,J=0;

        int nLen = alglib::sparsegetnrows(K);
        mat.row.resize(3*nLen);
        mat.col.resize(3*nLen);

        // Stores the non-zero indices per row and column

        while(alglib::sparseenumerate(K,I,J,i,j,val)){

            if(abs(val)<1e-16){
                alglib::sparseset(K,i,j,0);
            }
            else{
                mat.row[i].push_back(j);
                mat.col[j].push_back(i);
            }
        }
        return mat;
    }
}