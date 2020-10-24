#include "..\include\func.h"
using namespace std;

// Quadrature rule generator (to be modified)

namespace math{

    quadStruct quadrature(int dim){

        quadStruct quad;
        vector<double> weight;
        vector<darray> gRST;
        vector<double> w = {5.0/9,5.0/9,8.0/9};
        vector<double> g = {-sqrt(3.0/5),sqrt(3.0/5),0};

        if(dim==3){

            // Cube element quadrature rule
                
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    for(int k=0; k<3; k++){
                        
                        weight.push_back(w[i]*w[j]*w[k]);
                        double v1[] = {g[i],g[j],g[k]};
                        darray v2; v2.setcontent(3,v1);
                        gRST.push_back(v2);
                    }
                }
            }
        }

        if(dim==2){

            // Square face quadrature rule

            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                        
                    double v1[] = {g[i],g[j]};
                    weight.push_back(w[i]*w[j]);
                    darray v2; v2.setcontent(2,v1);
                    gRST.push_back(v2);
                }
            }
        }

        quad.weight = weight;
        quad.gRST = gRST;
        return quad;
    }

    // Converts a 3D square into 2D square

    vector<darray> to2D(vector<darray> &nXYZ){

        vector<darray> nXY;
        vector<double> n{0,0,0};
        darray v1; v1.setlength(3);
        darray v2; v2.setlength(3);
        darray v3; v3.setlength(3);
        double n21 = 0;
        double n31 = 0;

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
        }

        double d[2] = {sqrt(n21/n[0]),sqrt(n31/n[0])};
        double l[2] = {sqrt(n[1]-n21/n[0]),sqrt(n[2]-n31/n[0])};
        double xy[4][2] {{0,0},{sqrt(n[0]),0},{l[0],d[0]},{l[1],d[1]}};

        for(int i=0; i<4; i++){

            darray arr; arr.attach_to_ptr(2,xy[i]);
            nXY.push_back(arr);
        }
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

        for(int i=0; i<V.length(); i++){
            V(i) = 0;
        }
    }

    // Removes quasi-zero elements from a sparse matrix

    void clean(sparse &M){
    
        double val;
        alglib::ae_int_t i,j;
        alglib::ae_int_t end;
        alglib::ae_int_t start;

        while(alglib::sparseenumerate(M,start,end,i,j,val)){
            if(val*val<1e-25){alglib::sparserewriteexisting(M,i,j,0);}
        }
    }
}