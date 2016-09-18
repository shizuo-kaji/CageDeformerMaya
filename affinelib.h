/**
 * @file affinelib.h
 * @brief Library for 3D affine transformation
 * @brief An implementation of the paper "A concise parametrisation of affine transformation" by S. Kaji and H. Ochiai
 * @section LICENSE
 *                   the MIT License
 * @section Requirements
 *                   Eigen library,   (optional) MKL
 * @section CAUTION
 *           The convention used here is different from the one in the paper;
 *           we assume that matrices act on row vectors by right multiplication.
 *           ( that is, everything is transposed compared to the paper. )
 * @version 0.30
 * @date  Jun. 2016
 * @author Shizuo KAJI
 */


#pragma once

// uncomment if you use MKL
// #define  EIGEN_USE_MKL_ALL

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/StdVector>
#include <Eigen/Geometry>

#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>


/// threshold for small values to be regarded zero
#define EPSILON 10e-15
// for assert()
#define TOLERANCE 10e-5

/// 3x3 identity matrix
#define Id3 Matrix3d::Identity()
/// macro to print an Eigen object
#define PRINT_MAT(X) std::cout << #X << ":\n" << X << std::endl << std::endl

//
using namespace Eigen;
using namespace std;
/// For vecterization of Eigen objects
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Matrix4d);
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Vector4d);

// main body
namespace AffineLib{
    Matrix4d pad(const Matrix3d& m, const Vector3d& l, const double br=1.0)
    /** compose an affine matrix from linear matrix and translation vector
     * @param m 3x3-matrix
     * @param l 3-dim translation vector
     * @param br value of bottom right corner; set to 0.0 for log matrix
     * @return 4x4-affine transformation matrix
     */
    {
        Matrix4d aff;
        aff << m(0,0),m(0,1),m(0,2),0.0,
        m(1,0),m(1,1),m(1,2),0.0,
        m(2,0),m(2,1),m(2,2),0.0,
        l[0],l[1],l[2],br;
        return aff;
    }
    
    RowVector4d pad(const Vector3d& v)
    /** compose a 4-vector by concatenating 1 at the 4th coordinate
     * @param v 3-vector
     * @return 4-vector
     */
    {
        return RowVector4d(v[0],v[1],v[2],1);
    }

    Vector3d transPart(const Matrix4d& m){
        /** extract the translation vector from an affine matrix
         * @param m affine matrix
         * @return 3D translation vector
         */
        return Vector3d(m(3,0),m(3,1),m(3,2));
    }
    
    Matrix3d logSO(const Matrix3d& m)
    /** Log for a rotational matrix using Axis-Angle
     * @param m rotational matrix
     * @return primary value of log(m)
     */
    {
        assert( (m*m.transpose()-Id3).squaredNorm() < TOLERANCE );
        AngleAxisd X(m);
        Matrix3d A;
        A << 0,     -X.axis()[2], X.axis()[1],
        X.axis()[2], 0,        -X.axis()[0],
        -X.axis()[1], X.axis()[0], 0;
        return X.angle() * A;
    }

    Matrix3d logSOc(const Matrix3d& m, const Matrix3d& P)
    /** "Continuous" log for a rotational matrix
     * @param m rotational matrix
     * @param P anti-symmetric matrix
     * @return the branch of log(m) closest to P
     */
    {
        assert( (m*m.transpose()-Id3).squaredNorm() < TOLERANCE );
        AngleAxisd X(m);
        Matrix3d A;
        double theta=X.angle();
        Vector3d n=X.axis();
        Vector3d prevN(-P(1,2),P(0,2),-P(0,1));
        double prevTheta=prevN.norm();
        if(abs(sin(theta))==0.0){
            n=prevN;
        }
        A << 0,     -n[2], n[1],
        n[2], 0,        -n[0],
        -n[1], n[0], 0;
        if(n.dot(prevN)<0){
            A = -A;
            theta = -theta;
        }
        while(theta-prevTheta>M_PI){
            theta -= 2*M_PI;
        }
        while(prevTheta-theta>M_PI){
            theta += 2*M_PI;
        }
        return(theta*A);
    }
    
    Matrix3d expTaylor(const Matrix3d& m, const int deg = 50)
        /** exp by Taylor expansion
         * @param m 3x3 matrix
         * @param deg degree of Taylor expansion
         * @return exp(m)
         */
    {
        Matrix3d A = Id3;
        Matrix3d mPow = Id3;
        for(int i=1;i<deg+1;i++){
            mPow = m*mPow/i;
            A += mPow;
        }
        return A;
    }

    Matrix3d logTaylor(const Matrix3d& m, const int deg=50)
    /** log by Taylor expansion
     * @param m 3x3 matrix
     * @param deg degree of Taylor expansion
     * @return log(m)
     */
    {
        Matrix3d A = Matrix3d::Zero();
        Matrix3d mPow = -Id3;
        for(int i=1;i<deg+1;i++){
            mPow = (Id3-m)*mPow;
            A += mPow/i;
        }
        return A;
    }

    
    Matrix3d expSO(const Matrix3d& m)
    /** exp for an anti-symmetric matrix using Rodrigues' formula
     * @param m anti-symmetric matrix
     * @return exp(m)
     */
    {
        assert( ((m + m.transpose())).squaredNorm() < TOLERANCE );
        double norm2=m(0,1)*m(0,1) + m(0,2)*m(0,2) + m(1,2)*m(1,2);
        if(norm2<EPSILON){
            return Id3 + m + m*m/2.0;
        }else{
            double norm = sqrt(norm2);
            return Id3 + sin(norm)/norm * m + (1.0-cos(norm))/norm2 * m*m;
        }
    }
    
    Matrix4d expSE(const Matrix4d& mm)
    /** exp for the log of a rigid transformation (screw) matrix
     * @param mm the log of a rigid transformation matrix
     * @return exp(mm)
     */
    {
        Matrix3d m = mm.block(0,0,3,3);
        assert( ((m + m.transpose())).squaredNorm() < TOLERANCE );
        Vector3d v;
        v << mm(3,0), mm(3,1), mm(3,2);
        double norm2=m(0,1)*m(0,1) + m(0,2)*m(0,2) + m(1,2)*m(1,2);
        Matrix3d A,ans;
        if(norm2<EPSILON){
            return (Matrix4d::Identity() + mm + mm*mm/2.0);
        }else{
            double norm = sqrt(norm2);
            ans = Id3 + sin(norm)/norm * m + (1.0-cos(norm))/norm2 * m*m;
            A = Id3 + (1.0-cos(norm))/norm2 * m + (norm-sin(norm))/(norm*norm2) * m*m ;
            return pad(ans, A.transpose()*v);
        }
    }
    
    Matrix4d logSEc(const Matrix4d& mm, const Matrix4d& P = Matrix4d::Zero())
    /** "Continuous" log for a rigid transformation (screw) matrix
     * @param mm rigid transformation matrix
     * @param P log matrix
     * @return branch of log(mm) closest to P
     */
    {
        Matrix3d m = mm.block(0,0,3,3);
        assert( ((m * m.transpose()) - Id3).squaredNorm() < TOLERANCE );
        Vector3d v = transPart(mm);
        Matrix3d X = logSOc(m, P.block(0,0,3,3));
        double theta=sqrt(X(1,2)*X(1,2) + X(0,2)*X(0,2) + X(0,1)*X(0,1));
        Matrix3d A;
        Vector3d l;
        if(abs(1-cos(theta))<EPSILON){
            double prevTheta = sqrt(P(1,2)*P(1,2) + P(0,2)*P(0,2) + P(0,1)*P(0,1));
            if( prevTheta < EPSILON ){
                l = Vector3d::Zero();
            }else{
                l = theta / prevTheta * transPart(P);
            }
        }else if(abs(1+cos(theta))<EPSILON){
            A = Id3 - 0.5 * X + 1.0/(theta*theta) * X*X;
            l = A.transpose()*v;
        }else{
            A = Id3 - 0.5 * X + (1.0/(theta*theta) - 0.5*(1.0+cos(theta))/(sin(theta)*theta)) * X * X;
            l = A.transpose()*v;
        }
        return(pad(X, l,0.0));
    }
    
    
    Matrix3d logDiag(const Matrix3d& U, const Vector3d& s)
    /** log for a diagonalised matrix  m = U diag(s) U^t with positive eigenvalues
     * @param U rotational diagonalising matrix
     * @param s eigenvalues (must be all positive)
     * @return log(m)
     */
    {
        assert( s[0] > 0 && s[1] > 0 && s[2] > 0);
        Vector3d d = s.array().log();
        return(U* d.asDiagonal() *U.transpose());
    }
    
    Matrix3d expDiag(const Matrix3d& U, const Vector3d& s)
    /** exp for a diagonalised matrix  m = U diag(s) U^{-1}
     * @param U diagonalising matrix
     * @param s eigenvalues
     * @return exp(m)
     */
    {
        Vector3d d = s.array().exp();
        return(U* d.asDiagonal() *U.inverse());
    }
    
    Matrix3d expSymDiag(const Matrix3d& m)
    /** exp for a symmetric matrix by diagonalization
     * @param m symmetric matrix
     * @return exp(m)
     */
    {
        assert( ((m - m.transpose())).squaredNorm() < TOLERANCE );
        if (m.squaredNorm() < EPSILON){
            return Id3+m+0.5*m*m;
        }
        SelfAdjointEigenSolver<Matrix3d> eigensolver;
        eigensolver.compute(m);
        Vector3d s(eigensolver.eigenvalues());
        Matrix3d U(eigensolver.eigenvectors());
        s = s.array().exp();
        return(U * s.asDiagonal() * U.transpose());
    }
    
    Matrix3d logSymDiag(const Matrix3d& m)
    /** log for a symmetric matrix by diagonalization
     * @param m symmetric matrix
     * @return log(m)
     */
    {
        assert( ((m - m.transpose())).squaredNorm() < TOLERANCE );
        if ((m-Id3).squaredNorm() < EPSILON){
            return m-Id3-0.5*(m-Id3)*(m-Id3);
        }
        SelfAdjointEigenSolver<Matrix3d> eigensolver;
        eigensolver.compute(m);
        Vector3d s(eigensolver.eigenvalues());
        Matrix3d U(eigensolver.eigenvectors());
        s = s.array().log();
        return(U * s.asDiagonal() * U.transpose());
    }
    
    Matrix3d expSym(const Matrix3d& m, Vector3d e=Vector3d::Zero())
    /** exp for a symmetric matrix by spectral decomposition
     * @param m symmetric matrix
     * @param e if eigenvalues of m are provided, we use it
     * @return exp(m)
     */
    {
        assert( ((m - m.transpose())).squaredNorm() < TOLERANCE );
        if(e == Vector3d::Zero()){
            // compute eigenvalues if not given
            // eigenvalues are sorted in increasing order.
            SelfAdjointEigenSolver<Matrix3d> eigensolver;
            eigensolver.computeDirect(m, EigenvaluesOnly);
            e = eigensolver.eigenvalues();
            if(abs(e[0])<TOLERANCE){
                eigensolver.compute(m, EigenvaluesOnly);
                e = eigensolver.eigenvalues();
            }
        }
        Matrix3d A = m-e[1]*Id3;
        if (A.squaredNorm() < EPSILON){
            return exp(e[1])*(Id3+A+0.5*A*A);
        }
        double x(e[0]-e[1]),y(e[2]-e[1]);
        double t2ex = abs(x)>EPSILON ? (exp(x)-1-x)/(x*x) : 0.5+x/6+x*x/24;
        double t2ey = abs(y)>EPSILON ? (exp(y)-1-y)/(y*y) : 0.5+y/6+y*y/24;
        
        double b = 1- x*y*(t2ex-t2ey)/(x-y);
        double c = (x*t2ex- y*t2ey)/(x-y);
        Matrix3d ans(exp(e[1])*( Id3 + b*A + c*A*A));
        return (ans+ans.transpose())/2;
    }
    
    
    Matrix3d logSym(const Matrix3d& m, Vector3d& lambda)
    /** log for a positive definite symmetric matrix by spectral decomposition
     * @param m symmetric matrix
     * @param lambda returns eigen values for log(m)
     * @return log(m)
     */
    {
        assert( ((m - m.transpose())).squaredNorm() < TOLERANCE );
        // compute eigenvalues only
        // eigenvalues are sorted in the increasing order.
        SelfAdjointEigenSolver<Matrix3d> eigensolver;
        eigensolver.computeDirect(m, EigenvaluesOnly);
        Vector3d e(eigensolver.eigenvalues());
        if(abs(e[0])<TOLERANCE){
            eigensolver.compute(m, EigenvaluesOnly);
            e = eigensolver.eigenvalues();
        }
        assert(e[0] > 0 && e[1] > 0 && e[2] > 0);
        lambda = e.array().log();
        double x(e[0]/e[1]),y(e[2]/e[1]);
        double t2lx = abs(x-1)>EPSILON ? (log(x)-x+1)/(x-1) : -(x-1)/2+(x-1)*(x-1)/3;
        double t2ly = abs(y-1)>EPSILON ? (log(y)-y+1)/(y-1) : -(y-1)/2+(y-1)*(y-1)/3;
        double a,c;
        if(abs(x-y)>EPSILON){
            a = -1 + (y*t2lx - x*t2ly)/(x-y);
            c = (t2lx - t2ly)/(x-y);
        }else{
            a = -1/6 + x*y/3;
            c = -0.5 + (x+y)/3;
        }
        Matrix3d ans((a+log(e[1]))*Id3 - (a+c)/e[1]*m + c/(e[1]*e[1])*m*m);
        return (ans+ans.transpose())/2;
    }
    
    
    Matrix3d frechetSO(const std::vector<Matrix3d> &m, const std::vector<double> &w, const int max_step=10)
    /** the Frechet mean for rotations
     * @param m array of rotation matrices to be averaged
     * @param w array of weights
     * @param max_step max steps for iteration
     * @return weighted Frechet mean
     */
    {
        assert(m.size() == w.size());
        if(m.empty()) return(Id3);
        Matrix3d Z = m[0];
        for(int i=0;i<max_step;i++){
            Matrix3d W = Matrix3d::Zero();
            Matrix3d ZI = Z.transpose();
            for(int j=0;j<m.size();j++){
                W += w[j] * logSO(ZI * m[j]);
            }
            W = (W-W.transpose())/2;
            if(W.squaredNorm()<EPSILON) break;
            Z = Z * expSO(W);
        }
        return(Z);
    }
    
    
    Matrix3d frechetSym(const std::vector<Matrix3d> &m, const std::vector<double> &w, const int max_step=10)
    /** the Frechet mean for symmetric matrices
     * @param m array of positive definite symmetric matrices to be averaged
     * @param w array of weights
     * @param max_step max steps for iteration
     * @return weighted Frechet mean
     */
    {
        assert(m.size() == w.size());
        if(m.empty()) return(Id3);
        Vector3d e;
        Matrix3d Z = m[0];
        for(int i=0;i<max_step;i++){
            Matrix3d W = Matrix3d::Zero();
            Matrix3d ZI = Z.inverse();
            for(int j=0;j<m.size();j++){
                W += w[j] * logSym(ZI * m[j], e);
            }
            W = (W+W.transpose())/2;
            if(W.squaredNorm()<EPSILON) break;
            Z = Z * expSO(W);
        }
        return(Z);
    }
    
    void polarDiag(const Matrix3d& m, Matrix3d& U, Vector3d& s, Matrix3d& R)
    /** Polar decomposition m = U diag(s) U^T R by diagonalisation
     * @param m matrix to be decomposed
     * @param U diagonaliser of symmetric part
     * @param s singular values
     * @param R rotation part
     */
    {
        assert(m.determinant()>0);
        Matrix3d A= m*m.transpose();
        SelfAdjointEigenSolver<Matrix3d> eigensolver;
        eigensolver.computeDirect(A);
        s = eigensolver.eigenvalues();
        U = Matrix3d(eigensolver.eigenvectors());
        s = s.array().sqrt();
        Vector3d si = s.array().inverse();
        R = U * si.asDiagonal() * U.transpose() * m;
    }
    
    void polarN(const MatrixXf& m, MatrixXf& S, MatrixXf& R)
    /** Polar decomposition m = S R for square matrix of any size
     * @param m matrix to be decomposed
     * @param S symmetric part
     * @param R rotation part
     */
    {
        assert(m.determinant()>0);
        MatrixXf A= m*m.transpose();
        long N = A.rows();
        SelfAdjointEigenSolver<MatrixXf> eigensolver;
        eigensolver.compute(A, EigenvaluesOnly);
        VectorXf s = eigensolver.eigenvalues();
        MatrixXf VM = MatrixXf::Zero(N,N);
        double ss;
        for(int i=0;i<N; i++){
            ss = 1.0;
            for(int j=0;j<N; j++){
                VM(i,j) = ss;
                ss *= s[i];
            }
        }
        // compute logS
        VectorXf logs = s.array().log();
        VectorXf a = VM.colPivHouseholderQr().solve(logs);
        MatrixXf logS = MatrixXf::Zero(N,N);
        MatrixXf AA = MatrixXf::Identity(N,N);
        for(int i=0;i<N; i++){
            logS += a[i] * AA;
            AA *= A;
        }
        logS /= 2.0;
        // compute S
        s = logs/2.0;
        for(int i=0;i<N; i++){
            ss = 1.0;
            for(int j=0;j<N; j++){
                VM(i,j) = ss;
                ss *= s[i];
            }
        }
        VectorXf exps = s.array().exp();
        VectorXf b = VM.colPivHouseholderQr().solve(exps);
        S = MatrixXf::Zero(N,N);
        AA = MatrixXf::Identity(N,N);
        for(int i=0;i<N; i++){
            S += b[i] * AA;
            AA *= logS;
        }
        // compute R
        s = -s;
        for(int i=0;i<N; i++){
            ss = 1.0;
            for(int j=0;j<N; j++){
                VM(i,j) = ss;
                ss *= s[i];
            }
        }
        VectorXf expsinv = exps.array().inverse();
        VectorXf c = VM.colPivHouseholderQr().solve(expsinv);
        MatrixXf SINV = MatrixXf::Zero(N,N);
        AA = MatrixXf::Identity(N,N);
        for(int i=0;i<N; i++){
            SINV += c[i] * AA;
            AA *= -logS;
        }
        R = SINV * m;
    }
    
    
    void polarBySVD(const Matrix3d& m, Matrix3d& U, Vector3d& s, Matrix3d& R){
        /** Polar decomposition m = U diag(s) U^T R  by SVD
         * @param m matrix to be decomposed
         * @param U diagonaliser of symmetric part
         * @param s singular values
         * @param R rotation part
         */
        JacobiSVD<Matrix3d> svd(m, ComputeFullU | ComputeFullV);
        U = svd.matrixU();
        s = svd.singularValues();
        R = svd.matrixU() * svd.matrixV().transpose();
    }
    
    void polarByParam(const Matrix3d& m, Matrix3d& S, Matrix3d& R){
        /** Polar decomposition m = S R   by parametrisation map
         * @param m matrix to be decomposed
         * @param S shear part
         * @param R rotation part
         */
        Vector3d lambda=Vector3d::Zero().eval();
        Matrix3d logS = logSym(m*m.transpose(), lambda)/2.0;
        S = expSym(logS, lambda/2.0);
        R = expSym(-logS, -lambda/2.0) * m;
    }
    
    int polarHigham(const Matrix3d& m, Matrix3d& S, Matrix3d& R){
        /** Polar decomposition m = S R   by Higham's iterative method
         * @param m matrix to be decomposed
         * @param S shear part
         * @param R rotation part
         * @return number of iterations
         */
        Matrix3d Curr = m;
        Matrix3d Prev;
        int iter=0;
        do {
            assert(Curr.determinant() != 0.0);
            Matrix3d Ad = Curr.inverse().transpose();
            double nad = Ad.array().abs().colwise().sum().maxCoeff() * Ad.array().abs().rowwise().sum().maxCoeff();
            double na = Curr.array().abs().colwise().sum().maxCoeff() * Curr.array().abs().rowwise().sum().maxCoeff();
            double gamma = sqrt(sqrt(nad / na));
            //        std::cout << gamma << std::endl;
            Prev = Curr;
            Curr = (0.5*gamma)*Curr + (0.5/gamma) *Ad;
            iter++;
        } while ((Prev-Curr).lpNorm<1>() > EPSILON*Prev.lpNorm<1>());
        R = Curr;
        S = m * Curr.transpose();
        return iter;
    }
    

    void parametriseGL(const Matrix3d& m, Matrix3d& logS, Matrix3d& R)
    /** Parametrisation map for GL(3)
     * @param m linear matrix to be decomposed
     * @param logS log of shear part
     * @param R rotation part
     */
    {
        assert(m.determinant()>0);
        Vector3d lambda=Vector3d::Zero();
        logS = logSym(m*m.transpose(), lambda)/2.0;
        R = expSym(-logS, -lambda/2.0) * m;
    }
    
    template<typename T>
    T blendMat(const std::vector<T>& A, const std::vector<double>& weight){
    /** blend matrices
     * @param A array of matrices
     * @param weight array of weights
     * @return blended matrix
     */
        assert(A.size() == weight.size());
        T X=T::Zero();
        for(int i=0;i<A.size();i++){
            X += weight[i]*A[i];
        }
        return X;
    }
    
    template<typename T>
    T blendMatLin(const std::vector<T>& A, const std::vector<double>& weight){
        /** blend matrices; if weight doesn't sum up to one, the result will be complimented by identity
         * this is suitable for linear blending
         * @param A list of matrices
         * @param weight list of weights
         * @return blended matrix
         */
        assert(A.size() == weight.size());
        T X=T::Zero();
        double sum = 0.0;
        for(int i=0;i<A.size();i++){
            X += weight[i]*A[i];
            sum += weight[i];
        }
        X += (1.0-sum) * T::Identity();
        return X;
    }
    
    Vector4d blendQuat(const std::vector<Vector4d>& A, const std::vector<double>& weight){
        /** blend 4-vector (quaternion); if weight doesn't sum up to one, the result will be complimented by 1
         * this is suitable for linear blending of quaternions
         * @param A list of 4-vectors
         * @param weight list of weights
         * @return blended 4-vector
         */
        assert(A.size() == weight.size());
        Vector4d I(0,0,0,1);
        Vector4d X=Vector4d::Zero();
        double sum = 0.0;
        for(int i=0;i<A.size();i++){
            X += weight[i] * A[i];
            sum += weight[i];
        }
        X += (1.0-sum) * I;
        return X.normalized();
    }
}

