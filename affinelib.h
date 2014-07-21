/**
 * @file affinelib.h
 * @brief Library for 3D affine transformation
 * @section LICENSE The MIT License
 * @section  requirements:  Eigen library,   (optional) MKL
 * @section  CAUTION: the convention here is different from the one in the paper;
 *           we assume that matrices act on row vectors by right multiplication.
 *           ( that is, everything is transposed compared to the paper )
 * @version 0.20
 * @date  Nov. 2013
 * @author Shizuo KAJI
 */


#pragma once

// uncomment if you use MKL
// #define  EIGEN_USE_MKL_ALL

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/StdVector>
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>


// threshold for being zero
#define EPSILON 0.000001
// Identity matrix
#define E Matrix3d::Identity()
// for debug print
#define PRINT_MAT(X) std::cout << #X << ":\n" << X << std::endl << std::endl

using namespace Eigen;
using namespace std;
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Matrix4d);

namespace AffineLib{
    Matrix4d affine(const Matrix3d& m, const Vector3d& l)
    /** compose affine matrix from linear matrix and translation vector
     * @param m 3x3-matrix
     * @param l 3-dim translation vector
     * @return 4x4-affine transformation matrix
     */
    {
        Matrix4d aff;
        aff << m(0,0),m(0,1),m(0,2),0.0,
        m(1,0),m(1,1),m(1,2),0.0,
        m(2,0),m(2,1),m(2,2),0.0,
        l[0],l[1],l[2],1.0;
        return aff;
    }
    
    Matrix4d affine0(const Matrix3d& m, const Vector3d l)
    /** compose 4x4-matrix in log space
     * @param m 3x3-matrix
     * @param l 3-dim vector
     * @ return 4x4 matrix in log space
     */
    {
        Matrix4d aff;
        aff << m(0,0),m(0,1),m(0,2),0.0,
        m(1,0),m(1,1),m(1,2),0.0,
        m(2,0),m(2,1),m(2,2),0.0,
        l[0],l[1],l[2],0.0;
        return aff;
    }
    
    Vector3d transPart(const Matrix4d& m){
        return Vector3d(m(3,0),m(3,1),m(3,2));
    }
    
    Matrix3d logSO(const Matrix3d& m)
    /** Log for rotational matrix using Rodrigues' formula
     * @param m rotational matrix
     * @return primary value of log(m)
     */
    {
        assert( ((m * m.transpose()) - E).squaredNorm() < EPSILON );
        double tr=(m(0,0)+m(1,1)+m(2,2)-1.0f)/2.0;
        double theta;
        Matrix3d ans=Matrix3d::Zero();
        if(tr>=1.0){
            return ans;
        }else if(tr<=-1.0){
            ans(0,1)=M_PI;
            return ans;
        }else{
            theta=acosf(tr);
        }
        if(theta<EPSILON){
            return ans;
        }else if(M_PI-theta<EPSILON){
            ans(0,1)=M_PI;
            return ans;
        }else{
            ans=0.5f*theta/sinf(theta) * (m-m.transpose());
            return ans;
        }
    }
    
    Matrix3d logSOc(const Matrix3d& m, double& prevTheta, Vector3d& prevN)
    /** "Continuous" log for rotational matrix
     * @param m rotational matrix
     * @param prevTheta rotational angle
     * @param prevN rotational axis vector
     * @return branch of log(m) with axis and angle closest to the given ones
     */
    {
        //    assert( ((m * m.transpose()) - E).squaredNorm() < EPSILON );
        double tr=(m(0,0)+m(1,1)+m(2,2)-1.0)/2.0;
        double theta;
        Matrix3d ans=Matrix3d::Zero();
        if(tr>=1.0){
            theta=0.0;
        }else if(tr<=-1.0){
            theta=M_PI;
        }else{
            theta=acosf(tr);
        }
        if(theta<EPSILON || M_PI-theta<EPSILON){
            ans << 0,prevN[0],prevN[1],-prevN[0],0,prevN[2],-prevN[1],-prevN[2],0;
        }else{
            ans=0.5f/sinf(theta) * (m-m.transpose());
        }
        if(ans(0,1)*prevN[0] + ans(0,2)*prevN[1] + ans(1,2)*prevN[2]<0){
            ans = -ans;
            theta = -theta;
        }
        while(theta-prevTheta>M_PI){
            theta -= 2*M_PI;
        }
        while(prevTheta-theta>M_PI){
            theta += 2*M_PI;
        }
        prevTheta=theta;
        prevN << ans(0,1), ans(0,2), ans(1,2);
        return(theta*ans);
    }
    
    Matrix3d expSO(const Matrix3d& m)
    /** exp for rotational matrix using Rodrigues' formula
     * @param m rotational matrix
     * @return exp(m)
     */
    {
        assert( ((m + m.transpose())).squaredNorm() < EPSILON );
        double norm2=m(0,1)*m(0,1) + m(0,2)*m(0,2) + m(1,2)*m(1,2);
        if(norm2<EPSILON){
            return E;
        }else{
            double norm = sqrtf(norm2);
            return E + sinf(norm)/norm * m + (1.0f-cosf(norm))/norm2 * m*m;
        }
    }
    
    Matrix4d logSEc(const Matrix4d& mm, double& prevTheta, Vector3d& prevN)
    /** "Continuous" log for rigid transformation (screw) matrix
     * @param mm rigid transformation matrix
     * @param prevTheta rotational angle
     * @param prevN rotational axis vector
     * @return branch of log(m) with axis and angle closest to the given ones
     */
    {
        Matrix3d m = mm.block(0,0,3,3);
        assert( ((m * m.transpose()) - E).squaredNorm() < EPSILON );
        Vector3d v(mm(3,0), mm(3,1), mm(3,2));
        Matrix3d X = logSOc(m, prevTheta, prevN);
        Matrix3d A;
        if(fabs(sinf(prevTheta))<EPSILON){
            A = E;
        }else{
            A = E - 0.5f * X + (1.0f/(prevTheta*prevTheta) - 0.5f*(1.0f+cosf(prevTheta))/(sinf(prevTheta)*prevTheta)) * X * X;
        }
        return(affine0(X, A.transpose()*v));
    }
    
    Matrix4d expSE(const Matrix4d& mm)
    /** exp for rigid transformation (screw) matrix
     * @param mm rigid transformation matrix
     * @return exp(mm)
     */
    {
        Matrix3d m = mm.block(0,0,3,3);
        assert( ((m + m.transpose())).squaredNorm() < EPSILON );
        Vector3d v;
        v << mm(3,0), mm(3,1), mm(3,2);
        double norm2=m(0,1)*m(0,1) + m(0,2)*m(0,2) + m(1,2)*m(1,2);
        Matrix3d A,ans;
        if(norm2<EPSILON){
            return affine(E,v);
        }else{
            double norm = sqrtf(norm2);
            ans = E + sinf(norm)/norm * m + (1.0f-cosf(norm))/norm2 * m*m;
            A = E + (1.0f-cosf(norm))/norm2 * m + (norm-sinf(norm))/(norm*norm2) * m*m ;
            //        A = sinf(norm)/norm * E + (1.0f-cosf(norm))/norm2 * ans;
            return affine(ans, A.transpose()*v);
        }
    }
    
    
    Matrix3d logDiag(const Matrix3d& U, const Vector3d& s)
    /** log for diagonalised matrix  m = U diag(s) U^t with positive eigenvalues
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
    /** exp for diagonalised matrix  m = U diag(s) U^{-1}
     * @param U diagonalising matrix
     * @param s eigenvalues
     * @return exp(m)
     */
    {
        Vector3d d = s.array().exp();
        return(U* d.asDiagonal() *U.inverse());
    }
    
    Matrix3d expSymD(const Matrix3d& m)
    /** exp for symmetric matrix by diagonalization (slower than expSym)
     * @param m symmetric matrix
     * @return exp(m)
     */
    {
        assert( ((m - m.transpose())).squaredNorm() < EPSILON );
        if (m.squaredNorm() < EPSILON){
            return E+m;
        }
        SelfAdjointEigenSolver<Matrix3d> eigensolver;
        eigensolver.computeDirect(m);
        Vector3d s(eigensolver.eigenvalues());
        Matrix3d U(eigensolver.eigenvectors());
        s = s.array().exp();
        return(U * s.asDiagonal() * U.transpose());
    }
    
    Matrix3d logSymD(const Matrix3d& m)
    /** log for symmetric matrix by diagonalization (slower than expSym)
     * @param m symmetric matrix
     * @return log(m)
     */
    {
        assert( ((m - m.transpose())).squaredNorm() < EPSILON );
        if ((m-E).squaredNorm() < EPSILON){
            return m-E;
        }
        SelfAdjointEigenSolver<Matrix3d> eigensolver;
        eigensolver.computeDirect(m);
        Vector3d s(eigensolver.eigenvalues());
        Matrix3d U(eigensolver.eigenvectors());
        s = s.array().log();
        return(U * s.asDiagonal() * U.transpose());
    }
    
    Matrix3d expSym(const Matrix3d& m, Vector3d e=Vector3d::Zero())
    /** exp for symmetric matrix by spectral decomposition
     * @param m symmetric matrix
     * @param e if eigenvalues of m are given, we use it
     * @return exp(m)
     */
    {
        assert( ((m - m.transpose())).squaredNorm() < EPSILON );
        if (m.squaredNorm() < EPSILON){
            return E + m;
        }
        if(e == Vector3d::Zero()){
            // compute eigenvalues
            // eigenvalues are sorted in increasing order.
            SelfAdjointEigenSolver<Matrix3d> eigensolver;
            eigensolver.computeDirect(m, EigenvaluesOnly);
            e = eigensolver.eigenvalues();
        }
        double a, b, c;
        double e12 = e[0] - e[1];
        double e23 = e[1] - e[2];
        double e13 = e[0] - e[2];
        // when some eigenvalues conside
        if(fabsf(e12)<EPSILON){
            if(fabsf(e23)<EPSILON){
                return expf(e[1]) * E;
            }else{
                a=expf(e[1]) / e23;
                b=expf(e[2]) / e23;
                return (a-b)*m + (e[1]*b-e[2]*a)*E;
            }
        }else{
            if(fabsf(e23)<EPSILON){
                a=expf(e[0]) / e12;
                b=expf(e[1]) / e12;
                return (a-b)*m + (e[0]*b-e[1]*a)*E;
            }
        }
        // when all eigenvalues are distinct
        a = expf(e[0]) / (e12*e13);
        b = -expf(e[1]) / (e23*e12);
        c = expf(e[2]) / (e13*e23);
        return (a + b + c) * m * m
        - (a * (e[1] + e[2]) + b * (e[2] + e[0]) + c * (e[0] + e[1])) * m
        + (a * e[1] * e[2] + b * e[2] * e[0] + c * e[0] * e[1]) * E;
    }
    
    Matrix3d logSym(const Matrix3d& m, Vector3d& lambda)
    /** log for positive symmetric matrix by spectral decomposition
     * @param m symmetric matrix
     * @param lambda returns eigen values for log(m)
     * @return log(m)
     */
    {
        assert( ((m - m.transpose())).squaredNorm() < EPSILON );
        if ((m-E).squaredNorm() < EPSILON){
            return m-E;
        }
        // compute eigenvalues only
        // eigenvalues are sorted in increasing order.
        SelfAdjointEigenSolver<Matrix3d> eigensolver;
        eigensolver.computeDirect(m, EigenvaluesOnly);
        Vector3d e;
        e = eigensolver.eigenvalues();
        assert(e[0] > 0 && e[1] > 0 && e[2] > 0);
        double a, b, c;
        double e12 = e[0] - e[1];
        double e23 = e[1] - e[2];
        double e13 = e[0] - e[2];
        lambda = e.array().log();
        // when some eigenvalues conside
        if(fabsf(e12)<EPSILON){
            if(fabsf(e23)<EPSILON){
                return lambda(0) * E;
            }else{
                a= lambda(1) / e23;
                b= lambda(2) / e23;
                return (a-b)*m + (e[1]*b-e[2]*a)*E;
            }
        }else{
            if(fabsf(e23)<EPSILON){
                a=lambda(0) / e12;
                b=lambda(1) / e12;
                return (a-b)*m + (e[0]*b-e[1]*a)*E;
            }
        }
        // when all eigenvalues are distinct
        a = lambda(0) / (e12*e13);
        b = -lambda(1) / (e23*e12);
        c = lambda(2) / (e13*e23);
        return (a + b + c) * m * m
        - (a * (e[1] + e[2]) + b * (e[2] + e[0]) + c * (e[0] + e[1])) * m
        + (a * e[1] * e[2] + b * e[2] * e[0] + c * e[0] * e[1]) * E;
    }
    
    
    Matrix3d frechetSO(const std::vector<Matrix3d> &m, const std::vector<double> &w, const int max_step=10)
    /** Frechet sum for rotations
     * @param m array of rotation matrices to be averaged
     * @param w array of weights
     * @param max_step max steps for iteration
     * @return weighted Frechet sum
     */
    {
        assert(m.size() == w.size());
        if(m.empty()) return(E);
        Matrix3d Z = m[0];
        for(int i=0;i<max_step;i++){
            Matrix3d W = Matrix3d::Zero();
            Matrix3d ZI = Z.transpose();
            for(int j=0;j<m.size();j++){
                W += w[j] * logSO(ZI * m[j]);
            }
            if(W.squaredNorm()<EPSILON) break;
            Z = Z * expSO(W);
        }
        return(Z);
    }
    
    Matrix4d frechetSE(const std::vector<Matrix4d> &m, const std::vector<double> &w, const int max_step=10)
    /** Frechet sum for rigid transformations
     * @param m array of rigid transformation matrices to be averaged
     * @param w array of weights
     * @param max_step max steps for iteration
     * @return weighted Frechet sum
     */
    {
        assert(m.size() == w.size());
        if(m.empty()) return(Matrix4d::Identity());
        Matrix4d Z = m[0];
        double theta = 0.0;
        Vector3d prevN = Vector3d::Zero();
        for(int i=0;i<max_step;i++){
            Matrix4d W = Matrix4d::Zero();
            Matrix4d ZI = Z.inverse();
            for(int j=0;j<m.size();j++){
                W += w[j] * logSEc(ZI * m[j], theta, prevN);
            }
            if(W.squaredNorm()<EPSILON) break;
            Z = Z * expSE(W);
        }
        return(Z);
    }
    
    Vector3d eigenvaluesSym(const Matrix3d &m)
    /** (Obsolete) Eigenvalues for symmetric matrix using Viete's formula
     * use Eigen's "SelfAdjointEigenSolver< MatrixType > & computeDirect" instead
     * @param m symmetric matrix
     * @return eigenvalues
     */
    {
        assert( ((m - m.transpose())).squaredNorm() < EPSILON );
        Vector3d w;
        double s1,s2,s3;
        // check if m is diagonal;
        if (m(0,1)*m(0,1)+m(0,2)*m(0,2)+m(1,2)*m(1,2)<EPSILON){
            s1 = m(0,0);
            s2 = m(1,1);
            s3 = m(2,2);
        }else{
            double a = - m.trace();
            double aa = a*a;
            double b = (aa - m.squaredNorm()) / 2;
            
            double p = aa / 3 - b;
            double q = 2 * aa * a / 27 - a * b / 3 - m.determinant();
            
            // since p,q are close to 0, the next statement is numerically unstable.
            if (p==0) {
                s1 = s2 = s3 = -a/3;
            }else {
                double r = sqrtf(4 * p / 3);
                double k = - 4 * q / (r * r * r);
                double theta;
                if (k>1.0){
                    theta = 0.0;
                }else if(k<-1.0){
                    theta = M_PI;
                }else{
                    theta = acosf(k);
                }
                s1 = r*cosf(theta / 3) - a / 3;
                s2 = r*cosf((theta +  2 * M_PI) / 3) - a / 3;
                s3 = r*cosf((theta + 4 * M_PI) / 3) - a / 3;
            }
        }
        // sort
        if (s1 < s2) {
            if (s3 < s1) std::swap(s1,s3);
        } else {
            if (s2 < s3) std::swap(s1,s2);
            else std::swap(s1,s3);
        }
        if(s3<s2) std::swap(s2,s3);
        w << s1, s2, s3;
        return w;
    }
    
    void polar(const Matrix3d& m, Matrix3d& U, Vector3d& s, Matrix3d& R)
    /** Polar decomposition m = U diag(s) U^T R
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
        Vector3d lambda=Vector3d::Zero();
        Matrix3d logS = logSym(m*m.transpose(), lambda)/2.0;
        S = expSym(logS, lambda/2.0);
        R = expSym(-logS, -lambda/2.0) * m;
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
}

