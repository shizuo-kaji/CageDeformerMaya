/**
 * @file affinelib.h
 * @brief Library for 3D affine transformation
 * @section LICENSE The MIT License
 * @section  requirements:  Eigen library,   (optional) MKL
 * @section  CAUTION: the convention here is different from the one in the paper;
 *           we assume that matrices act on row vectors by right multiplication.
 *           ( that is, everything is transposed compared to the paper )
 * @version 0.10
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
#define EPSILON 0.000001f
// Identity matrix
#define E Matrix3f::Identity()
// for debug print
#define PRINT_MAT(X) std::cout << #X << ":\n" << X << std::endl << std::endl

using namespace Eigen;
using namespace std;
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Matrix4f);

class Matrixlib{
public:
    Matrixlib() {};
    
    // Polar decomposition
    static  void polar(const Matrix3f& m, Matrix3f& U, Vector3f& s, Matrix3f& R);
    static void polarBySVD(const Matrix3f& m, Matrix3f& U, Vector3f& s, Matrix3f& R);
    static void polarByParam(const Matrix3f& m, Matrix3f& S, Matrix3f& R);
    static void parametriseGL(const Matrix3f& m, Matrix3f& logS, Matrix3f& R);
    static  Matrix3f logSO(const Matrix3f& m);
    static  Matrix3f logSOc(const Matrix3f& m, float& prevTheta, Vector3f& prevN);
    static  Matrix3f expSO(const Matrix3f& m);
    // log for Euclidian transformation
    static Matrix4f logSEc(const Matrix4f& m, float& prevTheta, Vector3f& prevN);
    // exp for Eucledian transformation
    static Matrix4f expSE(const Matrix4f& m);
    // log for diagonalised matrix
    static Matrix3f logDiag(const Matrix3f& U, const Vector3f& s);
    // exp for diagonalised matrix
    static Matrix3f expDiag(const Matrix3f& U, const Vector3f& s);
    // exp for symmetric matrix by diagonalization
    static  Matrix3f expSymD(const Matrix3f& mat);
    static Matrix3f logSymD(const Matrix3f& m);
    // exp for symmetric matrix by spectral decomp
    static  Matrix3f expSym(const Matrix3f& mat, Vector3f e=Vector3f::Zero());
    // log for symmetric matrix by spectral decomp
    static Matrix3f logSym(const Matrix3f& m, Vector3f& lambda);
    // Frechet sum for rotations
    static Matrix3f frechetSO(const std::vector<Matrix3f> &m, const std::vector<float> &w, const int max_step=10);
    // Frechet sum for rigid transformations
    static Matrix4f frechetSE(const std::vector<Matrix4f> &m, const std::vector<float> &w, const int max_step=10);
	// (obsolete) eigenvalues for symmetric matrix by Viete
	static Vector3f eigenvaluesSym(const Matrix3f &mat);
    // compose affine matrix from linear matrix and translation vector
    static Matrix4f affine(const Matrix3f& m, const Vector3f l=Vector3f::Zero());
    // compose 4x4-matrix in log space
    static Matrix4f affine0(const Matrix3f& m, const Vector3f l=Vector3f::Zero());
    static Vector3f transPart(const Matrix4f& m);

private:
};

void Matrixlib::polar(const Matrix3f& m, Matrix3f& U, Vector3f& s, Matrix3f& R)
/** Polar decomposition m = U diag(s) U^T R
 * @param m matrix to be decomposed
 * @param U diagonaliser of symmetric part
 * @param s singular values
 * @param R rotation part
 */
{
    assert(m.determinant()>0);
    Matrix3f A= m*m.transpose();
	SelfAdjointEigenSolver<Matrix3f> eigensolver;
	eigensolver.computeDirect(A);
    s = eigensolver.eigenvalues();
    U = Matrix3f(eigensolver.eigenvectors());
    s = s.array().sqrt();
    Vector3f si = s.array().inverse();
    R = U * si.asDiagonal() * U.transpose() * m;
}

void Matrixlib::polarBySVD(const Matrix3f& m, Matrix3f& U, Vector3f& s, Matrix3f& R){
    /** Polar decomposition m = U diag(s) U^T R  by SVD
     * @param m matrix to be decomposed
     * @param U diagonaliser of symmetric part
     * @param s singular values
     * @param R rotation part
     */
    JacobiSVD<Matrix3f> svd(m, ComputeFullU | ComputeFullV);
    U = svd.matrixU();
    s = svd.singularValues();
    R = svd.matrixU() * svd.matrixV().transpose();
}

void Matrixlib::polarByParam(const Matrix3f& m, Matrix3f& S, Matrix3f& R){
    /** Polar decomposition m = S R   by parametrisation map
     * @param m matrix to be decomposed
     * @param S shear part
     * @param R rotation part
     */
    Vector3f lambda=Vector3f::Zero();
    Matrix3f logS = Matrixlib::logSym(m*m.transpose(), lambda)/2.0;
    S = Matrixlib::expSym(logS, lambda/2.0);
    R = Matrixlib::expSym(-logS, -lambda/2.0) * m;
}

void Matrixlib::parametriseGL(const Matrix3f& m, Matrix3f& logS, Matrix3f& R)
/** Parametrisation map for GL(3)
 * @param m linear matrix to be decomposed
 * @param logS log of shear part
 * @param R rotation part
 */
{
    assert(m.determinant()>0);
    Vector3f lambda=Vector3f::Zero();
    logS = Matrixlib::logSym(m*m.transpose(), lambda)/2.0;
    R = Matrixlib::expSym(-logS, -lambda/2.0) * m;
}


Matrix3f Matrixlib::logSO(const Matrix3f& m)
/** Log for rotational matrix using Rodrigues' formula
 * @param m rotational matrix
 * @return primary value of log(m)
 */
{
    assert( ((m * m.transpose()) - E).squaredNorm() < EPSILON );
    float tr=(m(0,0)+m(1,1)+m(2,2)-1.0f)/2.0;
    float theta;
    Matrix3f ans=Matrix3f::Zero();
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

Matrix3f Matrixlib::logSOc(const Matrix3f& m, float& prevTheta, Vector3f& prevN)
/** "Continuous" log for rotational matrix
 * @param m rotational matrix
 * @param prevTheta rotational angle
 * @param prevN rotational axis vector
 * @return branch of log(m) with axis and angle closest to the given ones
 */
{
//    assert( ((m * m.transpose()) - E).squaredNorm() < EPSILON );
    float tr=(m(0,0)+m(1,1)+m(2,2)-1.0)/2.0;
    float theta;
    Matrix3f ans=Matrix3f::Zero();
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

Matrix3f Matrixlib::expSO(const Matrix3f& m)
/** exp for rotational matrix using Rodrigues' formula
 * @param m rotational matrix
 * @return exp(m)
 */
 {
    assert( ((m + m.transpose())).squaredNorm() < EPSILON );
    float norm2=m(0,1)*m(0,1) + m(0,2)*m(0,2) + m(1,2)*m(1,2);
    if(norm2<EPSILON){
        return E;
    }else{
        float norm = sqrtf(norm2);
        return E + sinf(norm)/norm * m + (1.0f-cosf(norm))/norm2 * m*m;
    }
}

Matrix4f Matrixlib::logSEc(const Matrix4f& mm, float& prevTheta, Vector3f& prevN)
/** "Continuous" log for rigid transformation (screw) matrix
 * @param mm rigid transformation matrix
 * @param prevTheta rotational angle
 * @param prevN rotational axis vector
 * @return branch of log(m) with axis and angle closest to the given ones
 */
{
    Matrix3f m = mm.block(0,0,3,3);
    assert( ((m * m.transpose()) - E).squaredNorm() < EPSILON );
    Vector3f v(mm(3,0), mm(3,1), mm(3,2));
    Matrix3f X = Matrixlib::logSOc(m, prevTheta, prevN);
    Matrix3f A;
    if(fabs(sinf(prevTheta))<EPSILON){
        A = E;
    }else{
        A = E - 0.5f * X + (1.0f/(prevTheta*prevTheta) - 0.5f*(1.0f+cosf(prevTheta))/(sinf(prevTheta)*prevTheta)) * X * X;
    }
	return(Matrixlib::affine0(X, A.transpose()*v));
}
    
Matrix4f Matrixlib::expSE(const Matrix4f& mm)
/** exp for rigid transformation (screw) matrix
 * @param mm rigid transformation matrix
 * @return exp(mm)
 */
{
    Matrix3f m = mm.block(0,0,3,3);
    assert( ((m + m.transpose())).squaredNorm() < EPSILON );
    Vector3f v;
    v << mm(3,0), mm(3,1), mm(3,2);
    float norm2=m(0,1)*m(0,1) + m(0,2)*m(0,2) + m(1,2)*m(1,2);
    Matrix3f A,ans;
    if(norm2<EPSILON){
        return Matrixlib::affine(E,v);
    }else{
        float norm = sqrtf(norm2);
        ans = E + sinf(norm)/norm * m + (1.0f-cosf(norm))/norm2 * m*m;
        A = E + (1.0f-cosf(norm))/norm2 * m + (norm-sinf(norm))/(norm*norm2) * m*m ;
//        A = sinf(norm)/norm * E + (1.0f-cosf(norm))/norm2 * ans;
        return Matrixlib::affine(ans, A.transpose()*v);
    }
}

    
    
Matrix3f Matrixlib::logDiag(const Matrix3f& U, const Vector3f& s)
/** log for diagonalised matrix  m = U diag(s) U^t with positive eigenvalues
 * @param U rotational diagonalising matrix
 * @param s eigenvalues (must be all positive)
 * @return log(m)
 */
{
    assert( s[0] > 0 && s[1] > 0 && s[2] > 0);
    Vector3f d = s.array().log();
    return(U* d.asDiagonal() *U.transpose());
}

Matrix3f Matrixlib::expDiag(const Matrix3f& U, const Vector3f& s)
/** exp for diagonalised matrix  m = U diag(s) U^{-1}
 * @param U diagonalising matrix
 * @param s eigenvalues
 * @return exp(m)
 */
{
    Vector3f d = s.array().exp();
    return(U* d.asDiagonal() *U.inverse());
}

Matrix3f Matrixlib::expSymD(const Matrix3f& m)
/** exp for symmetric matrix by diagonalization (slower than expSym)
 * @param m symmetric matrix
 * @return exp(m)
 */
{
    assert( ((m - m.transpose())).squaredNorm() < EPSILON );
    if (m.squaredNorm() < EPSILON){
        return E+m;
    }
	SelfAdjointEigenSolver<Matrix3f> eigensolver;
	eigensolver.computeDirect(m);
    Vector3f s(eigensolver.eigenvalues());
    Matrix3f U(eigensolver.eigenvectors());
    s = s.array().exp();
    return(U * s.asDiagonal() * U.transpose());
}

Matrix3f Matrixlib::logSymD(const Matrix3f& m)
/** log for symmetric matrix by diagonalization (slower than expSym)
 * @param m symmetric matrix
 * @return log(m)
 */
{
    assert( ((m - m.transpose())).squaredNorm() < EPSILON );
    if ((m-E).squaredNorm() < EPSILON){
        return m-E;
    }
	SelfAdjointEigenSolver<Matrix3f> eigensolver;
	eigensolver.computeDirect(m);
    Vector3f s(eigensolver.eigenvalues());
    Matrix3f U(eigensolver.eigenvectors());
    s = s.array().log();
    return(U * s.asDiagonal() * U.transpose());
}

Matrix3f Matrixlib::expSym(const Matrix3f& m, Vector3f e)
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
    if(e == Vector3f::Zero()){
        // compute eigenvalues
        // eigenvalues are sorted in increasing order.
        SelfAdjointEigenSolver<Matrix3f> eigensolver;
        eigensolver.computeDirect(m, EigenvaluesOnly);
        e = eigensolver.eigenvalues();
    }
	float a, b, c;
	float e12 = e[0] - e[1];
	float e23 = e[1] - e[2];
	float e13 = e[0] - e[2];
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

Matrix3f Matrixlib::logSym(const Matrix3f& m, Vector3f& lambda)
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
	SelfAdjointEigenSolver<Matrix3f> eigensolver;
	eigensolver.computeDirect(m, EigenvaluesOnly);
	Vector3f e;
	e = eigensolver.eigenvalues();
    assert(e[0] > 0 && e[1] > 0 && e[2] > 0);
	float a, b, c;
	float e12 = e[0] - e[1];
	float e23 = e[1] - e[2];
	float e13 = e[0] - e[2];
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


Matrix3f Matrixlib::frechetSO(const std::vector<Matrix3f> &m, const std::vector<float> &w, const int max_step)
/** Frechet sum for rotations
 * @param m array of rotation matrices to be averaged
 * @param w array of weights
 * @param max_step max steps for iteration
 * @return weighted Frechet sum
 */
{
    assert(m.size() == w.size());
    if(m.empty()) return(E);
    Matrix3f Z = m[0];
    for(int i=0;i<max_step;i++){
        Matrix3f W = Matrix3f::Zero();
        Matrix3f ZI = Z.transpose();
        for(int j=0;j<m.size();j++){
            W += w[j] * Matrixlib::logSO(ZI * m[j]);
        }
        if(W.squaredNorm()<EPSILON) break;
        Z = Z * expSO(W);
    }
    return(Z);
}

Matrix4f Matrixlib::frechetSE(const std::vector<Matrix4f> &m, const std::vector<float> &w, const int max_step)
/** Frechet sum for rigid transformations
 * @param m array of rigid transformation matrices to be averaged
 * @param w array of weights
 * @param max_step max steps for iteration
 * @return weighted Frechet sum
 */
{
    assert(m.size() == w.size());
    if(m.empty()) return(Matrix4f::Identity());
    Matrix4f Z = m[0];
    float theta = 0.0;
    Vector3f prevN = Vector3f::Zero();
    for(int i=0;i<max_step;i++){
        Matrix4f W = Matrix4f::Zero();
        Matrix4f ZI = Z.inverse();
        for(int j=0;j<m.size();j++){
            W += w[j] * Matrixlib::logSEc(ZI * m[j], theta, prevN);
        }
        if(W.squaredNorm()<EPSILON) break;
        Z = Z * expSE(W);
    }
    return(Z);
}

Vector3f Matrixlib::eigenvaluesSym(const Matrix3f &m)
/** (Obsolete) Eigenvalues for symmetric matrix using Viete's formula
 * use Eigen's "SelfAdjointEigenSolver< MatrixType > & computeDirect" instead
 * @param m symmetric matrix
 * @return eigenvalues
 */
{
    assert( ((m - m.transpose())).squaredNorm() < EPSILON );
    Vector3f w;
    float s1,s2,s3;
    // check if m is diagonal;
    if (m(0,1)*m(0,1)+m(0,2)*m(0,2)+m(1,2)*m(1,2)<EPSILON){
        s1 = m(0,0);
        s2 = m(1,1);
        s3 = m(2,2);
    }else{
        float a = - m.trace();
        float aa = a*a;
        float b = (aa - m.squaredNorm()) / 2;
        
        float p = aa / 3 - b;
        float q = 2 * aa * a / 27 - a * b / 3 - m.determinant();
        
        // since p,q are close to 0, the next statement is numerically unstable.
        if (p==0) {
            s1 = s2 = s3 = -a/3;
        }else {
            float r = sqrtf(4 * p / 3);
            float k = - 4 * q / (r * r * r);
            float theta;
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


Matrix4f Matrixlib::affine(const Matrix3f& m, const Vector3f l)
/** compose affine matrix from linear matrix and translation vector
 * @param m 3x3-matrix
 * @param l 3-dim translation vector
 * @return 4x4-affine transformation matrix
 */
 {
    Matrix4f aff;
    aff << m(0,0),m(0,1),m(0,2),0.0f,
    m(1,0),m(1,1),m(1,2),0.0f,
    m(2,0),m(2,1),m(2,2),0.0f,
    l(0),l(1),l(2),1.0f;
    return aff;
}

Matrix4f Matrixlib::affine0(const Matrix3f& m, const Vector3f l)
/** compose 4x4-matrix in log space
 * @param m 3x3-matrix
 * @param l 3-dim vector
 * @ return 4x4 matrix in log space
 */
{
    Matrix4f aff;
    aff << m(0,0),m(0,1),m(0,2),0.0f,
    m(1,0),m(1,1),m(1,2),0.0f,
    m(2,0),m(2,1),m(2,2),0.0f,
    l(0),l(1),l(2),0.0f;
    return aff;
}

Vector3f Matrixlib::transPart(const Matrix4f& m){
    return Vector3f(m(3,0),m(3,1),m(3,2));
}
