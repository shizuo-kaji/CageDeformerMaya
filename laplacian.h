/**
 * @file ARAP.h
 * @brief ARAP related functions
 * @section LICENSE The MIT License
 * @section requirements:  Eigen library
 * @version 0.10
 * @date  Aug. 2014
 * @author Shizuo KAJI
 */

#pragma once

#include <utility>
#include <Eigen/Sparse>

#include "deformerConst.h"

//#define _SuiteSparse
//#define _CERES

using namespace Eigen;

typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;

#ifdef _SuiteSparse
#include <Eigen/CholmodSupport>
typedef CholmodDecomposition<SpMat> SpSolver;
#else
typedef SimplicialLDLT<SpMat> SpSolver;
//typedef SimplicialLLT<SpMat> SpSolver;
//typedef SparseLU<SpMat> SpSolver;
#endif

#ifdef _CERES
#include "ceres/ceres.h"
#include "glog/logging.h"
#endif



class Laplacian {
public:
    int numTet;  // the number of tetrahedra
    int dim;   // the dimension of the system including ghost vertices
    double transWeight;
    SpSolver solver;
    SpMat constraintMat;
    SpMat laplacian;
    std::vector<int> tetList;
    std::vector<Matrix4d> tetMatrix,tetMatrixInverse;
    std::vector<double> tetWeight;
    std::vector< std::pair<int,double> > constraintWeight;  //  [i,w] = i-th vertex is constrained with weight w
    MatrixXd constraintVal;       // i-th row = value of i-th constraint
    MatrixXd Sol;
    Laplacian(): numTet(0), tetMatrix(0), tetMatrixInverse(0), tetWeight(0), constraintWeight(0), transWeight(0) {
    };
    int ARAPprecompute();
    void ARAPSolve(const std::vector<Matrix4d>& targetMat);
    void harmonicSolve();
    int cotanPrecompute();
    void computeTetMatrixInverse();
};


// construct the system of ARAP with soft constraints
int Laplacian::ARAPprecompute(){
    std::vector<T> tripletListMat(0);
    tripletListMat.reserve(numTet*16);
    Matrix4d Hlist;
    Matrix4d diag=Matrix4d::Identity();
    diag(3,3)=transWeight;
    for(int i=0;i<numTet;i++){
        Hlist=tetWeight[i] * tetMatrixInverse[i].transpose() * diag * tetMatrixInverse[i];
        for(int j=0;j<4;j++){
            for(int k=0;k<4;k++){
                tripletListMat.push_back(T(tetList[4*i+j],tetList[4*i+k],Hlist(j,k)));
            }
        }
    }
    SpMat mat(dim,dim);
    mat.setFromTriplets(tripletListMat.begin(), tripletListMat.end());
    // set soft constraint
    int numConstraints = constraintWeight.size();
    std::vector<T> tripletListF(0),tripletListC(0);
    for(int i=0;i<numConstraints;i++){
        tripletListC.push_back(T( constraintWeight[i].first, i, constraintWeight[i].second));
        tripletListF.push_back(T( i, constraintWeight[i].first, 1));
    }
    constraintMat.resize(dim,numConstraints);
    constraintMat.setZero();
    constraintMat.setFromTriplets(tripletListC.begin(), tripletListC.end());
    SpMat F(numConstraints,dim);
    F.setFromTriplets(tripletListF.begin(), tripletListF.end());
    // mat = (L^T,C_M)*(L \\ C_F),   C_M = constraintWeight * C_F^T
    mat += numTet * constraintMat * F;
    solver.compute(mat);
    if(solver.info() != Success){
        //std::string error_mes = solver.lastErrorMessage();
        MGlobal::displayInfo("Cleanup the mesh first: Mesh menu => Cleanup => Remove zero edges, faces");
        return ERROR_ARAP_PRECOMPUTE;
    }
    return 0;
}

// solve the ARAP system
void Laplacian::ARAPSolve(const std::vector<Matrix4d>& targetMat){
    Matrix4d Glist;
    Matrix4d diag=Matrix4d::Identity();
    diag(3,3)=transWeight;
    MatrixXd G = MatrixXd::Zero(dim,3);
    for(int i=0;i<numTet;i++){
        Glist= tetWeight[i] * tetMatrixInverse[i].transpose() * diag * targetMat[i];
        for(int k=0;k<3;k++){
            for(int j=0;j<4;j++){
                G(tetList[4*i+j],k) += Glist(j,k);
            }
        }
    }
    // set soft constraint
    // (H^T,C_M) * (G \\ constraintVal)
    G += numTet * constraintMat * constraintVal;
    Sol = solver.solve(G);
}

// harmonic weighting
void Laplacian::harmonicSolve(){
    MatrixXd G = numTet * constraintMat * constraintVal;
    Sol = solver.solve(G);
}

// harmonic weighting with cotan laplacian
int Laplacian::cotanPrecompute(){
    std::vector<T> tripletListMat(0);
    tripletListMat.reserve(numTet*9);
    for(int i=0;i<numTet;i++){
        std::vector<Vector3d> l(3);
        l[0] << tetMatrix[i](2,0)-tetMatrix[i](1,0), tetMatrix[i](2,1)-tetMatrix[i](1,1), tetMatrix[i](2,2)-tetMatrix[i](1,2);
        l[1] << tetMatrix[i](0,0)-tetMatrix[i](2,0), tetMatrix[i](0,1)-tetMatrix[i](2,1), tetMatrix[i](0,2)-tetMatrix[i](2,2);
        l[2] << tetMatrix[i](1,0)-tetMatrix[i](0,0), tetMatrix[i](1,1)-tetMatrix[i](0,1), tetMatrix[i](1,2)-tetMatrix[i](0,2);
        double w = tetWeight[i]/(l[1].cross(l[2]).norm());
        for(int j=0;j<3;j++){
            int j1=(j+1) % 3;
            int j2=(j+2) % 3;
            tripletListMat.push_back(T(tetList[4*i+j],tetList[4*i+j],   w*l[j].dot(l[j1]+l[j2])));
            tripletListMat.push_back(T(tetList[4*i+j],tetList[4*i+j1], -w*l[j].dot(l[j1])));
            tripletListMat.push_back(T(tetList[4*i+j],tetList[4*i+j2], -w*l[j].dot(l[j2])));
        }
    }
    laplacian.resize(dim, dim);
    laplacian.setZero();
    laplacian.setFromTriplets(tripletListMat.begin(), tripletListMat.end());
    // set soft constraint
    int numConstraints = constraintWeight.size();
    std::vector<T> tripletListF(0),tripletListC(0);
    for(int i=0;i<numConstraints;i++){
        tripletListC.push_back(T( constraintWeight[i].first, i, constraintWeight[i].second));
        tripletListF.push_back(T( i, constraintWeight[i].first, 1));
    }
    constraintMat.resize(dim,numConstraints);
    constraintMat.setZero();
    constraintMat.setFromTriplets(tripletListC.begin(), tripletListC.end());
    SpMat F(numConstraints,dim);
    F.setFromTriplets(tripletListF.begin(), tripletListF.end());
    SpMat mat = laplacian.transpose() * laplacian + numTet * constraintMat * F;
    solver.compute(mat);
    if(solver.info() != Success){
        //std::string error_mes = solver.lastErrorMessage();
        MGlobal::displayInfo("Cleanup the mesh first: Mesh menu => Cleanup => Remove zero edges, faces");
        return ERROR_ARAP_PRECOMPUTE;
    }
    return 0;
}

void Laplacian::computeTetMatrixInverse(){
    tetMatrixInverse.resize(numTet);
    for(int i=0;i<numTet;i++){
        tetMatrixInverse[i] = tetMatrix[i].inverse().eval();
    }
}







