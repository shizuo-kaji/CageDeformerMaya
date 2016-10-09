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

#include <map>
#include <Eigen/Sparse>

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


#define ERROR_ARAP_PRECOMPUTE 1

// ARAP precompute with soft constraint
int ARAPprecompute(const std::vector<Matrix4d>& tetMatrixInverse, const std::vector<int>& tetList,
                    const std::vector<double>& tetWeight, std::vector< std::map<int,double> >& constraint,
                                   double transWeight, int dim, SpMat& constraintMat, SpSolver& solver){
    int numTet = (int)tetList.size()/4;
    std::vector<T> tripletListMat(0);
    tripletListMat.reserve(numTet*16);
    Matrix4d Hlist;
	Matrix4d diag=Matrix4d::Identity();
	diag(3,3)=transWeight;
	for(int i=0;i<numTet;i++){
		Hlist=tetWeight[i] * tetMatrixInverse[i].transpose()*diag*tetMatrixInverse[i];
		for(int j=0;j<4;j++){
			for(int k=0;k<4;k++){
                tripletListMat.push_back(T(tetList[4*i+j],tetList[4*i+k],Hlist(j,k)));
			}
		}
	}
    SpMat mat(dim, dim);
    mat.setFromTriplets(tripletListMat.begin(), tripletListMat.end());
    // set soft constraint
    int cur=0;
    tripletListMat.resize(0);
    std::vector<T> tripletListF(0);
    std::map<int, double>::iterator iter;
    for(int i=0;i<constraint.size();i++){
        for(iter = constraint[i].begin(); iter != constraint[i].end(); iter++){
            tripletListMat.push_back(T( iter->first, cur, iter->second));
            tripletListF.push_back(T( cur, iter->first, 1));
            cur++;
        }
    }
    constraintMat.resize(dim,cur);
    constraintMat.setFromTriplets(tripletListMat.begin(), tripletListMat.end());
    SpMat F(cur,dim);
    F.setFromTriplets(tripletListF.begin(), tripletListF.end());
    // mat = (H^T,C_M)*(H \\ C_F),   C_M = constraintWeight * C_F^T
    mat += numTet * constraintMat * F;
    solver.compute(mat);
    if(solver.info() != Success){
        //std::string error_mes = solver.lastErrorMessage();
        MGlobal::displayInfo("Cleanup the mesh first: Mesh menu => Cleanup => Remove zero edges, faces");
        return ERROR_ARAP_PRECOMPUTE;
    }
    return 0;
}

//
void ARAPSolve(const std::vector<Matrix4d>& A, const std::vector<Matrix4d>& tetMatrixInverse,
                                  const std::vector<int>& tetList,
               const std::vector<double>& tetWeight, const std::vector<Vector3d>& constraintVector,
               double transWeight, int dim, const SpMat& constraintMat, const SpSolver& solver,
               MatrixXd& Sol){
    int numTet = (int)tetList.size()/4;
    Matrix4d Glist;
    Matrix4d diag=Matrix4d::Identity();
    diag(3,3)=transWeight;
    MatrixXd G = MatrixXd::Zero(dim,3);
    for(int i=0;i<numTet;i++){
        Glist= tetWeight[i] * tetMatrixInverse[i].transpose() * diag * A[i];
        for(int k=0;k<3;k++){
            for(int j=0;j<4;j++){
                G(tetList[4*i+j],k) += Glist(j,k);
            }
        }
    }
    // set soft constraint
    std::vector<T> constraintList(0);
    int numConstraint = (int) constraintVector.size();
    constraintList.reserve(numConstraint*3);
    for(int i=0;i<numConstraint;i++){
        constraintList.push_back(T(i,0,constraintVector[i][0]));
        constraintList.push_back(T(i,1,constraintVector[i][1]));
        constraintList.push_back(T(i,2,constraintVector[i][2]));
    }
    SpMat S(numConstraint,3);
    S.setFromTriplets(constraintList.begin(), constraintList.end());
    // FS = (H^T,C_M) * (G \\ S)
    SpMat FS = numTet * constraintMat * S;
    G += MatrixXd(FS);
    Sol = solver.solve(G);
}

// ARAP precompute with approx. hard constraint
int ARAPprecomputeHard(const std::vector<Matrix4d>& tetMatrixInverse, const std::vector<int>& tetList,
                   const std::vector<double>& tetWeight, std::vector< std::map<int,double> >& constraint,
                   double transWeight, int dim, SpMat& constraintMat, SpSolver& solver){
    int numTet = (int)tetList.size()/4;
    std::vector<T> tripletListMat(0);
    tripletListMat.reserve(numTet*16);
    Matrix4d Hlist;
	Matrix4d diag=Matrix4d::Identity();
	diag(3,3)=transWeight;
	for(int i=0;i<numTet;i++){
		Hlist=tetWeight[i] * tetMatrixInverse[i].transpose()*diag*tetMatrixInverse[i];
		for(int j=0;j<4;j++){
			for(int k=0;k<4;k++){
                tripletListMat.push_back(T(tetList[4*i+j],tetList[4*i+k],Hlist(j,k)));
			}
		}
	}
    // set hard constraint
    int cur=0;
    tripletListMat.resize(0);
    std::vector<T> tripletListF(0);
    std::map<int, double>::iterator iter;
    for(int i=0;i<constraint.size();i++){
        for(iter = constraint[i].begin(); iter != constraint[i].end(); iter++){
            tripletListMat.push_back(T( iter->first, dim+cur, 1));
            tripletListF.push_back(T( dim+cur, iter->first, 1));
            cur++;
        }
    }
    SpMat mat(dim+cur,dim+cur);
    mat.setFromTriplets(tripletListMat.begin(), tripletListMat.end());
    //
    solver.compute(mat);
    if(solver.info() != Success){
        //        std::string error_mes = solver.lastErrorMessage();
        MGlobal::displayInfo("Cleanup the mesh first: Mesh menu => Cleanup => Remove zero edges, faces");
        return ERROR_ARAP_PRECOMPUTE;
    }
    return 0;
}


/// harmonic weighting
int harmonicWeight(int numPts, const std::vector<Matrix4d>& tetMatrix, const std::vector<int>& tList, const std::vector<double>& tetWeight,
                   std::vector< std::map<int,double> >& weightConstraint, std::vector<double>& weightConstraintValue,
                   std::vector< std::vector<double> >& ptsWeight){
    int num = (int)tList.size()/4;
    int dim = numPts + num;
    int numPrb = (int)weightConstraint.size();
    std::vector<Matrix4d> PI(num);
    for(int i=0;i<num;i++){
        PI[i]=tetMatrix[i].inverse().eval();
    }
    // LHS
    SpSolver weightSolver;
    SpMat weightConstraintMat;
    if(ARAPprecompute(PI, tList, tetWeight, weightConstraint, 0, dim, weightConstraintMat, weightSolver)){
        MGlobal::displayInfo("An error occurred in Harmonic field computation.");
        return 1;
    }
    // RHS
    std::vector<T> tripletListMat(0);
    std::map<int, double>::iterator iter;
    int cur=0;
    for(int i=0;i<numPrb;i++){
        for(iter = weightConstraint[i].begin(); iter != weightConstraint[i].end(); iter++){
            tripletListMat.push_back(T( cur, i, weightConstraintValue[cur]));
            cur++;
        }
    }
    SpMat F(cur,numPrb);
    F.setFromTriplets(tripletListMat.begin(), tripletListMat.end());
    SpMat G = num * weightConstraintMat * F;
    // solve
    SpMat Sol = weightSolver.solve(G);
    ptsWeight.resize(numPrb);
    for (int i=0;i<numPrb; i++){
        ptsWeight[i].resize(numPts);
        for(int j=0;j<numPts; j++){
            ptsWeight[i][j] = Sol.coeff(j,i);
        }
    }
    return 0;
}

/// harmonic weighting with cotan laplacian
int harmonicCotan(int numPts, const std::vector<Matrix4d>& P, const std::vector<int>& tetList, const std::vector<double>& tetWeight,
                   std::vector< std::map<int,double> >& constraint, std::vector<double>& constraintValue,
                   std::vector< std::vector<double> >& ptsWeight){
    int num = (int)tetList.size()/4;
    int numPrb = (int)constraint.size();
    // LHS
    std::vector<T> tripletListMat(0);
    tripletListMat.reserve(num*9);
    for(int i=0;i<num;i++){
        std::vector<Vector3d> l(3);
        l[0] << P[i](2,0)-P[i](1,0), P[i](2,1)-P[i](1,1), P[i](2,2)-P[i](1,2);
        l[1] << P[i](0,0)-P[i](2,0), P[i](0,1)-P[i](2,1), P[i](0,2)-P[i](2,2);
        l[2] << P[i](1,0)-P[i](0,0), P[i](1,1)-P[i](0,1), P[i](1,2)-P[i](0,2);
        double w = tetWeight[i]/(l[1].cross(l[2]).norm());
        for(int j=0;j<3;j++){
            int j1=(j+1) % 3;
            int j2=(j+2) % 3;
            tripletListMat.push_back(T(tetList[4*i+j],tetList[4*i+j],   w*l[j].dot(l[j1]+l[j2])));
            tripletListMat.push_back(T(tetList[4*i+j],tetList[4*i+j1], -w*l[j].dot(l[j2])));
            tripletListMat.push_back(T(tetList[4*i+j],tetList[4*i+j2], -w*l[j].dot(l[j1])));
        }
    }
    SpMat mat(numPts, numPts);
    mat.setFromTriplets(tripletListMat.begin(), tripletListMat.end());
    mat = mat.transpose() * mat;
    // set soft constraint
    int cur=0;
    std::vector<T> tripletListC(0), tripletListF(0), tripletListFR(0);
    std::map<int, double>::iterator iter;
    for(int i=0;i<constraint.size();i++){
        for(iter = constraint[i].begin(); iter != constraint[i].end(); iter++){
            tripletListC.push_back(T( iter->first, cur, iter->second));
            tripletListF.push_back(T( cur, iter->first, 1));
            tripletListFR.push_back(T( cur, i, constraintValue[cur]));
            cur++;
        }
    }
    SpMat constraintMat(numPts,cur);
    constraintMat.setFromTriplets(tripletListC.begin(), tripletListC.end());
    SpMat F(cur,numPts);
    F.setFromTriplets(tripletListF.begin(), tripletListF.end());
    mat += num * constraintMat * F;
    SpSolver solver;
    solver.compute(mat);
    if(solver.info() != Success){
        //std::string error_mes = solver.lastErrorMessage();
        MGlobal::displayInfo("Cleanup the mesh first: Mesh menu => Cleanup => Remove zero edges, faces");
        return ERROR_ARAP_PRECOMPUTE;
    }
    // RHS
    SpMat FR(cur,numPrb);
    FR.setFromTriplets(tripletListFR.begin(), tripletListFR.end());
    SpMat G = num * constraintMat * FR;
    // solve
    SpMat Sol = solver.solve(G);
    ptsWeight.resize(numPrb);
    for (int i=0;i<numPrb; i++){
        ptsWeight[i].resize(numPts);
        for(int j=0;j<numPts; j++){
            ptsWeight[i][j] = Sol.coeff(j,i);
        }
    }
    return 0;
}



