/**
 * @file cageDeformerARAP.h
 * @brief Cage Deformer with ARAP mod plugin for Maya
 * @section LICENSE The MIT License
 * @section requirements:  Eigen 3:  http://eigen.tuxfamily.org/
 * @section Autodesk Maya: http://www.autodesk.com/products/autodesk-maya/overview
 * @section (included) AffineLib: https://github.com/shizuo-kaji/AffineLib
 * @section limitation: the target mesh needs to be enough "clean" for ARAP
 * @version 0.15
 * @date  20/Apr/2014
 * @author Shizuo KAJI
 */

#pragma once

#pragma comment(linker, "/export:initializePlugin /export:uninitializePlugin")

#include <maya/MFnPlugin.h>

#include <numeric>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Sparse>

#include "affinelib.h"
#include "tetrise.h"
#include "deformerConst.h"
#include "MeshMaya.h"
#include "ARAP.h"

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

using namespace Eigen;

class CageDeformerNode : public MPxDeformerNode{
public:
    CageDeformerNode(): isError(0) {};
    virtual MStatus deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex );
    static  void*   creator();
    static  MStatus initialize();
 
    static MTypeId      id;
    static MString      nodeName;
    static MObject      aARAP;
    static MObject      aCageMesh;
    static MObject      aTetMode;
    static MObject      aCageMode;
    static MObject      aBlendMode;
	static MObject		aRotationConsistency;
	static MObject		aFrechetSum;
    static MObject      aConstraintMode;
    static MObject      aConstraintWeight;
    static MObject      aNormExponent;
    static MObject      aNormaliseTet;
    static MObject      aTransWeight;
    static MObject      aSymmetricFace;
    static MObject      aWeightMode;
    static MObject      aConstraintRadius;
    static MObject      aMaxDist;
    static MObject      aIteration;
    
private:
	void tetMatrixC(const MPointArray& p, const MIntArray& triangles, std::vector<Matrix4d>& m);
	void arapHI(const std::vector<Matrix4d>& PI, const std::vector<int>& tetList);
	void arapG(const std::vector< Matrix4d>& At, const std::vector<Matrix4d>& PI,
               const std::vector<int>& tetList, const std::vector< Matrix4d>& aff, MatrixXd& G);
    MObject     initCageMesh;
	int dim;
    short isError;
    std::vector<Matrix4d> logSE;   // for rotation consistency
    std::vector<Matrix3d> logR;   // for rotation consistency
    std::vector< std::vector<double> > w;   // weight
    // cage
    std::vector<Matrix4d> cageMatrixI;   // inverses of initial cage matrix
    std::vector<int> cageFaceList;
    std::vector<Vector3d> cageTetCenter;
    std::vector< int > cageTetList;
    std::vector< edge > cageEdgeList;
    std::vector<vertex> cageVertexList;
    std::vector<Vector3d> initCagePts;
    // ARAP
	std::vector<Matrix4d> PI;      // inverse of tet matrix
    std::vector< int > tetList;
    std::vector< edge > edgeList;
    std::vector<vertex> vertexList;
    std::vector<int> faceList;
    std::vector<Vector3d> pts;
    std::vector<double> tetWeight;
    std::vector<double> cageTetWeight;
    std::vector< std::map<int,double> > constraint;
    SpSolver solver;   // ARAP solver
    SpMat constraintMat;                // ARAP constraint matrix
    //  find affine transformations for tetrahedra
    std::vector<Matrix4d> cageMatrix, SE, logAff,Aff;
    std::vector<Matrix3d> R,logS,S,logGL;
    std::vector<Vector3d> L;
    std::vector<Vector4d> quat;
    std::vector<Matrix4d> A,blendedSE;
    std::vector<Matrix3d> blendedR, blendedS;
    std::vector<Vector4d> blendedL;
    std::vector<Vector3d> new_pts;
    std::vector<Matrix4d> Q;  // temporary
};

