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

#pragma comment(lib, "OpenMaya.lib")
#pragma comment(lib, "OpenMayaRender.lib")
#pragma comment(lib, "OpenMayaAnim.lib")
#pragma comment(lib, "Foundation.lib")
#pragma comment(linker, "/export:initializePlugin /export:uninitializePlugin")

#include <maya/MFnPlugin.h>

#include <numeric>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>

#include "affinelib.h"

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

using namespace Eigen;

class CageDeformerNode : public MPxDeformerNode
{
public:
    CageDeformerNode(): cageMode(99) {};
    virtual MStatus deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex );
    static  void*   creator();
    static  MStatus initialize();
 
    static MTypeId      id;
    static MString      nodeName;
    static MObject      aCageMesh;
    static MObject      aCageMode;
    static MObject      aBlendMode;
	static MObject		aRotationConsistency;
	static MObject		aFrechetSum;
    static MObject      aConstraintMode;
    static MObject      aConstraintWeight;

private:
	void tetMatrix(const MPointArray& p, const MIntArray& triangles, short cageMode, std::vector<Matrix4d>& m);
    void tetMatrixC(const MPointArray& p, const MIntArray& triangles, std::vector<Matrix4d>& m, std::vector<Vector3d>& tetCenter);
	void arapHI(const std::vector<Matrix4d>& PI, const MIntArray& triangles);
	void arapG(const std::vector< Matrix4d>& At, const std::vector<Matrix4d>& PI,
               const MIntArray& triangles, const std::vector<Matrix4d>& aff, MatrixXd& G);
    MObject     initCageMesh;
    MPointArray initCagePoints;
    MPointArray pts;        //  target mesh points
	MIntArray triangles;    // cage triangles
	int numPrb,numPts,numTet,numConstraint;
    short cageMode;
    short constraintMode;
    double constraintWeight;
    std::vector<Matrix4d> initMatrix;    // initial matrix for cage
    std::vector< std::vector<double> > w;  // weight for the mesh face tetrahedra
	std::vector<Vector3d> prevNs;         // for continuous log
	std::vector<double> prevThetas;        //  for continuous log
    // ARAP
	std::vector<Matrix4d> PI;      // inverse of tet matrix
	std::vector<Vector3d> tetCenter;     // barycenter of mesh tetrahedra
    std::vector<int> constraintTet;      // constrained tet (nearest to some cage face or point)
    std::vector<RowVector4d> constraintVector;   // initial position of constraint points
    SimplicialLDLT<SpMat> solver;
//    SparseLU<SpMat> solver;    // ARAP solver
    SpMat F;                // ARAP constraint matrix
	MIntArray meshTriangles;   // target mesh triangles
};

