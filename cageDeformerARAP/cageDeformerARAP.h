/**
 * @file cageDeformerARAP.h
 * @brief Cage Deformer with ARAP mod plugin for Maya
 * @section LICENSE The MIT License
 * @section  requirements:  Eigen library, Maya, Matrixlib
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

#include "affinelib.h"

typedef Eigen::SparseMatrix<float> SpMat;
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
	void tetMatrix(const MPointArray& p, const MIntArray& triangles, short cageMode, std::vector<Matrix4f>& m);
    void tetMatrixC(const MPointArray& p, const MIntArray& triangles, std::vector<Matrix4f>& m, std::vector<Vector3f>& tetCenter);
	void arapHI(const std::vector<Matrix4f>& PI, const MIntArray& triangles);
	void arapG(const std::vector< Matrix4f>& At, const std::vector<Matrix4f>& PI,
               const MIntArray& triangles, const std::vector<Matrix4f>& aff, MatrixXf& G);
    MObject     initCageMesh;
    MPointArray initCagePoints;
    MPointArray pts;        //  target mesh points
	MIntArray triangles;    // cage triangles
	int numPrb,numPts,numTet,numConstraint;
    short cageMode;
    short constraintMode;
    float constraintWeight;
    std::vector<Matrix4f> initMatrix;    // initial matrix for cage
    std::vector< std::vector<float> > w;  // weight for the mesh face tetrahedra
	std::vector<Vector3f> prevNs;         // for continuous log
	std::vector<float> prevThetas;        //  for continuous log
    // ARAP
	std::vector<Matrix4f> PI;      // inverse of tet matrix
	std::vector<Vector3f> tetCenter;     // barycenter of mesh tetrahedra
    std::vector<int> constraintTet;      // constrained tet (nearest to some cage face or point)
    std::vector<RowVector4f> constraintVector;   // initial position of constraint points
    SparseLU<SpMat> solver;    // ARAP solver
    SpMat F;                // ARAP constraint matrix
	MIntArray meshTriangles;   // target mesh triangles
};

