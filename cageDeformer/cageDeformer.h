#pragma once

#pragma comment(linker, "/export:initializePlugin /export:uninitializePlugin")

#include <maya/MFnPlugin.h>

#include <numeric>
#include <unsupported/Eigen/MatrixFunctions>

#include "affinelib.h"
#include "tetrise.h"
#include "MeshMaya.h"
#include "ARAP.h"
#include "deformerConst.h"

using namespace Eigen;

class CageDeformerNode : public MPxDeformerNode{
public:
    CageDeformerNode()  {};
    virtual MStatus deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex );
    static  void*   creator();
    static  MStatus initialize();
 
    static MTypeId      id;
    static MString      nodeName;
    static MObject      aCageMesh;
    static MObject      aCageMode;
    static MObject      aBlendMode;
    static MObject      aSymmetricFace;
	static MObject		aRotationConsistency;
	static MObject		aFrechetSum;
    static MObject      aNormExponent;
    static MObject      aNormaliseTet;
    static MObject      aWeightMode;
    static MObject		aEffectRadius;
    static MObject      aReconstructCage;
    static MObject      aNormaliseWeight;
    static MObject      aAreaWeighted;
    
private:
	void MVC(const std::vector<Vector3d>& pts, const std::vector<Vector3d>& cagePoints,
             const std::vector<int>& cageFaceList, std::vector< std::vector<double> >& w);
    MObject     initCageMesh;
    std::vector<Matrix4d> cageMatrixI;   // inverses of initial cage matrix
    std::vector< std::vector<double> > w;   // weight for the target mesh points
    std::vector<Matrix4d> logSE;   // for rotation consistency
    std::vector<Matrix3d> logR;   // for rotation consistency
    std::vector<int> cageFaceList;
    std::vector<Vector3d> cageTetCenter;
    std::vector< int > cageTetList;
    std::vector< edge > cageEdgeList;
    std::vector<vertex> cageVertexList;
    std::vector<Vector3d> initCagePts;
    std::vector<double> cageTetWeight;
    //  find affine transformations for tetrahedra
    std::vector<Matrix4d> cageMatrix, SE, logAff,Aff;
    std::vector<Matrix3d> R,logS,S,logGL;
    std::vector<Vector3d> L;
    std::vector<Vector4d> quat;
};

