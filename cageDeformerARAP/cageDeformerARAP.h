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

#include "../affinelib.h"
#include "../tetrise.h"
#include "../deformerConst.h"
#include "../MeshMaya.h"
#include "../laplacian.h"
#include "../blendAff.h"
#include "../distance.h"

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
    static MObject		aEffectRadius;
    static MObject      aConstraintRadius;
    static MObject      aIteration;
    static MObject      aNormaliseWeight;
    static MObject      aAreaWeighted;
    static MObject      aNeighbourWeighting;
    static MObject      aPositiveWeight;
    
private:
    MObject     initCageMesh;
    BlendAff    B;
    Laplacian mesh;
    Distance D;
    int numCageTet;
    short isError;
    std::vector< std::vector<double> > w;   // weight
    // cage
    std::vector<Matrix4d> cageMatrixI;   // inverses of initial cage matrix
    std::vector<int> cageFaceList;
    std::vector< int > cageTetList;
    std::vector< edge > cageEdgeList;
    std::vector<vertex> cageVertexList;
    std::vector<Vector3d> initCagePts;
    std::vector<double> cageTetWeight;
    // mesh
    std::vector< edge > edgeList;
    std::vector<vertex> vertexList;
    std::vector<int> faceList;
    std::vector<Vector3d> pts;
    std::vector<T> constraint;
    //  find affine transformations for tetrahedra
    std::vector<Matrix4d> cageMatrix;
    std::vector<Matrix4d> A,blendedSE;
    std::vector<Matrix3d> blendedR, blendedS;
    std::vector<Vector4d> blendedL;
    std::vector<Vector3d> new_pts;
    std::vector<Matrix4d> Q;  // temporary
    std::vector<double> dummy_weight;

};

