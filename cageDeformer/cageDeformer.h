#pragma once

#pragma comment(linker, "/export:initializePlugin /export:uninitializePlugin")

#include <maya/MFnPlugin.h>

#include <numeric>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Geometry>

#include "../affinelib.h"
#include "../tetrise.h"
#include "../MeshMaya.h"
#include "../laplacian.h"
#include "../blendAff.h"
#include "../distance.h"
#include "../deformerConst.h"

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
    static MObject      aNeighbourWeighting;
    static MObject      aPositiveWeight;
    
private:
    MObject     initCageMesh;
    BlendAff    B;
    Distance    D;
    int         numCageTet;
    std::vector<Matrix4d> cageMatrixI;   // inverses of initial cage matrix
    std::vector< std::vector<double> > w;   // weight for the target mesh points
    std::vector<int> cageFaceList;
    std::vector<Vector3d> cageTetCenter;
    std::vector< int > cageTetList;
    std::vector< edge > cageEdgeList;
    std::vector<vertex> cageVertexList;
    std::vector<Vector3d> initCagePts;
    std::vector<double> cageTetWeight;
    std::vector<Matrix4d> cageMatrix;
};

