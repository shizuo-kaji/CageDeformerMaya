#pragma once

#pragma comment(lib, "OpenMaya.lib")
#pragma comment(lib, "OpenMayaRender.lib")
#pragma comment(lib, "OpenMayaAnim.lib")
#pragma comment(lib, "Foundation.lib")
#pragma comment(linker, "/export:initializePlugin /export:uninitializePlugin")

#include <maya/MFnPlugin.h>

#include <numeric>
#include <unsupported/Eigen/MatrixFunctions>

#include "affinelib.h"

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

private:
	void tetMatrix(const MPointArray& p, const MIntArray& triangles, short cageMode, std::vector<Matrix4f>& m);
	void MVC(const MPointArray& pts, const MPointArray& cagePoints, const MIntArray& triangles, std::vector< std::vector<float> >& w);
    MObject     initCageMesh;
    MPointArray initCagePoints;
	MIntArray triangles;    // cage triangles
	int numTet;             // number of tetrahedra in cage
    short cageMode;
    std::vector<Matrix4f> initInvMatrix;   // inverses of initial cage matrix
    std::vector< std::vector<float> > w;   // weight for the target mesh points
	std::vector<Vector3f> prevNs;          //  for continuous log
	std::vector<float> prevThetas;         //  for continuous log
};

