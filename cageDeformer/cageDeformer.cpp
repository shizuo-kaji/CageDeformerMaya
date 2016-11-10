/**
 * @file cageDeformer.cpp
 * @brief Cage Deformer plugin for Maya
 * @section LICENSE The MIT License
 * @section requirements:  Eigen 3:  http://eigen.tuxfamily.org/
 * @section Autodesk Maya: http://www.autodesk.com/products/autodesk-maya/overview
 * @section (included) AffineLib: https://github.com/shizuo-kaji/AffineLib
 * @version 0.20
 * @date  19/Sep/2016
 * @author Shizuo KAJI
 */

#include "StdAfx.h"
#include "cageDeformer.h"



using namespace Eigen;
using namespace AffineLib;
using namespace Tetrise;

MTypeId CageDeformerNode::id( 0x00000200 );
MString CageDeformerNode::nodeName("cageDeformer");
MObject CageDeformerNode::aCageMesh;
MObject CageDeformerNode::aCageMode;
MObject CageDeformerNode::aSymmetricFace;
MObject CageDeformerNode::aBlendMode;
MObject CageDeformerNode::aRotationConsistency;
MObject CageDeformerNode::aFrechetSum;
MObject CageDeformerNode::aNormExponent;
MObject CageDeformerNode::aWeightMode;
MObject CageDeformerNode::aNormaliseTet;
MObject CageDeformerNode::aReconstructCage;
MObject CageDeformerNode::aEffectRadius;
MObject CageDeformerNode::aNormaliseWeight;
MObject CageDeformerNode::aAreaWeighted;
MObject CageDeformerNode::aNeighbourWeighting;
MObject CageDeformerNode::aPositiveWeight;


void* CageDeformerNode::creator() { return new CageDeformerNode; }
 
MStatus CageDeformerNode::deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex )
{
    MStatus status;
    MThreadUtils::syncNumOpenMPThreads();    // for OpenMP
    // load cage mesh
    MObject oCageMesh = data.inputValue( aCageMesh ).asMesh();
    short blendMode = data.inputValue(aBlendMode).asShort();
    short weightMode = data.inputValue( aWeightMode ).asShort();
	B.rotationConsistency = data.inputValue( aRotationConsistency ).asBool();
	bool frechetSum = data.inputValue( aFrechetSum ).asBool();
    bool symmetricFace = data.inputValue( aSymmetricFace ).asBool();
    bool normaliseTet = data.inputValue( aNormaliseTet ).asBool();
    bool areaWeighted = data.inputValue( aAreaWeighted ).asBool();
    double effectRadius = data.inputValue( aEffectRadius ).asDouble();
    double normExponent = data.inputValue( aNormExponent ).asDouble();
    if ( oCageMesh.isNull() || blendMode == BM_OFF)
        return MS::kSuccess;
    short cageMode = data.inputValue(aCageMode).asShort();
    MFnMesh fnCageMesh( oCageMesh, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MPointArray cagePoints;
    fnCageMesh.getPoints( cagePoints,  MSpace::kWorld );
    int numCagePts = cagePoints.length();
    std::vector<Vector3d> cagePts(numCagePts);
    for(int i=0;i<numCagePts;i++){
        cagePts[i] << cagePoints[i].x, cagePoints[i].y, cagePoints[i].z;
    }
    // save initial cage state
    bool cageChanged = false;
    if (initCagePts.size() != numCagePts){
        initCageMesh = oCageMesh;
        initCagePts = cagePts;
        cageChanged = true;
    }
    
    // load target mesh pts
    MPointArray Mpts;
    itGeo.allPositions(Mpts);
    int numPts = Mpts.length();
    std::vector<Vector3d> pts(numPts);
    for(int j=0; j<numPts; j++ ){
        Mpts[j] *= localToWorldMatrix;
        pts[j] << Mpts[j].x, Mpts[j].y, Mpts[j].z;
    }
    
    // when cage mode is changed
    if(!data.isClean(aReconstructCage) || cageChanged){
        std::vector<int> cageFaceCount(0);
        
        if( cageMode == TM_FACE || cageMode == TM_EDGE || cageMode == TM_VERTEX || cageMode == TM_VFACE){
            // set cage tetrahedra
            makeFaceList( initCageMesh, cageFaceList, cageFaceCount, symmetricFace);
            std::vector<Matrix4d> initCageMatrix;
            makeVertexList(initCageMesh, cageVertexList);
            makeEdgeList(cageFaceList, cageEdgeList);
            makeTetList(cageMode, numCagePts, cageFaceList, cageEdgeList, cageVertexList, cageTetList);
            makeTetMatrix(cageMode, initCagePts, cageTetList, cageFaceList, cageEdgeList, cageVertexList, initCageMatrix, cageTetWeight, normaliseTet);
            makeTetCenterList(cageMode, initCagePts, cageTetList, cageTetCenter);
            numCageTet = (int) cageTetList.size()/4;
            cageMatrixI.resize(numCageTet);
            for(int i=0;i<numCageTet;i++){
                cageMatrixI[i] = initCageMatrix[i].inverse().eval();
            }
        }else{
            numCageTet = numCagePts;
        }

        if(!areaWeighted){
            cageTetWeight.assign(numCageTet,1.0);
        }
        // weight computation
        D.setNum(numCageTet, numPts, 0);
        D.computeCageDistPts(cageMode, pts, initCagePts, cageTetList);
        w.resize(numPts);
        for(int j=0;j<numPts;j++){
            w[j].resize(numCageTet);
        }
        // first, compute the distance between cage transforms and mesh points
        switch (cageMode){
            case TM_FACE:
            {
                for(int i=0;i<numCageTet;i++){
                    cageTetWeight[i] /= cageFaceCount[i];
                }
                break;
            }
            case TM_EDGE:
            {
                for(int i=0;i<cageEdgeList.size();i++){
                    for(int k=0;k<2;k++){
                        cageTetWeight[2*i+k] /= cageFaceCount[cageEdgeList[i].faces[k]];
                    }
                }
                break;
            }
            case TM_VERTEX:
            case TM_VFACE:
            {
                int cur=0;
                for(int i=0;i<numCagePts;i++){
                    for(int k=0;k<cageVertexList[i].connectedTriangles.size()/2;k++){
                        cageTetWeight[cur] /= cageVertexList[i].connectedTriangles.size()/2;
                        cur++;
                    }
                }
                break;
            }
            default:
            {
                cageTetWeight.assign(numCageTet,1.0);
                break;
            }
        }

        //
        if(cageMode == CM_MVC){
            D.MVC(pts, initCagePts, cageFaceList, w);
        }else if(weightMode == WM_INV_DISTANCE){
            for(int j=0; j<numPts; j++ ){
                for(int i=0;i<numCageTet;i++){
                    w[j][i] = cageTetWeight[i] / pow(D.distPts[i][j],normExponent);
                }
            }
        }else if(weightMode == WM_CUTOFF_DISTANCE){
            double delta;  // avoid under determined system for MLS
            delta = (cageMode == CM_MLS_AFF || cageMode == CM_MLS_SIM || cageMode == CM_MLS_RIGID) ? EPSILON : 0;
            for(int j=0; j<numPts; j++ ){
                for( int i=0; i<numCageTet; i++){
                    w[j][i] = (D.distPts[i][j] > effectRadius)
                    ? delta : cageTetWeight[i]*pow((effectRadius-D.distPts[i][j])/effectRadius,normExponent);
                }
            }
        }else if(weightMode & WM_HARMONIC){
            Laplacian harmonicWeighting;
            makeFaceTet(data, input, inputGeom, mIndex, pts, harmonicWeighting.tetList, harmonicWeighting.tetMatrix, harmonicWeighting.tetWeight);
            harmonicWeighting.numTet = (int)harmonicWeighting.tetList.size()/4;
            std::vector<T> weightConstraint(numCageTet);
            // the vertex closest to the probe is given probeWeight
            D.findClosestPts();
            for(int i=0;i<numCageTet;i++){
                weightConstraint[i]=T(i,D.closestPts[i],cageTetWeight[i]);
            }
            // vertices within effectRadius are given probeWeight
            if( data.inputValue( aNeighbourWeighting ).asBool() ){
                for(int i=0;i<numCageTet;i++){
                    for(int j=0;j<numPts;j++){
                        if(D.distPts[i][j]<effectRadius){
                            weightConstraint.push_back(T(i,j,cageTetWeight[i]));
                        }
                    }
                }
            }
            // set boundary condition for weight computation
            int numConstraint=weightConstraint.size();
            harmonicWeighting.constraintWeight.resize(numConstraint);
            harmonicWeighting.constraintVal.resize(numConstraint,numCageTet);
            harmonicWeighting.constraintVal.setZero();
            for(int i=0;i<numConstraint;i++){
                harmonicWeighting.constraintVal(i,weightConstraint[i].row())=weightConstraint[i].value();
                harmonicWeighting.constraintWeight[i] = std::make_pair(weightConstraint[i].col(), weightConstraint[i].value());
            }
            // clear tetWeight
            if(!areaWeighted){
                harmonicWeighting.tetWeight.assign(harmonicWeighting.numTet,1.0);
            }
            // solve the laplace equation
            int isError;
            if( weightMode == WM_HARMONIC_ARAP ){
                harmonicWeighting.computeTetMatrixInverse();
                harmonicWeighting.dim = numPts + harmonicWeighting.numTet;
                isError = harmonicWeighting.ARAPprecompute();
            }else if(weightMode == WM_HARMONIC_COTAN){
                harmonicWeighting.dim = numPts;
                isError = harmonicWeighting.cotanPrecompute();
            }
            if(isError>0) return MS::kFailure;
            harmonicWeighting.harmonicSolve();
            for(int i=0;i<numCageTet;i++){
                for(int j=0;j<numPts;j++){
                    w[j][i] = harmonicWeighting.Sol(j,i);
                }
            }
        }
        // normalise weights
        if (data.inputValue( aPositiveWeight ).asBool()){
            for(int j=0;j<numPts;j++){
                for (int i = 0; i < numCageTet; i++){
                    w[j][i] = max(w[j][i], 0.0);
                }
            }
        }
        bool normaliseWeight = data.inputValue( aNormaliseWeight ).asBool();
        for(int j=0;j<numPts;j++){
            double sum = std::accumulate(w[j].begin(), w[j].end(), 0.0);
            if ( (sum > 1 || normaliseWeight || cageMode & CM_MLS) && sum>0){
                for (int i = 0; i < numCageTet; i++){
                    w[j][i] /= sum;
                }
            }
        }
        data.setClean(aReconstructCage);
    }
    
    // transform vertices
	if(cageMode == TM_FACE || cageMode == TM_EDGE || cageMode == TM_VERTEX || cageMode == TM_VFACE){
        //  find affine transformations for tetrahedra
        B.setNum(numCageTet);
        std::vector<double> dummy_weight;
        makeTetMatrix(cageMode, cagePts, cageTetList, cageFaceList, cageEdgeList,
                      cageVertexList, cageMatrix, dummy_weight, normaliseTet);
        for(int i=0; i<numCageTet; i++)
            B.Aff[i]=cageMatrixI[i]*cageMatrix[i];
        B.parametrise(blendMode);
        
        #pragma omp parallel for
        for(int j=0; j<numPts; j++ ){
            Matrix4d mat;
            // blend matrix
            if(blendMode == BM_SRL){
                Matrix3d RR,SS=expSym(blendMat(B.logS, w[j]));
                Vector3d l=blendMat(B.L, w[j]);
                RR = frechetSum ? frechetSO(B.R, w[j]) : expSO(blendMat(B.logR, w[j]));
                mat = pad(SS*RR, l);
            }else if(blendMode == BM_SSE){
                Matrix4d RR;
                Matrix3d SS=expSym(blendMat(B.logS, w[j]));
                RR = expSE(blendMat(B.logSE, w[j]));
                mat = pad(SS,Vector3d::Zero()) * RR;
            }else if(blendMode == BM_LOG3){
                Matrix3d RR=blendMat(B.logGL, w[j]).exp();
                Vector3d l=blendMat(B.L, w[j]);
                mat = pad(RR, l);
            }else if(blendMode == BM_LOG4){
                mat=blendMat(B.logAff, w[j]).exp();
            }else if(blendMode == BM_SQL){
                Vector4d q=blendQuat(B.quat,w[j]);
                Vector3d l=blendMat(B.L, w[j]);
                Matrix3d SS=blendMatLin(B.S,w[j]);
                Quaternion<double> Q(q);
                Matrix3d RR = Q.matrix().transpose();
                mat = pad(SS*RR, l);
            }else if(blendMode == BM_AFF){
                mat = blendMatLin(B.Aff,w[j]);
            }
            // apply matrix
            RowVector4d p = pad(pts[j]) * mat;
            Mpts[j].x = p[0];
            Mpts[j].y = p[1];
            Mpts[j].z = p[2];
            Mpts[j] *= localToWorldMatrix.inverse();
        }
	}else if(cageMode == CM_MVC){
		for(int j=0; j<numPts; j++ ){
			Mpts[j] = MPoint::origin;
			for(int i=0;i<numCagePts;i++){
				Mpts[j] +=  cagePoints[i] * w[j][i];
            }
            Mpts[j] *= localToWorldMatrix.inverse();
		}
    }else if(cageMode & CM_MLS){
        for(int j=0; j<numPts; j++ ){
            // barycentre of the original and the current cage points
            Vector3d icenter = Vector3d::Zero();
            Vector3d center = Vector3d::Zero();
            for(int i=0;i<numCagePts;i++){
                icenter += w[j][i] * initCagePts[i];
                center += w[j][i] * cagePts[i];
            }
            // determine the closest affine matrix
            MatrixXd p(3,numCagePts), q(3,numCagePts);   // relative coordinates of the cage
            Matrix3d M,S;
            for(int i=0;i<numCagePts;i++){
                p.col(i) = (w[j][i]) * (initCagePts[i]-icenter);
                q.col(i) = (w[j][i]) * (cagePts[i]-center);
            }
            JacobiSVD<Matrix3d> svd(p*q.transpose(), ComputeFullU | ComputeFullV);
            M = (p*p.transpose()).inverse() * (p*q.transpose());
            Matrix3d sgn = Matrix3d::Identity();
            if(M.determinant()<0){
                sgn(2,2) = -1;
            }
            if(cageMode == CM_MLS_SIM){
                M = svd.matrixU() * sgn * svd.matrixV().transpose();
                M *= svd.singularValues().sum()/(p*p.transpose()).trace();
            }else if(cageMode == CM_MLS_RIGID){
                M = svd.matrixU() * sgn * svd.matrixV().transpose();
            }
            RowVector4d v = pad(pts[j]) * pad(M, center-icenter);
            Mpts[j].x = v[0];
            Mpts[j].y = v[1];
            Mpts[j].z = v[2];
            Mpts[j] *= localToWorldMatrix.inverse();
        }
    }
    itGeo.setAllPositions(Mpts);
 
    return MS::kSuccess;
}


// maya plugin initialization
MStatus CageDeformerNode::initialize(){
    MFnTypedAttribute tAttr;
    MFnNumericAttribute nAttr;
    MFnEnumAttribute eAttr;
 
    // this attr will be dirtied when weight recomputation is needed
    aReconstructCage = nAttr.create( "reconstructCage", "reconstructCage", MFnNumericData::kBoolean, true );
    nAttr.setHidden(true);
    nAttr.setStorable(false);
    nAttr.setKeyable(false);
    addAttribute( aReconstructCage );

    aCageMesh = tAttr.create( "cageMesh", "cm", MFnData::kMesh );
    addAttribute( aCageMesh );
    attributeAffects( aCageMesh, outputGeom );
 
    aCageMode = eAttr.create( "cageMode", "cgm", TM_VERTEX );
    eAttr.addField( "face", TM_FACE );
    eAttr.addField( "edge", TM_EDGE );
    eAttr.addField( "vertex", TM_VERTEX );
    eAttr.addField( "vface", TM_VFACE );
    eAttr.addField( "MVC", CM_MVC );
    eAttr.addField( "MLS-AFF", CM_MLS_AFF );
    eAttr.addField( "MLS-SIM", CM_MLS_SIM );
    eAttr.addField( "MLS-RIGID", CM_MLS_RIGID );
    addAttribute( aCageMode );
    attributeAffects( aCageMode, outputGeom );
    attributeAffects( aCageMode, aReconstructCage );
    
    aSymmetricFace = nAttr.create( "symmetricFace", "sf", MFnNumericData::kBoolean, true );
    nAttr.setStorable(true);
    addAttribute( aSymmetricFace );
    attributeAffects( aSymmetricFace, outputGeom );
    attributeAffects( aSymmetricFace, aReconstructCage );

    aBlendMode = eAttr.create( "blendMode", "bm", BM_SRL );
    eAttr.addField( "expSO+expSym", BM_SRL );
    eAttr.addField( "expSE+expSym", BM_SSE );
    eAttr.addField( "logmatrix3", BM_LOG3 );
    eAttr.addField( "logmatrix4", BM_LOG4 );
    eAttr.addField( "quat+linear", BM_SQL );
    eAttr.addField( "linear", BM_AFF );
    eAttr.addField( "off", BM_OFF );
    eAttr.setStorable(true);
    addAttribute( aBlendMode );
    attributeAffects( aBlendMode, outputGeom );

	aRotationConsistency = nAttr.create( "rotationConsistency", "rc", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aRotationConsistency );
    attributeAffects( aRotationConsistency, outputGeom );
    
	aNormaliseTet = nAttr.create( "normaliseTet", "nr", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aNormaliseTet );
    attributeAffects( aNormaliseTet, outputGeom );
    attributeAffects( aNormaliseTet, aReconstructCage );

	aFrechetSum = nAttr.create( "frechetSum", "fs", MFnNumericData::kBoolean, 0 );
    nAttr.setStorable(true);
    addAttribute( aFrechetSum );
    attributeAffects( aFrechetSum, outputGeom );

    aNormaliseWeight = nAttr.create( "normaliseWeight", "nw", MFnNumericData::kBoolean, true );
    nAttr.setStorable(true);
    addAttribute( aNormaliseWeight );
    attributeAffects( aNormaliseWeight, outputGeom );
    attributeAffects( aNormaliseWeight, aReconstructCage );

    aAreaWeighted = nAttr.create( "areaWeighted", "aw", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aAreaWeighted );
    attributeAffects( aAreaWeighted, outputGeom );
    attributeAffects( aAreaWeighted, aReconstructCage );

    aPositiveWeight = nAttr.create( "positiveWeight", "posw", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aPositiveWeight );
    attributeAffects( aPositiveWeight, outputGeom );
    attributeAffects( aPositiveWeight, aReconstructCage  );

    aWeightMode = eAttr.create( "weightMode", "wtm", WM_INV_DISTANCE );
    eAttr.addField( "inverse", WM_INV_DISTANCE );
    eAttr.addField( "cut-off", WM_CUTOFF_DISTANCE );
//    eAttr.addField( "harmonic-arap", WM_HARMONIC );
    eAttr.addField( "harmonic-cotan", WM_HARMONIC_COTAN );
    eAttr.setStorable(true);
    addAttribute( aWeightMode );
    attributeAffects( aWeightMode, outputGeom );
    attributeAffects( aWeightMode, aReconstructCage );

    aEffectRadius = nAttr.create("effectRadius", "er", MFnNumericData::kDouble, 8.0);
    nAttr.setMin( EPSILON );
    nAttr.setStorable(true);
    addAttribute( aEffectRadius );
    attributeAffects( aEffectRadius, outputGeom );
    attributeAffects( aEffectRadius, aReconstructCage  );
    
    aNormExponent = nAttr.create("normExponent", "ne", MFnNumericData::kDouble, 2.0);
    nAttr.setStorable(true);
	addAttribute( aNormExponent );
	attributeAffects( aNormExponent, outputGeom );
	attributeAffects( aNormExponent, aReconstructCage );
    
    aNeighbourWeighting = nAttr.create( "neighbourWeighting", "nghbrw", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aNeighbourWeighting );
    attributeAffects( aNeighbourWeighting, outputGeom );
    attributeAffects( aNeighbourWeighting, aReconstructCage );
    
    
    return MS::kSuccess;
}
 
MStatus initializePlugin( MObject obj )
{
    MStatus status;
    MFnPlugin plugin( obj, "Shizuo KAJI", "0.1", "Any");
 
    status = plugin.registerNode( CageDeformerNode::nodeName, CageDeformerNode::id, CageDeformerNode::creator, CageDeformerNode::initialize, MPxNode::kDeformerNode );
    CHECK_MSTATUS_AND_RETURN_IT( status );
 
    return status;
}
 
MStatus uninitializePlugin( MObject obj )
{
    MStatus   status;
    MFnPlugin plugin( obj );
 
    status = plugin.deregisterNode( CageDeformerNode::id );
    CHECK_MSTATUS_AND_RETURN_IT( status );
 
    return status;
}
