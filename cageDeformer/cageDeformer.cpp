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

// cage mode
#define CM_MVC 20
#define CM_MLS 30



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

void* CageDeformerNode::creator() { return new CageDeformerNode; }
 
MStatus CageDeformerNode::deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex )
{
    MStatus status;
    MThreadUtils::syncNumOpenMPThreads();    // for OpenMP
    // load cage mesh
    MObject oCageMesh = data.inputValue( aCageMesh ).asMesh();
    short blendMode = data.inputValue(aBlendMode).asShort();
    short weightMode = data.inputValue( aWeightMode ).asShort();
	bool rotationCosistency = data.inputValue( aRotationConsistency ).asBool();
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
    int numCageTet = (int) cageTetList.size()/4;

    
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
            cageTetWeight.clear();
            cageTetWeight.resize(numCageTet,1.0);
        }
        // weight computation
        w.resize(numPts);
        std::vector< std::vector<double> > distPts(numPts);
        for(int j=0;j<numPts;j++){
            distPts[j].resize(numCageTet);
            w[j].resize(numCageTet);
        }
        // first, compute the distance between cage transforms and mesh points
        switch (cageMode){
            case TM_FACE:
            {
                for(int i=0;i<numCageTet;i++){
                    cageTetWeight[i] /= cageFaceCount[i];
                }
                for(int j=0;j<numPts;j++){
                    for(int i=0;i<numCageTet;i++){
                        Vector3d a=initCagePts[cageTetList[4*i]];
                        Vector3d b=initCagePts[cageTetList[4*i+1]];
                        Vector3d c=initCagePts[cageTetList[4*i+2]];
                        distPts[j][i] = distPtTri(pts[j], a,b,c);
                    }
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
                for(int j=0;j<numPts;j++){
                    for(int i=0;i<numCageTet;i++){
                        Vector3d a=initCagePts[cageTetList[4*i]];
                        Vector3d b=initCagePts[cageTetList[4*i+1]];
                        distPts[j][i] = distPtLin(pts[j], a,b);
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
                for(int j=0;j<numPts;j++){
                    for(int i=0;i<numCageTet;i++){
                        distPts[j][i] = (pts[j]-cageTetCenter[i]).norm();
                    }
                }
                break;
            }
            case CM_MLS:
            {
                for(int j=0;j<numPts;j++){
                    for(int i=0;i<numCagePts;i++){
                        distPts[j][i] = (pts[j]-initCagePts[i]).norm();
                    }
                }
                break;
            }
        }

        //
        if(cageMode == CM_MVC){
            MVC(pts, initCagePts, cageFaceList, w);
        }else if(weightMode == WM_INV_DISTANCE){
            for(int j=0; j<numPts; j++ ){
                for(int i=0;i<numCageTet;i++){
                    w[j][i] = cageTetWeight[i] / pow(distPts[j][i],normExponent);
                }
            }
        }else if(weightMode == WM_CUTOFF_DISTANCE){
            double delta;  // avoid under determined system for MLS
            delta = (cageMode == CM_MLS) ? EPSILON : 0;
            for(int j=0; j<numPts; j++ ){
                for( int i=0; i<numCageTet; i++){
                    w[j][i] = (distPts[j][i] > effectRadius)
                    ? delta : cageTetWeight[i]*pow((effectRadius-distPts[j][i])/effectRadius,normExponent);
                }
            }
        }else if(weightMode == WM_HARMONIC || weightMode == WM_HARMONIC_NEIGHBOUR){
            std::vector<int> fList,tList;
            std::vector< std::vector<double> > ptsWeight(numCageTet), w_tet(numCageTet);
            std::vector<Matrix4d> P;
            std::vector<double> fWeight;
            int d=makeFaceTet(data, input, inputGeom, mIndex, pts, fList, tList, P, fWeight);
            std::vector< std::map<int,double> > weightConstraint(numCageTet);
            std::vector<double> weightConstraintValue(0);
            for(int i=0;i<numCageTet;i++){
                weightConstraint[i].clear();
            }
            if( weightMode == WM_HARMONIC_NEIGHBOUR ){
                for(int i=0;i<numCageTet;i++){
                    for(int j=0;j<numPts;j++){
                        if(distPts[i][j]<effectRadius){
                            weightConstraint[i][j] = 1;
                            weightConstraintValue.push_back(cageTetWeight[i]);
                        }
                    }
                }
            }else if( weightMode == WM_HARMONIC){
                std::vector<int> closestPts(numCageTet);
                for(int i=0;i<numCageTet;i++){
                    weightConstraint[i].clear();
                    closestPts[i] = 0;
                    double min_d = HUGE_VAL;
                    for(int j=0;j<numPts;j++){
                        if( distPts[j][i] < min_d){
                            min_d = distPts[j][i];
                            closestPts[i] = j;
                        }
                    }
                }
                for(int i=0;i<numCageTet;i++){
                    weightConstraint[i][closestPts[i]] = 1;
                    weightConstraintValue.push_back(cageTetWeight[i]);
                }
            }
            if(!areaWeighted){
                fWeight.clear();
                fWeight.resize(P.size(),1.0);
            }
            int isError = harmonicWeight(numPts, P, tList, fWeight, weightConstraint, weightConstraintValue, ptsWeight);
            if(isError>0) return MS::kFailure;
            for(int i=0;i<numCageTet;i++){
                for(int j=0;j<numPts;j++){
                    w[j][i] = ptsWeight[i][j];
                }
            }
        }
        // normalise weights
        bool normaliseWeight = data.inputValue( aNormaliseWeight ).asBool();
        for(int j=0;j<numPts;j++){
            double sum = std::accumulate(w[j].begin(), w[j].end(), 0.0);
            if (sum > 1 || normaliseWeight || cageMode == CM_MLS){
                for (int i = 0; i < numCageTet; i++){
                    w[j][i] /= sum;
                }
            }
        }
        data.setClean(aReconstructCage);
    }
    
    // transform vertices
	if(cageMode == TM_FACE || cageMode == TM_EDGE || cageMode == TM_VERTEX || cageMode == TM_VFACE){
        // clear previous rotation
        if( ! rotationCosistency || numCageTet != logR.size() || numCageTet != logSE.size()){
            logSE.clear();
            logSE.resize(numCageTet, Matrix4d::Zero().eval());
            logR.clear();
            logR.resize(numCageTet, Matrix3d::Zero().eval());
        }
        //  find affine transformations for tetrahedra
        SE.resize(numCageTet);
        logAff.resize(numCageTet); Aff.resize(numCageTet);
        R.resize(numCageTet); logS.resize(numCageTet); S.resize(numCageTet); logGL.resize(numCageTet);
        L.resize(numCageTet); quat.resize(numCageTet);
        std::vector<double> dummy_weight;
        makeTetMatrix(cageMode, cagePts, cageTetList, cageFaceList, cageEdgeList,
                      cageVertexList, cageMatrix, dummy_weight, normaliseTet);
        for(int i=0; i<numCageTet; i++)
            Aff[i]=cageMatrixI[i]*cageMatrix[i];
        // parametrisation
        if(blendMode == BM_SRL || blendMode == BM_SSE || blendMode == BM_SQL){
            for(int i=0;i<numCageTet;i++){
                parametriseGL(Aff[i].block(0,0,3,3), logS[i] ,R[i]);
                L[i] = transPart(Aff[i]);
                if(blendMode == BM_SRL){
                    logR[i]=logSOc(R[i], logR[i]);
                }else if(blendMode == BM_SSE){
                    SE[i]=pad(R[i], L[i]);
                    logSE[i]=logSEc(SE[i], logSE[i]);
                }else if(blendMode == BM_SQL){
                    Quaternion<double> Q(R[i].transpose());
                    quat[i] << Q.x(), Q.y(), Q.z(), Q.w();
                    S[i]=expSym(logS[i]);
                }
            }
        }else if(blendMode == BM_LOG3){
            for(int i=0;i<numCageTet;i++){
                logGL[i] = Aff[i].block(0,0,3,3).log();
                L[i] = transPart(Aff[i]);
            }
        }else if(blendMode == BM_LOG4){
            for(int i=0;i<numCageTet;i++){
                logAff[i] = Aff[i].log();
            }
        }
        
        #pragma omp parallel for
        for(int j=0; j<numPts; j++ ){
            Matrix4d mat;
            // blend matrix
            if(blendMode == BM_SRL){
                Matrix3d RR,SS=expSym(blendMat(logS, w[j]));
                Vector3d l=blendMat(L, w[j]);
                if(frechetSum){
                    RR = frechetSO(R, w[j]);
                }else{
                    RR = expSO(blendMat(logR, w[j]));
                }
                mat = pad(SS*RR, l);
            }else if(blendMode == BM_SSE){
                Matrix4d RR;
                Matrix3d SS=expSym(blendMat(logS, w[j]));
                RR = expSE(blendMat(logSE, w[j]));
                mat = pad(SS,Vector3d::Zero()) * RR;
            }else if(blendMode == BM_LOG3){
                Matrix3d RR=blendMat(logGL, w[j]).exp();
                Vector3d l=blendMat(L, w[j]);
                mat = pad(RR, l);
            }else if(blendMode == BM_LOG4){
                mat=blendMat(logAff, w[j]).exp();
            }else if(blendMode == BM_SQL){
                Vector4d q=blendQuat(quat,w[j]);
                Vector3d l=blendMat(L, w[j]);
                Matrix3d SS=blendMatLin(S,w[j]);
                Quaternion<double> Q(q);
                Matrix3d RR = Q.matrix().transpose();
                mat = pad(SS*RR, l);
            }else if(blendMode == BM_AFF){
                mat = blendMatLin(Aff,w[j]);
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
    }else if(cageMode == CM_MLS){
        for(int j=0; j<numPts; j++ ){
            // barycentre of the original and the current cage points
            Vector3d icenter = Vector3d::Zero();
            Vector3d center = Vector3d::Zero();
            for(int i=0;i<numCagePts;i++){
                icenter += w[j][i] * initCagePts[i];
                center += w[j][i] * cagePts[i];
            }
            // determine the closest affine matrix
            std::vector<Vector3d> p(numCagePts), q(numCagePts);   // relative coordinates of the cage
            Matrix3d M,P,Q;
            P = Matrix3d::Zero();
            Q = Matrix3d::Zero();
            for(int i=0;i<numCagePts;i++){
                p[i] = initCagePts[i]-icenter;
                q[i] = cagePts[i]-center;
                P += w[j][i] * (p[i] * p[i].transpose());
                Q += w[j][i] * p[i] * q[i].transpose();
            }
            M = P.inverse() * Q;

            RowVector4d v = pad(pts[j]) * pad(M,center-icenter);
            Mpts[j].x = v[0];
            Mpts[j].y = v[1];
            Mpts[j].z = v[2];
            Mpts[j] *= localToWorldMatrix.inverse();
        }
    }
    itGeo.setAllPositions(Mpts);
 
    return MS::kSuccess;
}

// mean value coordinate
void CageDeformerNode::MVC(const std::vector<Vector3d>& pts, const std::vector<Vector3d>& cagePts,
                           const std::vector<int>& cageFaceList, std::vector< std::vector<double> >& w)
{
	int numPts=(int) pts.size();
	int numCagePts=(int) cagePts.size();
	int numFaces=(int) cageFaceList.size()/3;
    w.resize(numPts);
#pragma omp parallel for
	for(int j=0; j<numPts; j++ ){
        w[j].resize(numCagePts);
		std::vector<double> mu(numCagePts), a(3), b(3);
		std::vector<Vector3d> e(3),n(3);
		for(int i=0;i<numFaces;i++){
			for(int k=0;k<3;k++)
				e[k]=(pts[j]-cagePts[cageFaceList[3*i+k]]).normalized();
			for(int k=0;k<3;k++)
				n[k]=(e[(k+1)%3].cross(e[(k+2)%3])).normalized();
			for(int k=0;k<3;k++){
				a[k]=n[(k+1)%3].dot(n[(k+2)%3]);
				b[k]=acos(e[(k+1)%3].dot(e[(k+2)%3]));
			}
			for(int k=0;k<3;k++)
                mu[cageFaceList[3*i+k]] -= (b[k]+b[(k+2)%3]*a[(k+1)%3]+b[(k+1)%3]*a[(k+2)%3])/(2.0*e[k].dot(n[k]));
		}
		double smu=0.0;
		for(int i=0;i<numCagePts;i++){
			mu[i] /= (pts[j]-cagePts[i]).norm();
			smu += mu[i];
		}
		for(int i=0;i<numCagePts;i++)
			w[j][i] = mu[i]/smu;
	}
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
 
    aCageMode = eAttr.create( "cageMode", "cgm", TM_FACE );
    eAttr.addField( "face", TM_FACE );
    eAttr.addField( "edge", TM_EDGE );
    eAttr.addField( "vertex", TM_VERTEX );
    eAttr.addField( "vface", TM_VFACE );
    eAttr.addField( "MVC", CM_MVC );
    eAttr.addField( "MLS", CM_MLS );
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

    aWeightMode = eAttr.create( "weightMode", "wtm", WM_INV_DISTANCE );
    eAttr.addField( "inverse", WM_INV_DISTANCE );
    eAttr.addField( "cut-off", WM_CUTOFF_DISTANCE );
    eAttr.addField( "harmonic-closest", WM_HARMONIC );
    eAttr.addField( "harmonic-neighbour", WM_HARMONIC_NEIGHBOUR );
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
