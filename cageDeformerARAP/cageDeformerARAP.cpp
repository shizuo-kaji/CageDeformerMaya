/**
 * @file cageDeformerARAP.cpp
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

#include "StdAfx.h"
#include "cageDeformerARAP.h"


// cage mode
#define CM_MVC 20
#define CM_MLS 30


using namespace Eigen;
using namespace AffineLib;
using namespace Tetrise;

MTypeId CageDeformerNode::id( 0x00000202 );
MString CageDeformerNode::nodeName("cageDeformerARAP");
MObject CageDeformerNode::aARAP;
MObject CageDeformerNode::aCageMesh;
MObject CageDeformerNode::aCageMode;
MObject CageDeformerNode::aBlendMode;
MObject CageDeformerNode::aTetMode;
MObject CageDeformerNode::aRotationConsistency;
MObject CageDeformerNode::aFrechetSum;
MObject CageDeformerNode::aConstraintMode;
MObject CageDeformerNode::aConstraintWeight;
MObject CageDeformerNode::aTransWeight;
MObject CageDeformerNode::aNormaliseTet;
MObject CageDeformerNode::aSymmetricFace;
MObject CageDeformerNode::aNormExponent;
MObject CageDeformerNode::aIteration;
MObject CageDeformerNode::aWeightMode;
MObject CageDeformerNode::aEffectRadius;
MObject CageDeformerNode::aConstraintRadius;
MObject CageDeformerNode::aNormaliseWeight;
MObject CageDeformerNode::aAreaWeighted;

double isDegenerate(MPoint a,MPoint b,MPoint c,MPoint d){
    /// check linear independency
    return ((b-a)^(c-a)) * (d-a);
}

void* CageDeformerNode::creator() { return new CageDeformerNode; }

MStatus CageDeformerNode::deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex )
{
    /// main
    MStatus status;
    MThreadUtils::syncNumOpenMPThreads();    // for OpenMP
    // load attributes
    MObject oCageMesh = data.inputValue( aCageMesh ).asMesh();
    short blendMode = data.inputValue(aBlendMode).asShort();
	bool rotationCosistency = data.inputValue( aRotationConsistency ).asBool();
	bool frechetSum = data.inputValue( aFrechetSum ).asBool();
    short tetMode = data.inputValue(aTetMode).asShort();
    short constraintMode = data.inputValue(aConstraintMode).asShort();
    short numIter = data.inputValue( aIteration ).asShort();
    double constraintWeight = data.inputValue( aConstraintWeight ).asDouble();
    double constraintRadius = data.inputValue( aConstraintRadius ).asDouble();
    double transWeight = data.inputValue( aTransWeight ).asDouble();
    short weightMode = data.inputValue( aWeightMode ).asShort();
    double effectRadius = data.inputValue( aEffectRadius ).asDouble();
    bool symmetricFace = data.inputValue( aSymmetricFace ).asBool();
    bool normaliseTet = data.inputValue( aNormaliseTet ).asBool();
    double normExponent = data.inputValue( aNormExponent ).asDouble();
    bool areaWeighted = data.inputValue( aAreaWeighted ).asBool();
    bool arap_flag = false;
    // load cage
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
    // if the cage topology is changed
    if (initCagePts.size() != numCagePts){
        initCageMesh = oCageMesh;
        initCagePts = cagePts;
        arap_flag = true;
    }
    int numCageTet = (int) cageTetList.size()/4;
    // load mesh points
    MPointArray Mpts;
    itGeo.allPositions(Mpts);
    int numPts=Mpts.length();
    int numTet=(int)tetList.size()/4;


    // when cage mode is changed
    if(!data.isClean(aARAP) || arap_flag){
        // prepare mesh tetrahedra
        pts.resize(numPts);
        for(int j=0; j<numPts; j++ ){
            Mpts[j] *= localToWorldMatrix;
            pts[j] << Mpts[j].x, Mpts[j].y, Mpts[j].z;
        }
        std::vector<Matrix4d> P;
        getMeshData(data, input, inputGeom, mIndex, tetMode, pts, tetList, faceList, edgeList, vertexList, P, tetWeight);
        dim = removeDegenerate(tetMode, numPts, tetList, faceList, edgeList, vertexList, P);
        makeTetMatrix(tetMode, pts, tetList, faceList, edgeList, vertexList, P, tetWeight);
        //makeTetCenterList(tetMode, pts, tetList, tetCenter);
        numTet = (int)tetList.size()/4;
        PI.resize(numTet);
        for(int i=0;i<numTet;i++){
            PI[i] = P[i].inverse().eval();
        }
        
        // prepare cage tetrahedra
        std::vector<int> cageFaceCount(0);
        if (cageMode != CM_MLS){
            // face list
            makeFaceList( initCageMesh, cageFaceList, cageFaceCount, symmetricFace);
            std::vector<Matrix4d> initCageMatrix;
            makeFaceList(initCageMesh, cageFaceList, cageFaceCount);
            makeVertexList(initCageMesh, cageVertexList);
            makeEdgeList(cageFaceList, cageEdgeList);
            makeTetList(cageMode, numCagePts, cageFaceList, cageEdgeList, cageVertexList, cageTetList);
            makeTetMatrix(cageMode, initCagePts, cageTetList, cageFaceList, cageEdgeList,
                              cageVertexList, initCageMatrix, cageTetWeight, normaliseTet);
            makeTetCenterList(cageMode, initCagePts, cageTetList, cageTetCenter);
            numCageTet = (int) cageTetList.size()/4;
            cageMatrixI.resize(numCageTet);
            for(int i=0;i<numCageTet;i++){
                cageMatrixI[i] = initCageMatrix[i].inverse();
            }
        }else{
            numCageTet = numCagePts;
        }
        
        if(!areaWeighted){
            cageTetWeight.clear();
            cageTetWeight.resize(numCageTet,1.0);
            tetWeight.clear();
            tetWeight.resize(numTet,1.0);
        }
        // set effect weight for cage tet
        if(cageMode == TM_FACE){
            for(int i=0;i<numCageTet;i++){
                cageTetWeight[i] /= cageFaceCount[i];
            }
        }else if(cageMode == TM_EDGE){
            for(int i=0;i<cageEdgeList.size();i++){
                for(int k=0;k<2;k++){
                    cageTetWeight[2*i+k] /= cageFaceCount[cageEdgeList[i].faces[k]];
                }
            }
        }else if(cageMode == TM_VERTEX || cageMode == TM_VFACE){
            int cur=0;
            for(int i=0;i<numCagePts;i++){
                for(int k=0;k<cageVertexList[i].connectedTriangles.size()/2;k++){
                    cageTetWeight[cur] /= cageVertexList[i].connectedTriangles.size()/2;
                    cur++;
                }
            }
        }
    
    
        // pt to cage distance computation
        std::vector< std::vector<double> > dist(numPts);
        for(int j=0;j<numPts;j++){
            dist[j].resize(numCageTet);
        }
        if(cageMode == TM_FACE){
            for(int j=0;j<numPts;j++){
                for(int i=0;i<numCageTet;i++){
                    Vector3d a=initCagePts[cageTetList[4*i]];
                    Vector3d b=initCagePts[cageTetList[4*i+1]];
                    Vector3d c=initCagePts[cageTetList[4*i+2]];
                    dist[j][i] = sqrt(distPtTri(pts[j], a,b,c));
                }
            }
        }else if(cageMode == TM_EDGE){
            for(int j=0;j<numPts;j++){
                for(int i=0;i<numCageTet;i++){
                    Vector3d a=initCagePts[cageTetList[4*i]];
                    Vector3d b=initCagePts[cageTetList[4*i+1]];
                    dist[j][i] = sqrt(distPtLin(pts[j], a,b));
                }
            }
        }else if(cageMode == TM_VERTEX || cageMode == TM_VFACE){
            for(int j=0;j<numPts;j++){
                for(int i=0;i<numCageTet;i++){
                    dist[j][i] = (pts[j]-cageTetCenter[i]).norm();
                }
            }
        }else if(cageMode == CM_MLS){
            for(int j=0;j<numPts;j++){
                for(int i=0;i<numCageTet;i++){
                    dist[j][i] = (pts[j]-initCagePts[i]).norm();
                }
            }
        }
        // find closest point on mesh from each transformation
        std::vector<int> closest(numCageTet);
        for(int i=0;i<numCageTet;i++){
            closest[i] = 0;
            double min_d = HUGE_VAL;
            for(int j=0;j<numPts;j++){
                if( dist[j][i] < min_d){
                    min_d = dist[j][i];
                    closest[i] = j;
                }
            }
        }
        
        
        // find constraint points
        constraint.resize(numCageTet);
        if( constraintMode == CONSTRAINT_NEIGHBOUR ){
            for(int i=0;i<numCageTet;i++){
                constraint[i].clear();
                for(int j=0;j<numPts;j++){
                    if(dist[j][i]<constraintRadius){
                        constraint[i][j] = constraintWeight * cageTetWeight[i] *
                        pow((constraintRadius-dist[j][i])/constraintRadius,normExponent);
                    }
                }
            }
        }else if( constraintMode == CONSTRAINT_CLOSEST){
            for(int i=0;i<numCageTet;i++){
                constraint[i].clear();
                constraint[i][closest[i]] = constraintWeight*cageTetWeight[i];
            }
        }
        
        // weight computation
        w.resize(numTet);
        if(weightMode == WM_INV_DISTANCE){
            for(int j=0;j<numTet;j++){
                w[j].resize(numCageTet);
                for( int i=0; i<numCageTet; i++){
                    w[j][i] = cageTetWeight[i]/pow(dist[tetList[4*j]][i],normExponent);
                }
            }
        }else if(weightMode == WM_CUTOFF_DISTANCE){
            double delta; // avoid under determined system for MLS
            delta = (cageMode == CM_MLS) ? EPSILON : 0;
            for(int j=0;j<numTet;j++){
                w[j].resize(numCageTet);
                for( int i=0; i<numCageTet; i++){
                    w[j][i] = (dist[tetList[4*j]][i] > effectRadius)
                    ? delta : cageTetWeight[i]*pow((effectRadius-dist[tetList[4*j]][i])/effectRadius,normExponent);
                }
            }
        }else if(weightMode == WM_HARMONIC || weightMode == WM_HARMONIC_NEIGHBOUR){
            std::vector<int> fList,tList;
            std::vector<double> fWeight;
            std::vector< std::vector<double> > ptsWeight(numCageTet), w_tet(numCageTet);
            std::vector<Matrix4d> P;
            int d=makeFaceTet(data, input, inputGeom, mIndex, pts, fList, tList, P, fWeight);
            std::vector< std::map<int,double> > weightConstraint(numCageTet);
            std::vector<double> weightConstraintValue(0);
            for(int i=0;i<numCageTet;i++){
                weightConstraint[i].clear();
            }
            if( weightMode == WM_HARMONIC_NEIGHBOUR ){
                for(int i=0;i<numCageTet;i++){
                    for(int j=0;j<numPts;j++){
                        if(dist[i][j]<effectRadius){
                            weightConstraint[i][j] = 1;
                            weightConstraintValue.push_back(cageTetWeight[i]);
                        }
                    }
                }
            }else if( weightMode == WM_HARMONIC){
                for(int i=0;i<numCageTet;i++){
                    weightConstraint[i][closest[i]] = 1;
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
                makeTetWeightList(tetMode, tetList, faceList, edgeList, vertexList, ptsWeight[i], w_tet[i]);
                for(int j=0;j<numTet;j++){
                    w[j][i] = w_tet[i][j];
                }
            }
        }
        // normalise weights
        bool normaliseWeight = data.inputValue( aNormaliseWeight ).asBool();
        for(int j=0;j<numTet;j++){
            double sum = std::accumulate(w[j].begin(), w[j].end(), 0.0);
            if (sum > 1 || normaliseWeight || cageMode == CM_MLS){
                for (int i = 0; i < numCageTet; i++){
                    w[j][i] /= sum;
                }
            }
        }
        
        // prepare ARAP solver
        isError = ARAPprecompute(PI, tetList, tetWeight, constraint, transWeight, dim, constraintMat, solver);
        status = data.setClean(aARAP);
    }
//////  END of precomputation ///////
    if(isError>0){
        return MS::kFailure;
    }
    

    
    
    // compute transformation
    if( ! rotationCosistency || numCageTet != logR.size() || numCageTet != logSE.size()){
        logSE.clear();
        logSE.resize(numCageTet, Matrix4d::Zero().eval());
        logR.clear();
        logR.resize(numCageTet, Matrix3d::Zero().eval());
    }
    //  find affine transformations for tetrahedra
    cageMatrix.resize(numCageTet); SE.resize(numCageTet);
    logAff.resize(numCageTet); Aff.resize(numCageTet);
    R.resize(numCageTet); logS.resize(numCageTet); S.resize(numCageTet); logGL.resize(numCageTet);
    L.resize(numCageTet); quat.resize(numCageTet); A.resize(numTet); blendedSE.resize(numTet);
    blendedR.resize(numTet); blendedS.resize(numTet); blendedL.resize(numTet);
    
    
    if(cageMode == CM_MLS){
        for(int j = 0; j < numTet; j++){
            // barycentre of the original and the current cage points
            Vector3d icenter = Vector3d::Zero();
            Vector3d center = Vector3d::Zero();
            for(int i=0;i<numCagePts;i++){
                icenter += w[j][i] * initCagePts[i];
                center += w[j][i] * cagePts[i];
            }
            // determine the closest affine matrix
            std::vector<Vector3d> p(numCagePts), q(numCagePts);   // relative coordinates of the cage
            Matrix3d M,PP,QQ;
            PP = Matrix3d::Zero();
            QQ = Matrix3d::Zero();
            for(int i=0;i<numCagePts;i++){
                p[i] = initCagePts[i]-icenter;
                q[i] = cagePts[i]-center;
                PP += w[j][i] * (p[i] * p[i].transpose());
                QQ += w[j][i] * p[i] * q[i].transpose();
            }
            M = PP.inverse() * QQ;
//            polarHigham(M, S[j], R[j]);
            A[j] = pad(M,center-icenter);
        }
        for(int i=0; i<numCageTet; i++){
            Aff[i]=pad(Matrix3d::Identity().eval(),cagePts[i]-initCagePts[i]);
        }
    }else{
        std::vector<double> dummy_weight;
        makeTetMatrix(cageMode, cagePts, cageTetList, cageFaceList, cageEdgeList,
                      cageVertexList, cageMatrix, dummy_weight, normaliseTet);
        for(int i=0; i<numCageTet; i++)
            Aff[i]=cageMatrixI[i]*cageMatrix[i];
        
        // compute parametrisation
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
        
        // compute blended matrices
        // prepare transform matrix for each simplex
#pragma omp parallel for
        for (int j = 0; j < numTet; j++){
            // blend matrix
            if (blendMode == BM_SRL){
                blendedS[j] = expSym(blendMat(logS, w[j]));
//                blendedS[j] = frechetSum ? frechetSym(S, w[j]) : expSym(blendMat(logS, w[j]));
                Vector3d l = blendMat(L, w[j]);
                blendedR[j] = frechetSum ? frechetSO(R, w[j]) : expSO(blendMat(logR, w[j]));
                A[j] = pad(blendedS[j]*blendedR[j], l);
            }
            else if (blendMode == BM_SSE){
                blendedS[j] = expSym(blendMat(logS, w[j]));
                blendedSE[j] = expSE(blendMat(logSE, w[j]));
                A[j] = pad(blendedS[j], Vector3d::Zero()) * blendedSE[j];
            }
            else if (blendMode == BM_LOG3){
                blendedR[j] = blendMat(logGL, w[j]).exp();
                Vector3d l = blendMat(L, w[j]);
                A[j] = pad(blendedR[j], l);
            }
            else if (blendMode == BM_LOG4){
                A[j] = blendMat(logAff, w[j]).exp();
            }
            else if (blendMode == BM_SQL){
                Vector4d q = blendQuat(quat, w[j]);
                Vector3d l = blendMat(L, w[j]);
                blendedS[j] = blendMatLin(S, w[j]);
                Quaternion<double> Q(q);
                blendedR[j] = Q.matrix().transpose();
                A[j] = pad(blendedS[j]*blendedR[j], l);
            }
            else if (blendMode == BM_AFF){
                A[j] = blendMatLin(Aff, w[j]);
            }
        }
    }
    
    // compute target vertices position
    MatrixXd Sol;
    
    // set constraint
    std::vector<Vector3d> constraintVector(0);
    constraintVector.reserve(numCageTet * numPts);
    RowVector4d cv;
    Vector3d cvv;
    for(int i=0;i<constraint.size();i++){
        std::map<int, double>::iterator iter;
        for(iter = constraint[i].begin(); iter != constraint[i].end(); iter++){
            cv = pad(pts[iter->first]) * Aff[i];
            cvv << cv[0], cv[1], cv[2];
            constraintVector.push_back(cvv);
        }
    }
    
    // iterate to determine vertices position
    for(int k=0;k<numIter;k++){
        // solve ARAP
        ARAPSolve(A, PI, tetList, tetWeight, constraintVector, transWeight, dim, constraintMat, solver, Sol);
        // set new vertices position
        new_pts.resize(numPts);
        for(int i=0;i<numPts;i++){
            new_pts[i][0]=Sol(i,0);
            new_pts[i][1]=Sol(i,1);
            new_pts[i][2]=Sol(i,2);
        }
        // if iteration continues
        if(k+1<numIter){
            std::vector<double> dummy_weight;
            makeTetMatrix(tetMode, new_pts, tetList, faceList, edgeList, vertexList, Q, dummy_weight);
            Matrix3d S,R,newS,newR;
            if(blendMode == BM_AFF || blendMode == BM_LOG4 || blendMode == BM_LOG3 || cageMode == CM_MLS){
                for(int i=0;i<numTet;i++){
                    polarHigham(A[i].block(0,0,3,3), blendedS[i], blendedR[i]);
                }
            }
#pragma omp parallel for
            for(int i=0;i<numTet;i++){
                polarHigham((PI[i]*Q[i]).block(0,0,3,3), newS, newR);
//                tetEnergy[i] = (newS-blendedS[i]).squaredNorm();
                A[i].block(0,0,3,3) = blendedS[i]*newR;
            }
        }
    }
    for(int i=0;i<numPts;i++){
        Mpts[i].x=Sol(i,0);
        Mpts[i].y=Sol(i,1);
        Mpts[i].z=Sol(i,2);
        Mpts[i] *= localToWorldMatrix.inverse();
    }
    itGeo.setAllPositions(Mpts);
    
    return MS::kSuccess;
}


// maya plugin initialization
MStatus CageDeformerNode::initialize()
{
    MFnTypedAttribute tAttr;
    MFnNumericAttribute nAttr;
    MFnEnumAttribute eAttr;
    
    // this attr will be dirtied when ARAP recomputation is needed
    aARAP = nAttr.create( "arap", "arap", MFnNumericData::kBoolean, true );
    nAttr.setStorable(false);
    nAttr.setKeyable(false);
    nAttr.setHidden(true);
    addAttribute( aARAP );
    
    aCageMesh = tAttr.create( "cageMesh", "cm", MFnData::kMesh );
    addAttribute( aCageMesh );
    attributeAffects( aCageMesh, outputGeom );
    
    aCageMode = eAttr.create( "cageMode", "cgm", TM_FACE );
    eAttr.addField( "face", TM_FACE );
    eAttr.addField( "edge", TM_EDGE );
    eAttr.addField( "vertex", TM_VERTEX );
    eAttr.addField( "vface", TM_VFACE );
//    eAttr.addField( "MVC", CM_MVC );
    eAttr.addField( "MLS", CM_MLS );
    addAttribute( aCageMode );
    attributeAffects( aCageMode, outputGeom );
    attributeAffects( aCageMode, aARAP);
    
    
    aConstraintMode = eAttr.create( "constraintMode", "ctm", CONSTRAINT_NEIGHBOUR );
    eAttr.addField( "neighbour",  CONSTRAINT_NEIGHBOUR);
    eAttr.addField( "closestPt", CONSTRAINT_CLOSEST );
    eAttr.setStorable(true);
    addAttribute( aConstraintMode );
    attributeAffects( aConstraintMode, outputGeom );
    attributeAffects( aConstraintMode, aARAP);
    
    aConstraintWeight = nAttr.create("constraintWeight", "cw", MFnNumericData::kDouble, 0.001);
    nAttr.setStorable(true);
	addAttribute( aConstraintWeight );
	attributeAffects( aConstraintWeight, outputGeom );
    attributeAffects( aConstraintWeight, aARAP);
    
    aConstraintRadius = nAttr.create("constraintRadius", "cr", MFnNumericData::kDouble, 2.0);
    nAttr.setStorable(true);
	addAttribute( aConstraintRadius );
	attributeAffects( aConstraintRadius, outputGeom );
    attributeAffects( aConstraintRadius, aARAP );
    
    aTransWeight = nAttr.create("translationWeight", "tw", MFnNumericData::kDouble, 0.0001);
    nAttr.setStorable(true);
	addAttribute( aTransWeight );
	attributeAffects( aTransWeight, outputGeom );
    attributeAffects( aTransWeight, aARAP );
    
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
    
    aNormaliseWeight = nAttr.create( "normaliseWeight", "nw", MFnNumericData::kBoolean, true );
    nAttr.setStorable(true);
    addAttribute( aNormaliseWeight );
    attributeAffects( aNormaliseWeight, outputGeom );
    attributeAffects( aNormaliseWeight, aARAP );
    
    aWeightMode = eAttr.create( "weightMode", "wtm", WM_INV_DISTANCE );
    eAttr.addField( "inverse", WM_INV_DISTANCE );
    eAttr.addField( "cut-off", WM_CUTOFF_DISTANCE );
    eAttr.addField( "harmonic-closest", WM_HARMONIC);
    eAttr.addField( "harmonic-neighbour", WM_HARMONIC_NEIGHBOUR);
    eAttr.setStorable(true);
    addAttribute( aWeightMode );
    attributeAffects( aWeightMode, outputGeom );
    attributeAffects( aWeightMode, aARAP );
    
    aTetMode = eAttr.create( "tetMode", "tm", TM_VERTEX);
    eAttr.addField( "face", TM_FACE );
    eAttr.addField( "edge", TM_EDGE );
    eAttr.addField( "vertex", TM_VERTEX );
    eAttr.addField( "vface", TM_VFACE );
    eAttr.setStorable(true);
    addAttribute( aTetMode );
    attributeAffects( aTetMode, outputGeom );
    attributeAffects( aTetMode, aARAP );

	aNormaliseTet = nAttr.create( "normaliseTet", "nr", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aNormaliseTet );
    attributeAffects( aNormaliseTet, outputGeom );
    attributeAffects( aNormaliseTet, aARAP );

    aSymmetricFace = nAttr.create( "symmetricFace", "sf", MFnNumericData::kBoolean, true );
    nAttr.setStorable(true);
    addAttribute( aSymmetricFace );
    attributeAffects( aSymmetricFace, outputGeom );
    attributeAffects( aSymmetricFace, aARAP );

    aAreaWeighted = nAttr.create( "areaWeighted", "aw", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aAreaWeighted );
    attributeAffects( aAreaWeighted, outputGeom );
    attributeAffects( aAreaWeighted, aARAP );
    
	aRotationConsistency = nAttr.create( "rotationConsistency", "rc", MFnNumericData::kBoolean, 0 );
    nAttr.setStorable(true);
    addAttribute( aRotationConsistency );
    attributeAffects( aRotationConsistency, outputGeom );
    
    aIteration = nAttr.create("iteration", "it", MFnNumericData::kShort, 1);
    nAttr.setStorable(true);
    addAttribute(aIteration);
    attributeAffects(aIteration, outputGeom);
    
    aFrechetSum = nAttr.create( "frechetSum", "fs", MFnNumericData::kBoolean, 0 );
    nAttr.setStorable(true);
    addAttribute( aFrechetSum );
    attributeAffects( aFrechetSum, outputGeom );
    
    aNormExponent = nAttr.create("normExponent", "ne", MFnNumericData::kDouble, 2.0);
    nAttr.setStorable(true);
	addAttribute( aNormExponent );
	attributeAffects( aNormExponent, outputGeom );
    attributeAffects( aNormExponent, aARAP );

    aEffectRadius = nAttr.create("effectRadius", "er", MFnNumericData::kDouble, 8.0);
    nAttr.setMin( EPSILON );
    nAttr.setStorable(true);
    addAttribute( aEffectRadius );
    attributeAffects( aEffectRadius, outputGeom );
    attributeAffects( aEffectRadius, aARAP );

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
