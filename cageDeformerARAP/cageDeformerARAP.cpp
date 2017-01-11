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
MObject CageDeformerNode::aNeighbourWeighting;
MObject CageDeformerNode::aPositiveWeight;


// check if points are coplanar
double isDegenerate(MPoint a,MPoint b,MPoint c,MPoint d){
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
	B.rotationConsistency = data.inputValue( aRotationConsistency ).asBool();
	bool frechetSum = data.inputValue( aFrechetSum ).asBool();
    short tetMode = data.inputValue(aTetMode).asShort();
    short constraintMode = data.inputValue(aConstraintMode).asShort();
    short numIter = data.inputValue( aIteration ).asShort();
    double constraintWeight = data.inputValue( aConstraintWeight ).asDouble();
    double constraintRadius = data.inputValue( aConstraintRadius ).asDouble();
    mesh.transWeight = data.inputValue( aTransWeight ).asDouble();
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
    // load mesh points
    MPointArray Mpts;
    itGeo.allPositions(Mpts);
    int numPts=Mpts.length();

    // when cage mode is changed
    if(!data.isClean(aARAP) || arap_flag){
        // prepare mesh tetrahedra
        pts.resize(numPts);
        for(int j=0; j<numPts; j++ ){
            Mpts[j] *= localToWorldMatrix;
            pts[j] << Mpts[j].x, Mpts[j].y, Mpts[j].z;
        }
        // make tetrahedral structure
        std::vector<Vector3d> tetCenter;
        getMeshData(data, input, inputGeom, mIndex, tetMode, pts, mesh.tetList, faceList, edgeList, vertexList, mesh.tetMatrix, mesh.tetWeight);
        mesh.dim = removeDegenerate(tetMode, numPts, mesh.tetList, faceList, edgeList, vertexList, mesh.tetMatrix);
        makeTetMatrix(tetMode, pts, mesh.tetList, faceList, edgeList, vertexList, mesh.tetMatrix, mesh.tetWeight);
        makeTetCenterList(tetMode, pts, mesh.tetList, tetCenter);
        mesh.numTet = (int)mesh.tetList.size()/4;
        mesh.computeTetMatrixInverse();
        
        // prepare cage tetrahedra
        std::vector<int> cageFaceCount(0);
        if (cageMode == TM_FACE || cageMode == TM_EDGE || cageMode == TM_VERTEX || cageMode == TM_VFACE){
            // face list
            makeFaceList( initCageMesh, cageFaceList, cageFaceCount, symmetricFace);
            std::vector<Matrix4d> initCageMatrix;
            makeFaceList(initCageMesh, cageFaceList, cageFaceCount);
            makeVertexList(initCageMesh, cageVertexList);
            makeEdgeList(cageFaceList, cageEdgeList);
            makeTetList(cageMode, numCagePts, cageFaceList, cageEdgeList, cageVertexList, cageTetList);
            makeTetMatrix(cageMode, initCagePts, cageTetList, cageFaceList, cageEdgeList,
                              cageVertexList, initCageMatrix, cageTetWeight, normaliseTet);
            numCageTet = (int) cageTetList.size()/4;
            cageMatrixI.resize(numCageTet);
            for(int i=0;i<numCageTet;i++){
                cageMatrixI[i] = initCageMatrix[i].inverse();
            }
        }else{
            numCageTet = numCagePts;
        }
        
        if(!areaWeighted){
            mesh.tetWeight.assign(mesh.numTet,1.0);
        }
        
        
        D.setNum(numCageTet, numPts, mesh.numTet);
        D.computeCageDistPts(cageMode, pts, initCagePts, cageTetList);
        D.computeCageDistTet(cageMode, tetCenter, initCagePts, cageTetList);
        D.findClosestPts();
        D.findClosestTet();
        
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
        }else{
            cageTetWeight.assign(numCageTet,1.0);
        }
    

        // find constraint points
        constraint.resize(3*numCageTet);
        std::vector<double> strength(numCageTet);
        for(int i=0;i<numCageTet;i++){
            strength[i] = constraintWeight*cageTetWeight[i];
            constraint[3*i] = T(i,mesh.tetList[4*D.closestTet[i]],strength[i]);
            constraint[3*i+1] = T(i,mesh.tetList[4*D.closestTet[i]+1],strength[i]);
            constraint[3*i+2] = T(i,mesh.tetList[4*D.closestTet[i]+2],strength[i]);
        }
        if( constraintMode == CONSTRAINT_NEIGHBOUR ){
            for(int i=0;i<numCageTet;i++){
                for(int j=0;j<numPts;j++){
                    if(D.distPts[i][j]<constraintRadius){
                        constraint.push_back(T(i,j,strength[i] * pow((constraintRadius-D.distPts[i][j])/constraintRadius,normExponent)));
                    }
                }
            }
        }
        int numConstraint=constraint.size();
        mesh.constraintWeight.resize(numConstraint);
        mesh.constraintVal.resize(numConstraint,numCageTet);
        for(int cur=0;cur<numConstraint;cur++){
            mesh.constraintWeight[cur] = std::make_pair(constraint[cur].col(), constraint[cur].value());
        }
        
        // weight computation
        w.resize(mesh.numTet);
        for(int j=0;j<mesh.numTet;j++){
            w[j].resize(numCageTet);
        }
        if(weightMode == WM_INV_DISTANCE){
            for(int j=0;j<mesh.numTet;j++){
                for( int i=0; i<numCageTet; i++){
                    w[j][i] = cageTetWeight[i]/pow(D.distTet[i][j],normExponent);
                }
            }
        }else if(weightMode == WM_CUTOFF_DISTANCE){
            double delta; // avoid under determined system for MLS
            delta = (cageMode & CM_MLS) ? EPSILON : 0;
            for(int j=0;j<mesh.numTet;j++){
                for( int i=0; i<numCageTet; i++){
                    w[j][i] = (D.distTet[i][j] > effectRadius)
                    ? delta : cageTetWeight[i]*pow((effectRadius-D.distTet[i][j])/effectRadius,normExponent);
                }
            }
        }else if(weightMode & WM_HARMONIC){
            Laplacian harmonicWeighting;
            makeFaceTet(data, input, inputGeom, mIndex, pts, harmonicWeighting.tetList, harmonicWeighting.tetMatrix, harmonicWeighting.tetWeight);
            harmonicWeighting.numTet = (int)harmonicWeighting.tetList.size()/4;
            std::vector<T> weightConstraint(numCageTet);
            // the vertex closest to the cage tet
            for(int i=0;i<numCageTet;i++){
                weightConstraint[i]=T(i,D.closestPts[i],cageTetWeight[i]);
            }
            // vertices within effectRadius
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
            if( weightMode == WM_HARMONIC_ARAP ){
                harmonicWeighting.computeTetMatrixInverse();
                harmonicWeighting.dim = numPts + harmonicWeighting.numTet;
                isError = harmonicWeighting.ARAPprecompute();
            }else if(weightMode == WM_HARMONIC_COTAN){
                harmonicWeighting.dim = numPts;
                isError = harmonicWeighting.cotanPrecompute();
            }
            if(isError>0) return MS::kFailure;
            std::vector< std::vector<double> > w_tet(numCageTet);
            harmonicWeighting.harmonicSolve();
            for(int i=0;i<numCageTet;i++){
                makeTetWeightList(tetMode, mesh.tetList, faceList, edgeList, vertexList, harmonicWeighting.Sol.col(i), w_tet[i]);
                for(int j=0;j<mesh.numTet;j++){
                    w[j][i] = w_tet[i][j];
                }
            }
        }
        // normalise weights
        if (data.inputValue( aPositiveWeight ).asBool()){
            for(int j=0;j<mesh.numTet;j++){
                for (int i = 0; i < numCageTet; i++){
                    w[j][i] = max(w[j][i], 0.0);
                }
            }
        }
        bool normaliseWeight = data.inputValue( aNormaliseWeight ).asBool();
        for(int j=0;j<mesh.numTet;j++){
            double sum = std::accumulate(w[j].begin(), w[j].end(), 0.0);
            if ((sum > 1 || normaliseWeight || cageMode & CM_MLS ) && sum>0){
                for (int i = 0; i < numCageTet; i++){
                    w[j][i] /= sum;
                }
            }
        }
        
        // prepare ARAP solver
        mesh.ARAPprecompute();
        status = data.setClean(aARAP);
    }
//////  END of precomputation ///////
    if(isError>0){
        return MS::kFailure;
    }
    

    B.setNum(numCageTet);
    //  find affine transformations for tetrahedra
    cageMatrix.resize(numCageTet); A.resize(mesh.numTet); blendedSE.resize(mesh.numTet);
    blendedR.resize(mesh.numTet); blendedS.resize(mesh.numTet); blendedL.resize(mesh.numTet);
    
    if(cageMode & CM_MLS){
        for(int j = 0; j < mesh.numTet; j++){
            // barycentre of the original and the current cage points
            Vector3d icenter = Vector3d::Zero();
            Vector3d center = Vector3d::Zero();
            for(int i=0;i<numCagePts;i++){
                icenter += w[j][i] * initCagePts[i];
                center += w[j][i] * cagePts[i];
            }
            // determine the closest transformation
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
            A[j] = pad(M, center-icenter);
            
//            if(cageMode == CM_MLS_AFF){
//                M = (q*p.transpose()) * (p*p.transpose()).inverse();
//            }else if(cageMode == CM_MLS_SIM){
//                M = umeyama(p,q).block(0,0,3,3);
//            }else if(cageMode == CM_MLS_RIGID){
//                M = umeyama(p,q,false).block(0,0,3,3);
//            }
//            A[j] = pad(M.transpose(), center-icenter);

        }
        for(int i=0; i<numCageTet; i++){
            B.Aff[i]=pad(Matrix3d::Identity().eval(),cagePts[i]-initCagePts[i]);
        }
    }else{
        makeTetMatrix(cageMode, cagePts, cageTetList, cageFaceList, cageEdgeList,
                      cageVertexList, cageMatrix, dummy_weight, normaliseTet);
        for(int i=0; i<numCageTet; i++){
            B.Aff[i]=cageMatrixI[i]*cageMatrix[i];
        }
        B.parametrise(blendMode);
        
        // compute blended matrices
#pragma omp parallel for
        for (int j = 0; j < mesh.numTet; j++){
            // blend matrix
            if (blendMode == BM_SRL){
                blendedS[j] = expSym(blendMat(B.logS, w[j]));
                Vector3d l = blendMat(B.L, w[j]);
                blendedR[j] = frechetSum ? frechetSO(B.R, w[j]) : expSO(blendMat(B.logR, w[j]));
                A[j] = pad(blendedS[j]*blendedR[j], l);
            }
            else if (blendMode == BM_SSE){
                blendedS[j] = expSym(blendMat(B.logS, w[j]));
                blendedSE[j] = expSE(blendMat(B.logSE, w[j]));
                A[j] = pad(blendedS[j], Vector3d::Zero()) * blendedSE[j];
            }
            else if (blendMode == BM_LOG3){
                blendedR[j] = blendMat(B.logGL, w[j]).exp();
                Vector3d l = blendMat(B.L, w[j]);
                A[j] = pad(blendedR[j], l);
            }
            else if (blendMode == BM_LOG4){
                A[j] = blendMat(B.logAff, w[j]).exp();
            }
            else if (blendMode == BM_SQL){
                Vector4d q = blendQuat(B.quat, w[j]);
                Vector3d l = blendMat(B.L, w[j]);
                blendedS[j] = blendMatLin(B.S, w[j]);
                Quaternion<double> Q(q);
                blendedR[j] = Q.matrix().transpose();
                A[j] = pad(blendedS[j]*blendedR[j], l);
            }
            else if (blendMode == BM_AFF){
                A[j] = blendMatLin(B.Aff, w[j]);
            }
        }
    }
    
    // set constraint
    int numConstraints = constraint.size();
    mesh.constraintVal.resize(numConstraints,3);
    RowVector4d cv;
    for(int cur=0;cur<numConstraints;cur++){
        cv = pad(pts[constraint[cur].col()]) * B.Aff[constraint[cur].row()];
        mesh.constraintVal(cur,0) = cv[0];
        mesh.constraintVal(cur,1) = cv[1];
        mesh.constraintVal(cur,2) = cv[2];
    }
    
    // iterate to determine vertices position
    for(int k=0;k<numIter;k++){
        // solve ARAP
        mesh.ARAPSolve(A);
        // set new vertices position
        new_pts.resize(numPts);
        for(int i=0;i<numPts;i++){
            new_pts[i][0]=mesh.Sol(i,0);
            new_pts[i][1]=mesh.Sol(i,1);
            new_pts[i][2]=mesh.Sol(i,2);
        }
        // if iteration continues
        if(k+1<numIter){
            makeTetMatrix(tetMode, new_pts, mesh.tetList, faceList, edgeList, vertexList, Q, dummy_weight);
            Matrix3d S,R,newS,newR;
            if(blendMode == BM_AFF || blendMode == BM_LOG4 || blendMode == BM_LOG3 || cageMode & CM_MLS){
                for(int i=0;i<mesh.numTet;i++){
                    polarHigham(A[i].block(0,0,3,3), blendedS[i], blendedR[i]);
                }
            }
#pragma omp parallel for
            for(int i=0;i<mesh.numTet;i++){
                polarHigham((mesh.tetMatrixInverse[i]*Q[i]).block(0,0,3,3), newS, newR);
//                tetEnergy[i] = (newS-blendedS[i]).squaredNorm();
                A[i].block(0,0,3,3) = blendedS[i]*newR;
            }
        }
    }
    for(int i=0;i<numPts;i++){
        Mpts[i].x=mesh.Sol(i,0);
        Mpts[i].y=mesh.Sol(i,1);
        Mpts[i].z=mesh.Sol(i,2);
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
    
    aCageMode = eAttr.create( "cageMode", "cgm", TM_VERTEX );
    eAttr.addField( "face", TM_FACE );
    eAttr.addField( "edge", TM_EDGE );
    eAttr.addField( "vertex", TM_VERTEX );
    eAttr.addField( "vface", TM_VFACE );
//    eAttr.addField( "MVC", CM_MVC );
    eAttr.addField( "MLS-AFF", CM_MLS_AFF );
    eAttr.addField( "MLS-SIM", CM_MLS_SIM );
    eAttr.addField( "MLS-RIGID", CM_MLS_RIGID );
    addAttribute( aCageMode );
    attributeAffects( aCageMode, outputGeom );
    attributeAffects( aCageMode, aARAP);
    
    
    aConstraintMode = eAttr.create( "constraintMode", "ctm", CONSTRAINT_CLOSEST );
    eAttr.addField( "neighbour",  CONSTRAINT_NEIGHBOUR);
    eAttr.addField( "closestPt", CONSTRAINT_CLOSEST );
    eAttr.setStorable(true);
    addAttribute( aConstraintMode );
    attributeAffects( aConstraintMode, outputGeom );
    attributeAffects( aConstraintMode, aARAP);
    
    aConstraintWeight = nAttr.create("constraintWeight", "cw", MFnNumericData::kDouble, 1e-10);
    nAttr.setStorable(true);
	addAttribute( aConstraintWeight );
	attributeAffects( aConstraintWeight, outputGeom );
    attributeAffects( aConstraintWeight, aARAP);
    
    aConstraintRadius = nAttr.create("constraintRadius", "cr", MFnNumericData::kDouble, 2.0);
    nAttr.setStorable(true);
	addAttribute( aConstraintRadius );
	attributeAffects( aConstraintRadius, outputGeom );
    attributeAffects( aConstraintRadius, aARAP );
    
    aTransWeight = nAttr.create("translationWeight", "tw", MFnNumericData::kDouble, 0);
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

    aPositiveWeight = nAttr.create( "positiveWeight", "posw", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aPositiveWeight );
    attributeAffects( aPositiveWeight, outputGeom );
    attributeAffects( aPositiveWeight, aARAP );

    aWeightMode = eAttr.create( "weightMode", "wtm", WM_INV_DISTANCE );
    eAttr.addField( "inverse", WM_INV_DISTANCE );
    eAttr.addField( "cut-off", WM_CUTOFF_DISTANCE );
//    eAttr.addField( "harmonic-arap", WM_HARMONIC);
    eAttr.addField( "harmonic-cotan", WM_HARMONIC_COTAN);
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

    aNeighbourWeighting = nAttr.create( "neighbourWeighting", "nghbrw", MFnNumericData::kBoolean, false );
    nAttr.setStorable(true);
    addAttribute( aNeighbourWeighting );
    attributeAffects( aNeighbourWeighting, outputGeom );
    attributeAffects( aNeighbourWeighting, aARAP );
    
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
