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

// (small) weight in ARAP for the translation part
#define TRANSWEIGHT 0.0001f

using namespace Eigen;
using namespace AffineLib;

MTypeId CageDeformerNode::id( 0x00000202 );
MString CageDeformerNode::nodeName("cageARAP");
MObject CageDeformerNode::aCageMesh;
MObject CageDeformerNode::aCageMode;
MObject CageDeformerNode::aBlendMode;
MObject CageDeformerNode::aRotationConsistency;
MObject CageDeformerNode::aFrechetSum;
MObject CageDeformerNode::aConstraintMode;
MObject CageDeformerNode::aConstraintWeight;

double distPtLin(Vector3d p,Vector3d a,Vector3d b){
    /// compute distance between a line segment and a point
    double t= (a-b).dot(p-b)/(a-b).squaredNorm();
    if(t>1){
        return (a-p).squaredNorm();
    }else if(t<0){
        return (b-p).squaredNorm();
    }else{
        return (t*(a-b)-(p-b)).squaredNorm();
    }
}

double distPtTri(Vector3d p, Matrix4d m){
    /// compute distance between a triangle and a point
    double s[4];
    Vector3d a,b,c,n;
    a << m(0,0), m(0,1), m(0,2);
    b << m(1,0), m(1,1), m(1,2);
    c << m(2,0), m(2,1), m(2,2);
    n << m(3,0), m(3,1), m(3,2);
    double k=(n-a).dot(a-p);
    if(k<0) return HUGE_VALF;
    s[0]=distPtLin(p,a,b);
    s[1]=distPtLin(p,b,c);
    s[2]=distPtLin(p,c,a);
    Matrix3d A;
    A << b(0)-a(0), c(0)-a(0), n(0)-a(0),
    b(1)-a(1), c(1)-a(1), n(1)-a(1),
    b(2)-a(2), c(2)-a(2), n(2)-a(2);
    Vector3d v = A.inverse()*(p-a);
    if(v(0)>0 && v(1)>0 && v(0)+v(1)<1){
        s[3]=k*k;
    }else{
        s[3] = HUGE_VALF;
    }
    return min(min(min(s[0],s[1]),s[2]),s[3]);
}

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
    // load cage mesh and other attributes
    MObject oCageMesh = data.inputValue( aCageMesh ).asMesh();
    short blendMode = data.inputValue(aBlendMode).asShort();
	bool rotationCosistency = data.inputValue( aRotationConsistency ).asBool();
	bool frechetSum = data.inputValue( aFrechetSum ).asBool();
    short newConstraintMode = data.inputValue(aConstraintMode).asShort();
    double newConstraintWeight = data.inputValue( aConstraintWeight ).asDouble();
    if ( oCageMesh.isNull() || blendMode == 99)
        return MS::kSuccess;
    short newCageMode = data.inputValue(aCageMode).asShort();
    MFnMesh fnCageMesh( oCageMesh, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MPointArray cagePoints;
    fnCageMesh.getPoints( cagePoints,  MSpace::kWorld );
    // save initial cage state
    if (initCagePoints.length() != cagePoints.length()){
        initCageMesh = oCageMesh;
        initCagePoints=cagePoints;
    }
    // when cage mode is changed
    if(newCageMode != cageMode || newConstraintMode != constraintMode || newConstraintWeight != constraintWeight)
    {
        cageMode = newCageMode;
        constraintMode = newConstraintMode;
        constraintWeight = newConstraintWeight;
	    std::vector<double> tetWeight;
        // read target mesh data
        MArrayDataHandle hInput = data.outputArrayValue( input, &status );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        status = hInput.jumpToElement( mIndex );
        CHECK_MSTATUS_AND_RETURN_IT( status );
        MObject oInputGeom = hInput.outputValue().child( inputGeom ).asMesh();
        MFnMesh inputMesh(oInputGeom);
        inputMesh.getPoints( pts );
		numPts=pts.length();
        for(int j=0; j<numPts; j++ )
            pts[j] *= localToWorldMatrix;
        MIntArray count;
        inputMesh.getTriangles( count, meshTriangles );
		numTet=meshTriangles.length()/3;
		std::vector<Matrix4d> P(numTet);
        tetCenter.resize(numTet);
        tetMatrixC(pts, meshTriangles, P, tetCenter);
        PI.resize(numTet);
		for(int i=0;i<numTet;i++)
			PI[i] = P[i].inverse();
        // prepare cage tetrahedra
        MFnMesh fnInitCageMesh( initCageMesh, &status );
        if(cageMode == 10 || cageMode == 11)  // face mode
        {
			if(cageMode == 10){       // triangulate faces by MAYA standard
                MIntArray count;
                fnInitCageMesh.getTriangles( count, triangles );
                tetWeight.resize(triangles.length()/3, 1.0f);
			}else if(cageMode ==11){  // trianglate faces with more than 3 edges in a symmetric way
				triangles.clear();
				MItMeshPolygon iter(initCageMesh);
				MIntArray tmp;
                MVector normal;
				tetWeight.reserve(4*iter.count());
                unsigned int l;
				for(unsigned int i=0; ! iter.isDone(); i++){
					iter.getVertices(tmp);
					l=tmp.length();
					if(l==3){
						tetWeight.push_back(1.0);
						triangles.append(tmp[0]);
						triangles.append(tmp[1]);
						triangles.append(tmp[2]);
					}else{
						for(unsigned int j=0;j<l;j++){
                            tetWeight.push_back((l-2.0)/l);
							triangles.append(tmp[j]);
							triangles.append(tmp[(j+1) % l]);
							triangles.append(tmp[(j+2) % l]);
						}
					}
					iter.next();
				}
            }
            // face mode compute init matrix
            numPrb=triangles.length()/3;
            initMatrix.resize(numPrb);
            tetMatrix(initCagePoints, triangles, cageMode, initMatrix);
            // compute weight
            w.resize(numTet);
            std::vector< std::vector<double> > idist(numTet);
            for(int j=0;j<numTet;j++){
                idist[j].resize(numPrb);
                w[j].resize(numPrb);
                double sidist = 0.0;
                for(int i=0;i<numPrb;i++){
                    idist[j][i] = tetWeight[i]/distPtTri(tetCenter[j],initMatrix[i]);
                    sidist += idist[j][i];
                }
                assert(sidist>0.0f);
                for(int i=0;i<numPrb;i++)
                    w[j][i] = idist[j][i] /sidist;
            }// face mode end
        }else if(cageMode == 0 || cageMode == 1){   // vertex mode
            triangles.clear();
            std::vector<int> tetCount(initCagePoints.length());
            MItMeshVertex iter(initCageMesh);
            for(int j=0; ! iter.isDone(); j++){
                MIntArray v;
                iter.getConnectedVertices(v);     // at each vertex, construct tetrahedra from connected edges
                int l=v.length();
                if(l==3){
                    if(isDegenerate(initCagePoints[j],initCagePoints[v[0]],initCagePoints[v[1]],initCagePoints[v[2]]) != 0){
                        tetCount[j]++;
                        triangles.append(j);
                        triangles.append(v[0]);
                        triangles.append(v[1]);
                        triangles.append(v[2]);
                    }
                }else{
                    for(int k=0;k<l;k++){
                        if(isDegenerate(initCagePoints[j],initCagePoints[v[k]],initCagePoints[v[(k+1) % l]],initCagePoints[v[(k+2) % l]]) != 0){
                            tetCount[j]++;
                            triangles.append(j);
                            triangles.append(v[k]);
                            triangles.append(v[(k+1) % l]);
                            triangles.append(v[(k+2) % l]);
                        }
                    }
                }
                iter.next();
            }
            numPrb=triangles.length()/4;
            initMatrix.resize(numPrb);
            tetMatrix(initCagePoints, triangles, cageMode, initMatrix);
            // vertex mode compute weight
            w.resize(numTet);
            std::vector< std::vector<double> > idist(numTet);
            tetWeight.resize(numPrb);
            for(int i=0;i<numPrb;i++)
                tetWeight[i]=1.0/(double)tetCount[triangles[4*i]];
            for(int j=0;j<numTet;j++){
                idist[j].resize(numPrb);
                w[j].resize(numPrb);
                double sidist = 0.0;
                for(int i=0;i<numPrb;i++){
                    Vector3d c(initCagePoints[triangles[4*i]].x,initCagePoints[triangles[4*i]].y,initCagePoints[triangles[4*i]].z);
                    idist[j][i] = tetWeight[i] / ((tetCenter[j]-c).squaredNorm());
                    sidist += idist[j][i];
                }
                assert(sidist>0.0f);
                for(int i=0;i<numPrb;i++)
                    w[j][i] = idist[j][i] /sidist;
            }
        }else if(cageMode == 5 || cageMode == 6 ){ // vertex averaged normal mode
            triangles.clear();
            std::vector<int> tetCount(initCagePoints.length());
            MItMeshVertex iter(initCageMesh);
            for(int j=0; ! iter.isDone(); j++){
                MIntArray v;
                iter.getConnectedVertices(v);
                int l=v.length();
                for(int k=0;k<l;k++){
                    tetCount[j]++;
                    triangles.append(j);
                    triangles.append(v[k]);
                    triangles.append(v[(k+1) % l]);
                }
                iter.next();
            }
            numPrb=triangles.length()/3;
            initMatrix.resize(numPrb);
            tetMatrix(initCagePoints, triangles, cageMode, initMatrix);
            // vertex mode compute weight
            w.resize(numTet);
            std::vector< std::vector<double> > idist(numTet);
            tetWeight.resize(numPrb);
            for(int i=0;i<numPrb;i++)
                tetWeight[i]=1.0/(double)tetCount[triangles[3*i]];
            for(int j=0;j<numTet;j++){
                idist[j].resize(numPrb);
                w[j].resize(numPrb);
                double sidist = 0.0;
                for(int i=0;i<numPrb;i++){
                    Vector3d c(initCagePoints[triangles[3*i]].x,initCagePoints[triangles[3*i]].y,initCagePoints[triangles[3*i]].z);
                    idist[j][i] = tetWeight[i] / ((tetCenter[j]-c).squaredNorm());
                    sidist += idist[j][i];
                }
                assert(sidist>0.0f);
                for(int i=0;i<numPrb;i++)
                    w[j][i] = idist[j][i] /sidist;
            }
        }// end of cage setup
        
        // find constraint points
        if(constraintMode == 1){
            numConstraint = numPrb;
        }else{
            numConstraint = 1;    // at least one constraint is necessary to determine global translation
        }
        constraintTet.resize(numConstraint);
        constraintVector.resize(numConstraint);
        // for each cage tetrahedra, constraint the point on the mesh with largest weight
        for(int i=0;i<numConstraint;i++){
            constraintTet[i] = 0;
            for(int j=1;j<numTet;j++){
                if(w[j][i] > w[constraintTet[i]][i]){
                    constraintTet[i] = j;
                }
            }
            constraintVector[i] << tetCenter[constraintTet[i]](0), tetCenter[constraintTet[i]](1), tetCenter[constraintTet[i]](2), 1.0;
        }
        // precompute arap solver
        arapHI(PI, meshTriangles);
    }
    // compute deformation
    if( ! rotationCosistency || numPrb != prevNs.size()){        // clear previous rotation
        prevThetas.clear();
        prevThetas.resize(numPrb, 0.0);
        prevNs.clear();
        prevNs.resize(numPrb, Vector3d::Zero());
    }
    //  find affine transformations for tetrahedra
    std::vector<Matrix4d> cageMatrix(numPrb), SE(numPrb), logSE(numPrb),logAff(numPrb),aff(numPrb);
    std::vector<Matrix3d> logR(numPrb),R(numPrb),logS(numPrb),logGL(numPrb);
    std::vector<Vector3d> L(numPrb);
    std::vector<Vector4d> quat(numPrb);
    tetMatrix(cagePoints, triangles, cageMode, cageMatrix);
    for(int i=0; i<numPrb; i++)
        aff[i]=initMatrix[i].inverse()*cageMatrix[i];
    // compute parametrisation
    if(blendMode == 0 || blendMode == 1 || blendMode == 5)  // polarexp or quaternion
    {
        for(unsigned int i=0;i<numPrb;i++){
            parametriseGL(aff[i].block(0,0,3,3), logS[i] ,R[i]);
            L[i] = transPart(aff[i]);
            if(blendMode == 0){  // Rotational log
                logR[i]=logSOc(R[i], prevThetas[i], prevNs[i]);
            }else if(blendMode == 1){ // Eucledian log
                SE[i]=affine(R[i], L[i]);
                logSE[i]=logSEc(SE[i], prevThetas[i], prevNs[i]);
            }else if(blendMode == 5){ // quaternion
                Quaternion<double> Q(R[i].transpose());
                quat[i] << Q.x(), Q.y(), Q.z(), Q.w();
            }
        }
    }else if(blendMode == 2){    //logmatrix3
        for(unsigned int i=0;i<numPrb;i++){
            logGL[i] = aff[i].block(0,0,3,3).log();
            L[i] = transPart(aff[i]);
        }
    }else if(blendMode == 3){   // logmatrix4
        for(unsigned int i=0;i<numPrb;i++){
            logAff[i] = aff[i].log();
        }
    }
    // compute blended matrices
#pragma omp parallel for
    std::vector<Matrix4d> At(numTet);
    for(int j=0; j<numTet; j++ ){
        if(blendMode==0){
            Matrix3d RR=Matrix3d::Zero();
            Matrix3d SS=Matrix3d::Zero();
            Vector3d l=Vector3d::Zero();
            for(unsigned int i=0; i<numPrb; i++){
                RR += w[j][i] * logR[i];
                SS += w[j][i] * logS[i];
                l += w[j][i] * L[i];
            }
            SS = expSym(SS);
            if(frechetSum){
                RR = frechetSO(R, w[j]);
            }else{
                RR = expSO(RR);
            }
            At[j] = affine(SS*RR, l);
        }else if(blendMode==1){    // rigid transformation
            Matrix4d EE=Matrix4d::Zero();
            Matrix3d SS=Matrix3d::Zero();
            for(unsigned int i=0; i<numPrb; i++){
                EE +=  w[j][i] * logSE[i];
                SS +=  w[j][i] * logS[i];
            }
            if(frechetSum){
                EE = frechetSE(SE, w[j]);
            }else{
                EE = expSE(EE);
            }
            At[j] = affine(expSym(SS),Vector3d::Zero())*EE;
        }else if(blendMode == 2){    //logmatrix3
            Matrix3d G=Matrix3d::Zero();
            Vector3d l=Vector3d::Zero();
            for(unsigned int i=0; i<numPrb; i++){
                G +=  w[j][i] * logGL[i];
                l += w[j][i] * L[i];
            }
            At[j] = affine(G.exp(), l);
        }else if(blendMode == 3){   // logmatrix4
            Matrix4d A=Matrix4d::Zero();
            for(unsigned int i=0; i<numPrb; i++)
                A +=  w[j][i] * logAff[i];
            At[j] = A.exp();
        }else if(blendMode == 5){ // quaternion
            Vector4d q=Vector4d::Zero();
            Matrix3d SS=Matrix3d::Zero();
            Vector3d l=Vector3d::Zero();
            for(unsigned int i=0; i<numPrb; i++){
                q += w[j][i] * quat[i];
                SS += w[j][i] * logS[i];
                l += w[j][i] * L[i];
            }
            SS = expSym(SS);
            Quaternion<double> Q(q);
            Matrix3d RR = Q.matrix().transpose();
            At[j] = affine(SS*RR, l);
        }else if(blendMode==10){
            At[j] = Matrix4d::Zero();
            for(unsigned int i=0; i<numPrb; i++){
                At[j] += w[j][i] * aff[i];
            }
        }
    }
    
    // compute target vertices position
    MatrixXd G=MatrixXd::Zero(numTet+numPts,3);
    arapG(At, PI, meshTriangles, aff, G);
    MatrixXd Sol = solver.solve(G);
    for(unsigned int i=0;i<numPts;i++){
        pts[i].x=Sol(i,0);
        pts[i].y=Sol(i,1);
        pts[i].z=Sol(i,2);
        pts[i] *= localToWorldMatrix.inverse();
    }
    itGeo.setAllPositions(pts);
    return MS::kSuccess;
}

void CageDeformerNode::tetMatrix(const MPointArray& p, const MIntArray& triangles, short cageMode, std::vector<Matrix4d>& m)
/// prepare cage tetrahedra from points and triangulation
{
    if(cageMode == 10 || cageMode == 11){
        MVector u, v, q ;
        for(unsigned int i=0;i<triangles.length()/3;i++){
            u=MVector(p[triangles[3*i+1]]-p[triangles[3*i]]);
            v=MVector(p[triangles[3*i+2]]-p[triangles[3*i]]);
			q=u^v;
			q.normalize();
            
            m[i] << p[triangles[3*i]].x, p[triangles[3*i]].y, p[triangles[3*i]].z, 1,
            p[triangles[3*i+1]].x, p[triangles[3*i+1]].y, p[triangles[3*i+1]].z, 1,
            p[triangles[3*i+2]].x, p[triangles[3*i+2]].y, p[triangles[3*i+2]].z, 1,
            (p[triangles[3*i]]+q).x, (p[triangles[3*i]]+q).y, (p[triangles[3*i]]+q).z,1;
        }
    }else if (cageMode == 0){
        for(unsigned int i=0;i<triangles.length()/4;i++){
            m[i] << p[triangles[4*i]].x, p[triangles[4*i]].y, p[triangles[4*i]].z, 1,
            p[triangles[4*i+1]].x, p[triangles[4*i+1]].y, p[triangles[4*i+1]].z, 1,
            p[triangles[4*i+2]].x, p[triangles[4*i+2]].y, p[triangles[4*i+2]].z, 1,
            p[triangles[4*i+3]].x, p[triangles[4*i+3]].y, p[triangles[4*i+3]].z, 1;
        }
    }else if (cageMode == 1){
        MVector u, v, q;
        for(unsigned int i=0;i<triangles.length()/4;i++){
            u=MVector(p[triangles[4*i+1]]-p[triangles[4*i]]);
            v=MVector(p[triangles[4*i+2]]-p[triangles[4*i]]);
            q=MVector(p[triangles[4*i+3]]-p[triangles[4*i]]);
            u.normalize();
            v.normalize();
            q.normalize();
            m[i] << p[triangles[4*i]].x, p[triangles[4*i]].y, p[triangles[4*i]].z, 1,
            u[0]+p[triangles[4*i]].x,u[1]+p[triangles[4*i]].y,u[2]+p[triangles[4*i]].z,1,
            v[0]+p[triangles[4*i]].x,v[1]+p[triangles[4*i]].y,v[2]+p[triangles[4*i]].z,1,
            q[0]+p[triangles[4*i]].x,q[1]+p[triangles[4*i]].y,q[2]+p[triangles[4*i]].z,1;
        }
    }else if (cageMode == 5 || cageMode == 6){
        std::vector<Vector3d> normal(p.length(),Vector3d::Zero());           //averaged normal vector for vertex mode
        for(unsigned int i=0;i<triangles.length()/3;i++){
            MVector q = MVector(p[triangles[3*i+1]]-p[triangles[3*i]]) ^ MVector(p[triangles[3*i+2]]-p[triangles[3*i]]);
            normal[triangles[3*i]] += Vector3d(q[0], q[1], q[2]);
        }
        if (cageMode == 5){
            for(unsigned int i=0;i<triangles.length()/3;i++){
                m[i] << p[triangles[3*i]].x, p[triangles[3*i]].y, p[triangles[3*i]].z, 1,
                p[triangles[3*i+1]].x, p[triangles[3*i+1]].y, p[triangles[3*i+1]].z, 1,
                p[triangles[3*i+2]].x, p[triangles[3*i+2]].y, p[triangles[3*i+2]].z, 1,
                p[triangles[3*i]].x + normal[triangles[3*i]][0], p[triangles[3*i]].y + normal[triangles[3*i]][1], p[triangles[3*i]].z + normal[triangles[3*i]][2],1;
            }
        }else if (cageMode == 6){
            MVector u,v,q;
            for(unsigned int i=0;i<triangles.length()/3;i++){
                u=MVector(p[triangles[3*i+1]]-p[triangles[3*i]]);
                v=MVector(p[triangles[3*i+2]]-p[triangles[3*i]]);
                u.normalize();
                v.normalize();
                normal[triangles[3*i]].normalize();
                m[i] << p[triangles[3*i]].x, p[triangles[3*i]].y, p[triangles[3*i]].z, 1,
                u[0]+p[triangles[3*i]].x,u[1]+p[triangles[3*i]].y,u[2]+p[triangles[3*i]].z,1,
                v[0]+p[triangles[3*i]].x,v[1]+p[triangles[3*i]].y,v[2]+p[triangles[3*i]].z,1,
                normal[triangles[3*i]][0]+p[triangles[3*i]].x,normal[triangles[3*i]][1]+p[triangles[3*i]].y,normal[triangles[3*i]][2]+p[triangles[3*i]].z,1;
            }
        }
    }
}


void CageDeformerNode::arapHI(const std::vector<Matrix4d>& PI, const MIntArray& tr)
/// pre-compute arap solver
{
    int dim = numTet + numPts;
    std::vector<T> tripletList;
    tripletList.reserve(numTet*16+numConstraint*6);
    Matrix4d Hlist;
	Matrix4d diag=Matrix4d::Identity();
	diag(3,3)=TRANSWEIGHT;
    int s,t;
	for(unsigned int i=0;i<numTet;i++){
		Hlist=PI[i].transpose()*diag*PI[i];
		for(unsigned int j=0;j<4;j++){
            if(j==3){
                s=numPts+i;
            }else{
                s=tr[3*i+j];
            }
			for(unsigned int k=0;k<4;k++){
                if(k==3){
                    t=numPts+i;
                }else{
                    t=tr[3*i+k];
                }
                tripletList.push_back(T(t,s,Hlist(j,k)));
			}
		}
	}
    //    // set hard constraint
    //    for(int i=0;i<numConstraint;i++){
    //        tripletList.push_back(T(dim+i,tr[3*constraintTet[i]],1.0/3.0));
    //        tripletList.push_back(T(dim+i,tr[3*constraintTet[i]+1],1.0/3.0));
    //        tripletList.push_back(T(dim+i,tr[3*constraintTet[i]+2],1.0/3.0));
    //        tripletList.push_back(T(tr[3*constraintTet[i]],dim+i,1.0/3.0));
    //        tripletList.push_back(T(tr[3*constraintTet[i]+1],dim+i,1.0/3.0));
    //        tripletList.push_back(T(tr[3*constraintTet[i]+2],dim+i,1.0/3.0));
    //    }
    //    SpMat mat(dim+numConstraint, dim+numConstraint);
    SpMat mat(dim, dim);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    // set soft constraint
    std::vector<T> constraintList;
    constraintList.reserve(numConstraint*3);
    F.resize(dim,numConstraint);
    F.setZero();
    for(int i=0;i<numConstraint;i++){
        constraintList.push_back(T(tr[3*constraintTet[i]],i,1.0/3.0));
        constraintList.push_back(T(tr[3*constraintTet[i]+1],i,1.0/3.0));
        constraintList.push_back(T(tr[3*constraintTet[i]+2],i,1.0/3.0));
    }
    F.setFromTriplets(constraintList.begin(), constraintList.end());
    mat += constraintWeight * F * F.transpose();
    solver.compute(mat);
}

void CageDeformerNode::arapG(const std::vector<Matrix4d>& At, const std::vector<Matrix4d>& PI,
                             const MIntArray& tr, const std::vector<Matrix4d>& aff, MatrixXd& G)
/// compute run-time arap vector
{
    //    int dim = numPts+numTet;
    Matrix4d Glist;
    Matrix4d diag=Matrix4d::Identity();
    diag(3,3)=TRANSWEIGHT;
    for(unsigned int i=0;i<numTet;i++){
        Glist=At[i].transpose()*diag*PI[i];
        for(unsigned int k=0;k<3;k++){
            for(unsigned int j=0;j<3;j++){
                G(tr[3*i+j],k) += Glist(k,j);
            }
            G(numPts+i,k) += Glist(k,3);
        }
    }
    // set hard constraint
    //	MatrixXf G=MatrixXf::Zero(dim+numConstraint,3);
    //    for(int i=0;i<numConstraint;i++){
    //        RowVector4d cv = constraintVector[i]*aff[i];
    //        G.block(dim+i,0,1,3) << cv(0), cv(1), cv(2);
    //    }
    // set soft constraint
    std::vector<T> constraintList;
    constraintList.reserve(numConstraint*3);
    SpMat S(numConstraint,3);
    for(int i=0;i<numConstraint;i++){
        RowVector4d cv = constraintVector[i]*aff[i];
        constraintList.push_back(T(i,0,cv(0)));
        constraintList.push_back(T(i,1,cv(1)));
        constraintList.push_back(T(i,2,cv(2)));
    }
    S.setFromTriplets(constraintList.begin(), constraintList.end());
    SpMat FS = constraintWeight * F * S;
    G += MatrixXd(FS);
}

void CageDeformerNode::tetMatrixC(const MPointArray& p, const MIntArray& tr, std::vector<Matrix4d>& m, std::vector<Vector3d>& tetCenter)
/// prepare facial tetrahedra for target mesh
{
    MVector u, v, q;
    for(unsigned int i=0;i<numTet;i++)
    {
        tetCenter[i] << (p[tr[3*i]].x+p[tr[3*i+1]].x+p[tr[3*i+2]].x)/3.0,
        (p[tr[3*i]].y+p[tr[3*i+1]].y+p[tr[3*i+2]].y)/3.0,
        (p[tr[3*i]].z+p[tr[3*i+1]].z+p[tr[3*i+2]].z)/3.0;
        u=p[tr[3*i+1]]-p[tr[3*i]];
        v=p[tr[3*i+2]]-p[tr[3*i]];
        q=u^v;
        q.normalize();
        
        m[i] << p[tr[3*i]].x, p[tr[3*i]].y, p[tr[3*i]].z, 1,
        p[tr[3*i+1]].x, p[tr[3*i+1]].y, p[tr[3*i+1]].z, 1,
        p[tr[3*i+2]].x, p[tr[3*i+2]].y, p[tr[3*i+2]].z, 1,
        q[0]+p[tr[3*i]].x,q[1]+p[tr[3*i]].y, q[2]+p[tr[3*i]].z,1;
    }
}


// maya plugin initialization
MStatus CageDeformerNode::initialize()
{
    MFnTypedAttribute tAttr;
    MFnNumericAttribute nAttr;
    MFnEnumAttribute eAttr;
    
    aCageMesh = tAttr.create( "cageMesh", "cm", MFnData::kMesh );
    addAttribute( aCageMesh );
    attributeAffects( aCageMesh, outputGeom );
    
    aCageMode = eAttr.create( "cageMode", "cgm", 10 );
    eAttr.addField( "vertex", 0 );
    eAttr.addField( "vertexNormalized", 1 );
    eAttr.addField( "vertex avg. normal", 5 );
    eAttr.addField( "vertex avg. normal normalized", 6 );
    eAttr.addField( "face", 10 );
    eAttr.addField( "faceSymmetrized", 11 );
    //    eAttr.addField( "init", 99 );
    addAttribute( aCageMode );
    attributeAffects( aCageMode, outputGeom );
    
    aConstraintMode = eAttr.create( "constraintMode", "constraint", 0 );
    eAttr.addField( "none", 0 );
    eAttr.addField( "allFaces", 1 );
    addAttribute( aConstraintMode );
    attributeAffects( aConstraintMode, outputGeom );
    
    aConstraintWeight = nAttr.create("constraintWeight", "cw", MFnNumericData::kDouble, 1.0);
    nAttr.setStorable(true);
	addAttribute( aConstraintWeight );
	attributeAffects( aConstraintWeight, outputGeom );
    
    aBlendMode = eAttr.create( "blendMode", "bm", 0 );
    eAttr.addField( "polarexp", 0 );
    eAttr.addField( "polarexpSE", 1 );
    eAttr.addField( "logmatrix3", 2 );
    eAttr.addField( "logmatrix4", 3 );
    eAttr.addField( "quaternion", 5 );
    eAttr.addField( "linear", 10 );
    eAttr.addField( "off", 99 );
    addAttribute( aBlendMode );
    attributeAffects( aBlendMode, outputGeom );
    
	aRotationConsistency = nAttr.create( "rotationConsistency", "rc", MFnNumericData::kBoolean, 0 );
    nAttr.setStorable(true);
    addAttribute( aRotationConsistency );
    attributeAffects( aRotationConsistency, outputGeom );
    
	aFrechetSum = nAttr.create( "frechetSum", "fs", MFnNumericData::kBoolean, 0 );
    nAttr.setStorable(true);
    addAttribute( aFrechetSum );
    attributeAffects( aFrechetSum, outputGeom );
    
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
