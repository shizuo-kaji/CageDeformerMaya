/**
 * @file cageDeformer.cpp
 * @brief Cage Deformer plugin for Maya
 * @section LICENSE The MIT License
 * @section requirements:  Eigen 3:  http://eigen.tuxfamily.org/
 * @section Autodesk Maya: http://www.autodesk.com/products/autodesk-maya/overview
 * @section (included) AffineLib: https://github.com/shizuo-kaji/AffineLib
 * @version 0.15
 * @date  3/Dec/2012
 * @author Shizuo KAJI
 */

#include "StdAfx.h"
#include "cageDeformer.h"

using namespace Eigen;
using namespace AffineLib;

MTypeId CageDeformerNode::id( 0x00000200 );
MString CageDeformerNode::nodeName("cage");
MObject CageDeformerNode::aCageMesh;
MObject CageDeformerNode::aCageMode;
MObject CageDeformerNode::aBlendMode;
MObject CageDeformerNode::aRotationConsistency;
MObject CageDeformerNode::aFrechetSum;

// compute distance between a line segment and a point
double distPtLin(Vector3d p,Vector3d a,Vector3d b){
    double t= (a-b).dot(p-b)/(a-b).squaredNorm();
    if(t>1){
        return (a-p).squaredNorm();
    }else if(t<0){
        return (b-p).squaredNorm();
    }else{
        return (t*(a-b)-(p-b)).squaredNorm();
    }
}

// compute distance between a triangle and a point
double distPtTri(MPoint pt, Matrix4d m){
    double s[4];
    Vector3d a,b,c,n,p;
    a << m(0,0), m(0,1), m(0,2);
    b << m(1,0), m(1,1), m(1,2);
    c << m(2,0), m(2,1), m(2,2);
    n << m(3,0), m(3,1), m(3,2);
    p << pt.x, pt.y, pt.z;
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

// check linear independency
double isDegenerate(MPoint a,MPoint b,MPoint c,MPoint d){
    return ((b-a)^(c-a)) * (d-a);
}

// main
void* CageDeformerNode::creator() { return new CageDeformerNode; }
 
MStatus CageDeformerNode::deform( MDataBlock& data, MItGeometry& itGeo, const MMatrix &localToWorldMatrix, unsigned int mIndex )
{
    MStatus status;
    MThreadUtils::syncNumOpenMPThreads();    // for OpenMP
    // load cage mesh
    MObject oCageMesh = data.inputValue( aCageMesh ).asMesh();
    short blendMode = data.inputValue(aBlendMode).asShort();
	bool rotationCosistency = data.inputValue( aRotationConsistency ).asBool();
	bool frechetSum = data.inputValue( aFrechetSum ).asBool();
    if ( oCageMesh.isNull() || blendMode == 99)
        return MS::kSuccess;
    // load attributes
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
    // load target mesh pts
    MPointArray pts;
    itGeo.allPositions(pts);
    int numPts = pts.length();
    for(int j=0; j<numPts; j++ )
        pts[j] *= localToWorldMatrix;
    // when cage mode is changed
    if(newCageMode != cageMode)
    {
        cageMode = newCageMode;
        MFnMesh fnInitCageMesh( initCageMesh, &status );
	    std::vector<double> tetWeight;
        if(cageMode == 10 || cageMode == 11)  // face mode
        {
			if(cageMode == 10){
                MIntArray count;
                fnInitCageMesh.getTriangles( count, triangles );
                tetWeight.resize(triangles.length()/3, 1.0);
			}else if(cageMode ==11){  // trianglate faces with more than 3 edges in a symmetric way
				triangles.clear();
				MItMeshPolygon iter(initCageMesh);
				MIntArray tmp;
                MVector normal;
				tetWeight.reserve(4*iter.count());
                int l;
				for(int i=0; ! iter.isDone(); i++){
					iter.getVertices(tmp);
					l=tmp.length();
					if(l==3){
						tetWeight.push_back(1.0);
						triangles.append(tmp[0]);
						triangles.append(tmp[1]);
						triangles.append(tmp[2]);
					}else{
						for(int j=0;j<l;j++){
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
            numTet=triangles.length()/3;
            std::vector<Matrix4d> initMatrix(numTet);
            initInvMatrix.resize(numTet);
            tetMatrix(initCagePoints, triangles, cageMode, initMatrix);
            for(int i=0; i<numTet; i++)
                initInvMatrix[i]=initMatrix[i].inverse();
            // compute weight for non-ARAP mode
            w.resize(numPts);
            std::vector< std::vector<double> > idist(numPts);
            for(int j=0;j<numPts;j++){
                idist[j].resize(numTet);
                w[j].resize(numTet);
                double sidist = 0.0;
                for(int i=0;i<numTet;i++){
                    idist[j][i] = tetWeight[i]/distPtTri(pts[j],initMatrix[i]);
                    sidist += idist[j][i];
                }
                assert(sidist>0.0);
                for(int i=0;i<numTet;i++)
                    w[j][i] = idist[j][i] /sidist;
            } // face mode end
        }else if(cageMode == 0 || cageMode == 1 ){ // vertex mode
            triangles.clear();
            std::vector<int> tetCount(initCagePoints.length());
            MItMeshVertex iter(initCageMesh);
            for(int j=0; ! iter.isDone(); j++){
                MIntArray v;
                iter.getConnectedVertices(v);
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
			numTet=triangles.length()/4;
            std::vector<Matrix4d> initMatrix(numTet);
            initInvMatrix.resize(numTet);
            tetMatrix(initCagePoints, triangles, cageMode, initMatrix);
            for(int i=0; i<numTet; i++)
                initInvMatrix[i]=initMatrix[i].inverse();
            w.resize(numPts);
            std::vector< std::vector<double> > idist(numPts);
			tetWeight.resize(numTet);
			for(int i=0;i<numTet;i++)
				tetWeight[i]=1.0/(double)tetCount[triangles[4*i]];
            for(int j=0;j<numPts;j++){
                idist[j].resize(numTet);
                w[j].resize(numTet);
                double sidist = 0.0;
                Vector3d p(pts[j].x,pts[j].y,pts[j].z);
                for(int i=0;i<numTet;i++){
                    Vector3d c(initCagePoints[triangles[4*i]].x,initCagePoints[triangles[4*i]].y,initCagePoints[triangles[4*i]].z);
                    idist[j][i] = tetWeight[i] / ((p-c).squaredNorm());
                    sidist += idist[j][i];
                }
                assert(sidist>0.0);
                for(int i=0;i<numTet;i++)
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
			numTet=triangles.length()/3;
            std::vector<Matrix4d> initMatrix(numTet);
            initInvMatrix.resize(numTet);
            tetMatrix(initCagePoints, triangles, cageMode, initMatrix);
            for(int i=0; i<numTet; i++)
            initInvMatrix[i]=initMatrix[i].inverse();
            w.resize(numPts);
            std::vector< std::vector<double> > idist(numPts);
			tetWeight.resize(numTet);
			for(int i=0;i<numTet;i++)
            tetWeight[i]=1.0/(double)tetCount[triangles[3*i]];
            for(int j=0;j<numPts;j++){
                idist[j].resize(numTet);
                w[j].resize(numTet);
                double sidist = 0.0;
                Vector3d p(pts[j].x,pts[j].y,pts[j].z);
                for(int i=0;i<numTet;i++){
                    Vector3d c(initCagePoints[triangles[3*i]].x,initCagePoints[triangles[3*i]].y,initCagePoints[triangles[3*i]].z);
                    idist[j][i] = tetWeight[i] / ((p-c).squaredNorm());
                    sidist += idist[j][i];
                }
                assert(sidist>0.0);
                for(int i=0;i<numTet;i++)
                w[j][i] = idist[j][i] /sidist;
            }
        }else if(cageMode == 20){ // MVC mode
            MIntArray count;
            fnInitCageMesh.getTriangles( count, triangles );
            numTet=triangles.length()/3;
			int numCagePoints=initCagePoints.length();
			for(int j=0; j<numPts; j++ )
				w[j].resize(numCagePoints);
			MVC(pts, initCagePoints, triangles, w);
		}
    }
    
    // transform vertices
	if(cageMode<20){
        // clear previous rotation
        if( ! rotationCosistency || numTet != prevNs.size()){
            prevThetas.clear();
            prevThetas.resize(numTet, 0.0);
            prevNs.clear();
            prevNs.resize(numTet, Vector3d::Zero());
        }
        //  find affine transformations for tetrahedra
        std::vector<Matrix4d> cageMatrix(numTet), SE(numTet), logSE(numTet),logAff(numTet),aff(numTet);
        std::vector<Matrix3d> logR(numTet),R(numTet),logS(numTet),logGL(numTet);
        std::vector<Vector3d> L(numTet);
        std::vector<Vector4d> quat(numTet);
        tetMatrix(cagePoints, triangles, cageMode, cageMatrix);
        for(int i=0; i<numTet; i++)
            aff[i]=initInvMatrix[i]*cageMatrix[i];
        // compute parametrisation
        if(blendMode == 0 || blendMode == 1 || blendMode == 5)  // polarexp or quaternion
        {
            for(int i=0;i<numTet;i++){
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
            for(int i=0;i<numTet;i++){
                logGL[i] = aff[i].block(0,0,3,3).log();
                L[i] = transPart(aff[i]);
            }
        }else if(blendMode == 3){   // logmatrix4
            for(int i=0;i<numTet;i++){
                logAff[i] = aff[i].log();
            }
        }
        // compute blended matrix
#pragma omp parallel for
		for(int j=0; j<numPts; j++ ){
            Matrix4d mat=Matrix4d::Zero();
            if(blendMode==0){
                Matrix3d RR=Matrix3d::Zero();
                Matrix3d SS=Matrix3d::Zero();
                Vector3d l=Vector3d::Zero();
                for(int i=0; i<numTet; i++){
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
                mat = affine(SS*RR, l);
            }else if(blendMode==1){    // rigid transformation
                Matrix4d EE=Matrix4d::Zero();
                Matrix3d SS=Matrix3d::Zero();
                for(int i=0; i<numTet; i++){
                    EE +=  w[j][i] * logSE[i];
                    SS +=  w[j][i] * logS[i];
                }
                if(frechetSum){
                    EE = frechetSE(SE, w[j]);
                }else{
                    EE = expSE(EE);
                }
                mat = affine(expSym(SS),Vector3d::Zero())*EE;
            }else if(blendMode == 2){    //logmatrix3
                Matrix3d G=Matrix3d::Zero();
                Vector3d l=Vector3d::Zero();
                for(int i=0; i<numTet; i++){
                    G +=  w[j][i] * logGL[i];
                    l += w[j][i] * L[i];
                }
                mat = affine(G.exp(), l);
            }else if(blendMode == 3){   // logmatrix4
                Matrix4d A=Matrix4d::Zero();
                for(int i=0; i<numTet; i++)
                    A +=  w[j][i] * logAff[i];
                mat = A.exp();
            }else if(blendMode == 5){ // quaternion
                Vector4d q=Vector4d::Zero();
                Matrix3d SS=Matrix3d::Zero();
                Vector3d l=Vector3d::Zero();
                for(int i=0; i<numTet; i++){
                    q += w[j][i] * quat[i];
                    SS += w[j][i] * logS[i];
                    l += w[j][i] * L[i];
                }
                SS = expSym(SS);
                Quaternion<double> Q(q);
                Matrix3d RR = Q.matrix().transpose();
                mat = affine(SS*RR, l);
            }else if(blendMode==10){
                for(int i=0; i<numTet; i++){
                    mat += w[j][i] * aff[i];
                }
            }
            // apply matrix
            MMatrix m;
            for(int i=0; i<4; i++)
                for(int k=0; k<4; k++)
                    m(i,k)=mat(i,k);
            pts[j] *= m*localToWorldMatrix.inverse();
        }
	}else if(cageMode == 20){      // MVC mode
		int numCagePoints = cagePoints.length();
		for(int j=0; j<numPts; j++ ){
			pts[j] = MPoint::origin;
			for(int i=0;i<numCagePoints;i++)
				pts[j] +=  cagePoints[i] * w[j][i];
		}
	}
    itGeo.setAllPositions(pts);
 
    return MS::kSuccess;
}

// prepare tetrahedra
void CageDeformerNode::tetMatrix(const MPointArray& p, const MIntArray& triangles, short cageMode, std::vector<Matrix4d>& m)
{
    int numTet;
    if(cageMode == 10 || cageMode == 11){
        MVector u, v, q ;
        numTet=triangles.length()/3;
        for(int i=0;i<numTet;i++){
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
        numTet=triangles.length()/4;
        for(int i=0;i<numTet;i++){
            m[i] << p[triangles[4*i]].x, p[triangles[4*i]].y, p[triangles[4*i]].z, 1,
            p[triangles[4*i+1]].x, p[triangles[4*i+1]].y, p[triangles[4*i+1]].z, 1,
            p[triangles[4*i+2]].x, p[triangles[4*i+2]].y, p[triangles[4*i+2]].z, 1,
            p[triangles[4*i+3]].x, p[triangles[4*i+3]].y, p[triangles[4*i+3]].z, 1;
        }
    }else if (cageMode == 1){
        numTet=triangles.length()/4;
        MVector u, v, q;
        for(int i=0;i<numTet;i++){
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
        numTet=triangles.length()/3;
        std::vector<Vector3d> normal(p.length(),Vector3d::Zero());           //averaged normal vector for vertex mode
        for(int i=0;i<numTet;i++){
            MVector q = MVector(p[triangles[3*i+1]]-p[triangles[3*i]]) ^ MVector(p[triangles[3*i+2]]-p[triangles[3*i]]);
            normal[triangles[3*i]] += Vector3d(q[0], q[1], q[2]);
        }
        if (cageMode == 5){
            for(int i=0;i<numTet;i++){
                m[i] << p[triangles[3*i]].x, p[triangles[3*i]].y, p[triangles[3*i]].z, 1,
                p[triangles[3*i+1]].x, p[triangles[3*i+1]].y, p[triangles[3*i+1]].z, 1,
                p[triangles[3*i+2]].x, p[triangles[3*i+2]].y, p[triangles[3*i+2]].z, 1,
                p[triangles[3*i]].x + normal[triangles[3*i]][0], p[triangles[3*i]].y + normal[triangles[3*i]][1], p[triangles[3*i]].z + normal[triangles[3*i]][2],1;
            }
        }else if (cageMode == 6){
            MVector u,v,q;
            for(int i=0;i<numTet;i++){
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


// compute mean value coordinate
void CageDeformerNode::MVC(const MPointArray& pts, const MPointArray& cagePoints, const MIntArray& triangles, std::vector< std::vector<double> >& w)
{
	int numPts=pts.length();
	int num=cagePoints.length();
	int numTriangles=triangles.length()/3;
	std::vector<Vector3d> q(num);
	for(int i=0;i<num;i++)
        q[i] << cagePoints[i].x, cagePoints[i].y, cagePoints[i].z;

#pragma omp parallel for
	for(int j=0; j<numPts; j++ ){
		std::vector<double> mu(num), a(3), b(3);
		Vector3d v;
		v << pts[j].x, pts[j].y, pts[j].z;
		std::vector<Vector3d> e(3),n(3);
		for(int i=0;i<numTriangles;i++){
			for(int k=0;k<3;k++)
				e[k]=(v-q[triangles[3*i+k]]).normalized();
			for(int k=0;k<3;k++)
				n[k]=(e[(k+1)%3].cross(e[(k+2)%3])).normalized();
			for(int k=0;k<3;k++){
				a[k]=n[(k+1)%3].dot(n[(k+2)%3]);
				b[k]=acos(e[(k+1)%3].dot(e[(k+2)%3]));
			}
			for(int k=0;k<3;k++)
                mu[triangles[3*i+k]] -= (b[k]+b[(k+2)%3]*a[(k+1)%3]+b[(k+1)%3]*a[(k+2)%3])/(2.0*e[k].dot(n[k]));
		}
		double smu=0.0;
		for(int i=0;i<num;i++){
			mu[i] /= (v-q[i]).norm();
			smu += mu[i];
		}
		for(int i=0;i<num;i++)
			w[j][i] = mu[i]/smu;
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
    eAttr.addField( "MVC", 20 );
//    eAttr.addField( "init", 99 );
    addAttribute( aCageMode );
    attributeAffects( aCageMode, outputGeom );

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
