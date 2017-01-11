#pragma once

#include <utility>
#include <vector>
#include <Eigen/Core>

#include "deformerConst.h"

using namespace Eigen;

typedef SparseMatrix<double> SpMat;
typedef Triplet<double> T;

class Distance {
    
public:
    std::vector< std::vector<double> > distPts, distTet;   // [i,j]-entry is the distance to the i-th handle to j-th element
    std::vector<int> closestPts, closestTet;               // [i]-entry is the index of the closest element to i-th handle
    int nHdl, nPts, nTet;
    Distance(){};
    Distance(int _nHandle, int _nPts, int _nTet) {
        setNum(_nHandle, _nPts, _nTet);
    };
    void setNum(int _nHandle, int _nPts, int _nTet);
    void findClosestPts();
    void findClosestTet();
    void computeCageDistPts(short cageMode, const std::vector<Vector3d>& pts, const std::vector<Vector3d>& cagePts, const std::vector<int>& cageTetList);
    void computeCageDistTet(short cageMode, const std::vector<Vector3d>& tetCenter, const std::vector<Vector3d>& cagePts, const std::vector<int>& cageTetList);
    void computeDistPts(const std::vector<Vector3d>& pts, const std::vector<Vector3d>& hdlPts);
    void computeDistTet(const std::vector<Vector3d>& tetCenter, const std::vector<Vector3d>& hdlPts);
    double distPtLin(Vector3d p,Vector3d a,Vector3d b);
    double distPtTri(Vector3d p,Vector3d a,Vector3d b,Vector3d c);
    void MVC(const std::vector<Vector3d>& pts, const std::vector<Vector3d>& cagePts,
                       const std::vector<int>& cageFaceList, std::vector< std::vector<double> >& w);
    void normaliseWeight(short mode, std::vector<double>& w);
};

// initialise
void Distance::setNum(int _nHandle, int _nPts, int _nTet){
    nHdl = _nHandle;
    nPts = _nPts;
    nTet = _nTet;
    distPts.resize(nHdl);
    distTet.resize(nHdl);
    closestPts.resize(nHdl);
    closestTet.resize(nHdl);
    for(int i=0;i<nHdl;i++){
        distPts[i].resize(nPts);
        distTet[i].resize(nTet);
    }
}

// distance between probe handles and mesh pts
void Distance::computeDistPts(const std::vector<Vector3d>& pts, const std::vector<Vector3d>& hdlPts){
    for(int i=0;i<nHdl;i++){
        for(int j=0;j<nPts;j++){
            distPts[i][j] = (pts[j]-hdlPts[i]).norm();
        }
    }
}

// distance between probe handles and mesh tet
void Distance::computeDistTet(const std::vector<Vector3d>& tetCenter, const std::vector<Vector3d>& hdlPts){
    for(int i=0;i<nHdl;i++){
        for(int j=0;j<nTet;j++){
            distTet[i][j] = (tetCenter[j]-hdlPts[i]).norm();
        }
    }
}

// distance between cage and mesh pts
void Distance::computeCageDistPts(short cageMode, const std::vector<Vector3d>& pts, const std::vector<Vector3d>& cagePts, const std::vector<int>& cageTetList){
    switch (cageMode){
        case TM_FACE:
        {
            for(int j=0;j<nPts;j++){
                for(int i=0;i<nHdl;i++){
                    Vector3d a=cagePts[cageTetList[4*i]];
                    Vector3d b=cagePts[cageTetList[4*i+1]];
                    Vector3d c=cagePts[cageTetList[4*i+2]];
                    distPts[i][j] = distPtTri(pts[j], a,b,c);
                }
            }
            break;
        }
        case TM_EDGE:
        {
            for(int j=0;j<nPts;j++){
                for(int i=0;i<nHdl;i++){
                    Vector3d a=cagePts[cageTetList[4*i]];
                    Vector3d b=cagePts[cageTetList[4*i+1]];
                    distPts[i][j] = distPtLin(pts[j], a,b);
                }
            }
            break;
        }
        case TM_VERTEX:
        case TM_VFACE:
        {
            for(int j=0;j<nPts;j++){
                for(int i=0;i<nHdl;i++){
                    distPts[i][j] = (pts[j]-cagePts[cageTetList[4*i]]).norm();
                }
            }
            break;
        }
        case CM_MLS_AFF:
        case CM_MLS_SIM:
        case CM_MLS_RIGID:
        {
            for(int j=0;j<nPts;j++){
                for(int i=0;i<nHdl;i++){
                    distPts[i][j] = (pts[j]-cagePts[i]).norm();
                }
            }
            break;
        }
    }
}

// distance between cage and mesh tet
void Distance::computeCageDistTet(short cageMode, const std::vector<Vector3d>& tetCenter, const std::vector<Vector3d>& cagePts, const std::vector<int>& cageTetList){
    switch (cageMode){
        case TM_FACE:
        {
            for(int j=0;j<nTet;j++){
                for(int i=0;i<nHdl;i++){
                    Vector3d a=cagePts[cageTetList[4*i]];
                    Vector3d b=cagePts[cageTetList[4*i+1]];
                    Vector3d c=cagePts[cageTetList[4*i+2]];
                    distTet[i][j] = distPtTri(tetCenter[j], a,b,c);
                }
            }
            break;
        }
        case TM_EDGE:
        {
            for(int j=0;j<nTet;j++){
                for(int i=0;i<nHdl;i++){
                    Vector3d a=cagePts[cageTetList[4*i]];
                    Vector3d b=cagePts[cageTetList[4*i+1]];
                    distTet[i][j] = distPtLin(tetCenter[j], a,b);
                }
            }
            break;
        }
        case TM_VERTEX:
        case TM_VFACE:
        {
            for(int j=0;j<nTet;j++){
                for(int i=0;i<nHdl;i++){
                    distTet[i][j] = (tetCenter[j]-cagePts[cageTetList[4*i]]).norm();
                }
            }
            break;
        }
        case CM_MLS_AFF:
        case CM_MLS_SIM:
        case CM_MLS_RIGID:
        {
            for(int j=0;j<nTet;j++){
                for(int i=0;i<nHdl;i++){
                    distTet[i][j] = (tetCenter[j]-cagePts[i]).norm();
                }
            }
            break;
        }
    }
}


// find closest point on mesh from each handle
void Distance::findClosestPts(){
    for(int i=0;i<nHdl;i++){
        closestPts[i] = 0;
        double min_d = HUGE_VAL;
        for(int j=0;j<nPts;j++){
            if( distPts[i][j] < min_d){
                min_d = distPts[i][j];
                closestPts[i] = j;
            }
        }
    }
}
// find closest tet on mesh from each handle
void Distance::findClosestTet(){
    for(int i=0;i<nHdl;i++){
        closestTet[i] = 0;
        double min_d = HUGE_VAL;
        for(int j=0;j<nTet;j++){
            if( distTet[i][j] < min_d){
                min_d = distTet[i][j];
                closestTet[i] = j;
            }
        }
    }
}


// compute distance between a line segment (ab) and a point p
double Distance::distPtLin(Vector3d p,Vector3d a,Vector3d b){
    double t= (a-b).dot(p-b)/(a-b).squaredNorm();
    if(t>1){
        return (a-p).norm();
    }else if(t<0){
        return (b-p).norm();
    }else{
        return (t*(a-b)-(p-b)).norm();
    }
}

// compute distance between a triangle (abc) and a point p
double Distance::distPtTri(Vector3d p, Vector3d a, Vector3d b, Vector3d c){
    /// if p is in the outer half-space, it returns HUGE_VAL
    double s[4];
    Vector3d n=(b-a).cross(c-a);
    if(n.squaredNorm()<EPSILON){
        return (p-a).norm();
    }
    double k=n.dot(a-p);
    if(k<0) return HUGE_VAL;
    s[0]=distPtLin(p,a,b);
    s[1]=distPtLin(p,b,c);
    s[2]=distPtLin(p,c,a);
    Matrix3d A;
    A << b(0)-a(0), c(0)-a(0), n(0)-a(0),
    b(1)-a(1), c(1)-a(1), n(1)-a(1),
    b(2)-a(2), c(2)-a(2), n(2)-a(2);
    Vector3d v = A.inverse()*(p-a);  // barycentric coordinate of p
    if(v(0)>0 && v(1)>0 && v(0)+v(1)<1){
        s[3]=k;
    }else{
        s[3] = HUGE_VAL;
    }
    return min(min(min(s[0],s[1]),s[2]),s[3]);
}

// mean value coordinate
void Distance::MVC(const std::vector<Vector3d>& pts, const std::vector<Vector3d>& cagePts,
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

// normalise weights
void Distance::normaliseWeight(short mode, std::vector<double>& w){
    if(mode == NM_NONE || mode == NM_LINEAR){
        double sum = std::accumulate(w.begin(), w.end(), 0.0);
        if ((sum > 1 || mode == NM_LINEAR) && sum != 0.0){
            for (int i = 0; i < w.size(); i++){
                w[i] /= sum;
            }
        }
    }else if(mode == NM_SOFTMAX){
        double sum = 0.0;
        for (int i = 0; i < w.size(); i++){
            sum += exp(w[i]);
        }
        for (int i = 0; i < w.size(); i++){
            w[i] = exp(w[i])/sum;
        }
    }
}
