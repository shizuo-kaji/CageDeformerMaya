/**
 * @file blendAff.h
 * @brief blending ffine transformations
 * @section LICENSE The MIT License
 * @section requirements:  Eigen library
 * @version 0.10
 * @date  Oct. 2016
 * @author Shizuo KAJI
 */

#pragma once

#include <map>
#include <Eigen/Sparse>
#include "affinelib.h"
#include "deformerConst.h"

using namespace Eigen;
using namespace AffineLib;

class BlendAff {
public:
    std::vector<Matrix3d> logR,R,logS,S,logGL;
    std::vector<Matrix4d> logSE,SE,logAff,Aff;
    std::vector<Vector3d> L;
    std::vector<Vector4d> quat;
    std::vector<Vector3d> centre;
    
    int num;
    bool rotationConsistency;
    
    BlendAff(int n=0){
        num = n;
        Aff.resize(num);
        centre.resize(num);
        L.resize(num);
    };
    
    void setNum(int n){
        num = n;
        Aff.resize(num);
        centre.resize(num);
        L.resize(num);
    }
    void parametrise(int mode);
    void clearRotation();
};


void BlendAff::parametrise(int mode){
    if(mode == BM_SRL || mode == BM_SSE || mode == BM_SQL){
        R.resize(num); logS.resize(num); S.resize(num);
        for(int i=0;i<num;i++){
            parametriseGL(Aff[i].block(0,0,3,3), logS[i] ,R[i]);
            L[i] = transPart(Aff[i]);
        }
        switch(mode){
            case BM_SRL:
                if(rotationConsistency){
                    logR.resize(num);
                }else{
                    logR.assign(num, Matrix3d::Zero().eval());
                }
                for(int i=0;i<num;i++){
                    logR[i]=logSOc(R[i], logR[i]);
                    S[i]=expSym(logS[i]);
                }
                break;
            case BM_SSE:
                if(rotationConsistency){
                    logSE.resize(num);
                }else{
                    logSE.assign(num, Matrix4d::Zero().eval());
                }
                SE.resize(num);
                for(int i=0;i<num;i++){
                    SE[i]=pad(R[i], L[i]);
                    logSE[i]=logSEc(SE[i], logSE[i]);
                }
                break;
            case BM_SQL:
                quat.resize(num);
                for(int i=0;i<num;i++){
                    Quaternion<double> Q(R[i].transpose());
                    quat[i] << Q.x(), Q.y(), Q.z(), Q.w();
                    S[i]=expSym(logS[i]);
                }
                break;
        }
    }else if(mode == BM_LOG3){
        logGL.resize(num);
        for(int i=0;i<num;i++){
            logGL[i] = Aff[i].block(0,0,3,3).log();
            L[i] = transPart(Aff[i]);
        }
    }else if(mode == BM_LOG4){
        logAff.resize(num);
        for(int i=0;i<num;i++){
            logAff[i] = Aff[i].log();
        }
    }
}


