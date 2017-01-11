/**
 * @file MeshMaya.h
 * @brief Mesh manipulation functions for Maya
 * @section LICENSE The MIT License
 * @section requirements:  Eigen library,   Maya
 * @version 0.10
 * @date  Aug. 2014
 * @author Shizuo KAJI
 */

#pragma once


#include "tetrise.h"
#include <set>

using namespace Tetrise;



// get face list
int makeFaceTet(MDataBlock& data, MObject& input, MObject& inputGeom, unsigned int mIndex, const std::vector<Vector3d>& pts,
                std::vector<int>& tetList, std::vector<Matrix4d>& tetMatrix, std::vector<double>& tetWeight){
    // returns total number of pts including ghost ones
    // read mesh data
    int numPts = (int) pts.size();
    MStatus status;
    MArrayDataHandle hInput = data.outputArrayValue( input, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    status = hInput.jumpToElement( mIndex );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MObject oInputGeom = hInput.outputValue().child( inputGeom ).asMesh();
    MFnMesh inputMesh(oInputGeom);
    // face list
    MIntArray count, triangles;
    inputMesh.getTriangles( count, triangles );
    std::vector<int> faceList(triangles.length());
    for(int i=0;i<triangles.length();i++){
        faceList[i]=triangles[i];
    }
    //
    std::vector<vertex> vList;
    std::vector<edge> eList;
    makeTetList(TM_FACE, numPts, faceList, eList, vList, tetList);
    makeTetMatrix(TM_FACE, pts, tetList, faceList, eList, vList, tetMatrix, tetWeight);
    return numPts + (int)tetList.size()/4;
}

// make face list
void makeFaceList(MObject& mesh, std::vector<int>& faceList, std::vector<int>& faceCount,
                  bool isSymmetric=false){
    if( isSymmetric ){
        MItMeshPolygon iter(mesh);
        MIntArray faceVertices;
        faceList.clear();
        for(int i=0; ! iter.isDone(); i++){
            iter.getVertices(faceVertices);
            int count=faceVertices.length();
            if(count==3){
                faceCount.push_back(1);
                faceList.push_back(faceVertices[0]);
                faceList.push_back(faceVertices[1]);
                faceList.push_back(faceVertices[2]);
            }else{
                for(int j=0;j<count;j++){
                    faceCount.push_back(count);
                    faceList.push_back(faceVertices[j]);
                    faceList.push_back(faceVertices[(j+1) % count]);
                    faceList.push_back(faceVertices[(j+2) % count]);
                }
            }
            iter.next();
        }
    }else{
        MFnMesh fnMesh(mesh);
        MIntArray count, triangles;
        fnMesh.getTriangles( count, triangles );
        faceCount.resize(triangles.length()/3,1);
        faceList.resize(triangles.length());
        for(int i=0;i<triangles.length();i++){
            faceList[i]=triangles[i];
        }
    }
}

// vertex list
void makeVertexList(MObject& mesh, std::vector<vertex>& vertexList){
    int numPts = MFnMesh(mesh).numVertices();
    MItMeshPolygon iter(mesh);
    MIntArray faceVertices;
    vertexList.resize(numPts);
    for(int i=0;i<numPts;i++){
        vertexList[i].index = i;
        vertexList[i].connectedTriangles.clear();
    }
    for( ; ! iter.isDone(); iter.next()){
        iter.getVertices(faceVertices);
        int count = (int) faceVertices.length();
        for(int j=0;j<count;j++){
            vertexList[faceVertices[j]].connectedTriangles.push_back(faceVertices[(j+1)%count]);
            vertexList[faceVertices[j]].connectedTriangles.push_back(faceVertices[(j+count-1)%count]);
        }
    }
}

// get mesh data
int getMeshData(MDataBlock& data, MObject& input, MObject& inputGeom, unsigned int mIndex,
                short tetMode, const std::vector<Vector3d>& pts, std::vector<int>& tetList,
                std::vector<int>& faceList, std::vector<edge>& edgeList,
                std::vector<vertex>& vertexList, std::vector<Matrix4d>& tetMat, std::vector<double>& tetWeight){
    // returns total number of pts including ghost ones
    // read mesh data
    int numPts = (int) pts.size();
    MStatus status;
    MArrayDataHandle hInput = data.outputArrayValue( input, &status );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    status = hInput.jumpToElement( mIndex );
    CHECK_MSTATUS_AND_RETURN_IT( status );
    MObject oInputGeom = hInput.outputValue().child( inputGeom ).asMesh();
    std::vector<int> faceCount;
    makeFaceList(oInputGeom, faceList, faceCount);
    makeVertexList(oInputGeom, vertexList);
    makeEdgeList(faceList, edgeList);
    int dim=makeTetList(tetMode, numPts, faceList, edgeList, vertexList, tetList);
    makeTetMatrix(tetMode, pts, tetList, faceList, edgeList, vertexList, tetMat, tetWeight);
    return dim;
}

// visualise vertex assigned values
void visualise(MDataBlock& data, MObject& outputGeom, std::vector<double>& ptsColour){
    // load target mesh output
    MStatus status;
    MArrayDataHandle outputArray = data.outputArrayValue(outputGeom , &status );
    MDataHandle hOutput = outputArray.inputValue(&status);
    MFnMesh outMesh(hOutput.data());
    MColorArray Colours;
    MIntArray Index;
    // set vertex colour
    for(int i=0;i<ptsColour.size();i++){
        MColor colour(MColor::kHSV, 0.0, ptsColour[i], 1.0);
        Colours.append( colour );
        Index.append(i);
    }
    outMesh.setVertexColors(Colours, Index);
}

// read array of matrix attributes and convert them to Eigen matrices
void readMatrixArray(MArrayDataHandle& handle, std::vector<Matrix4d>& m){
    int numPrb=handle.elementCount();
    m.resize(numPrb);
    MMatrix mat;
    for(int i=0;i<numPrb;i++){
        handle.jumpToArrayElement(i);
        mat=handle.inputValue().asMatrix();
        m[i] << mat(0,0), mat(0,1), mat(0,2), mat(0,3),
        mat(1,0), mat(1,1), mat(1,2), mat(1,3),
        mat(2,0), mat(2,1), mat(2,2), mat(2,3),
        mat(3,0), mat(3,1), mat(3,2), mat(3,3);
    }
}

// read array of vector attributes into Eigen vectors
void readVectorArray(MArrayDataHandle& handle, std::vector<Vector3d>& V){
    int num=handle.elementCount();
    V.resize(num);
    MVector v;
    for(int i=0;i<num;i++){
        handle.jumpToArrayElement(i);
        v=handle.inputValue().asVector();
        V[i] << v(0), v(1), v(2);
    }
}


////
void outputAttr(MDataBlock& data, MObject& attribute, std::vector<double>& values){
    MStatus status;
    MArrayDataBuilder builder(attribute, (int)values.size(), &status);
    for(int i=0;i<values.size();i++){
        MDataHandle outHandle = builder.addElement(i);
        outHandle.set(values[i]);
    }
    MArrayDataHandle outputArray = data.outputArrayValue(attribute);
    outputArray.set(builder);
}

////
void deleteAttr(MDataBlock& data, MObject& attribute, std::set<int>& indices){
    MStatus status;
    MArrayDataHandle outputArray = data.outputArrayValue(attribute);
    MArrayDataBuilder builder=outputArray.builder();
    std::set<int>::iterator iter;
    for(iter = indices.begin(); iter != indices.end(); iter++){
        builder.removeElement(*iter);
    }
    outputArray.set(builder);
}


