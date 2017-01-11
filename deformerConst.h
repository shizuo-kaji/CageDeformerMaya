//
//  deformerConst.h
//

#pragma once

#ifndef ProbeDeformer_deformerConst_h
#define ProbeDeformer_deformerConst_h

// parametrisation mode
#define BM_SRL 0    // sym^+ so(3) R^3
#define BM_SSE 1    // sym^+ se(3)
#define BM_LOG3 3   // gl(3) R^3
#define BM_LOG4 4   // aff(3)
#define BM_SQL 5    // Sym quat R^3
#define BM_SlRL 6    // Sym so(3) R^3
#define BM_AFF 10   // Aff(3)
#define BM_OFF -1

// weight normalisation mode
#define NM_NONE 0
#define NM_LINEAR 1
#define NM_SOFTMAX 2

// weight mode
#define WM_INV_DISTANCE 0
#define WM_CUTOFF_DISTANCE 1
#define WM_DRAW 2
#define WM_HARMONIC 16
#define WM_HARMONIC_ARAP 16
#define WM_HARMONIC_COTAN 17
#define WM_HARMONIC_TRANS 18

// tetrahedra construction mode
#define TM_FACE 0
#define TM_EDGE 1
#define TM_VERTEX 2
#define TM_VFACE 3   // for each vertex and an adjacent face, make a tet by adding the face normal to the vertex

// cage mode
#define CM_MVC 8
#define CM_MLS 16
#define CM_MLS_AFF 16
#define CM_MLS_SIM 17
#define CM_MLS_RIGID 18

// MLS mode
#define MLS_OFF -1
#define MLS_AFF 0
#define MLS_RIGID 1
#define MLS_SIM 2
#define MLS_POSITIVE_DET 4


// constraint mode
#define CONSTRAINT_NEIGHBOUR 0
#define CONSTRAINT_CLOSEST 1

// stiffness mode
#define SM_NONE 0
#define SM_PAINT 1
#define SM_LEARN 2

// visualisation mode
#define VM_OFF 0
#define VM_ENERGY 1
#define VM_EFFECT 2
#define VM_CONSTRAINT 3
#define VM_STIFFNESS 4

// error codes
#define ERROR_ARAP_PRECOMPUTE 1
#define INCOMPATIBLE_MESH 2
#define ERROR_ATTR 4
#define ERROR_MLS_SINGULAR 8


#endif
