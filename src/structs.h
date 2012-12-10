/*
 *  structs.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#ifndef INC_STRUCTS_H
#define INC_STRUCTS_H

#include <stdio.h>
#include "macros.h"

typedef struct _Element {
	int Nodes[4];
	int NeighborIndex[4];
} Element;

typedef struct _LagrangianPoint { 
	double X[3];
	double V[3];
    double vorticity[3];
	int    ElementIndex;
	int    AuxElementIndex;
	int    LeftDomain;
	double LeftDomainTime;
	double Scalar;
} LagrangianPoint;

typedef struct _FTLEPoint {
	LagrangianPoint Pt;
	int    HaveFTLE;
	double FTLEwT;
	double FTLEwoT;
} FTLEPoint;

typedef struct _Launch {
	double StartTime;
	double StopTime;
	int    Status;
} Launch;

typedef struct _CartMesh {
	double XMin;
	double XMax;
	int    XRes;
	double XDelta;
	
	double YMin;
	double YMax;
	int    YRes;
	double YDelta;
	
	double ZMin;
	double ZMax;
	int    ZRes;
	double ZDelta;
} CartMesh;

typedef struct _FileData {
	char  FilePath[LONGSTRING];
	FILE *FileID;
} FileData;

typedef struct _ReleaseLocation { 
	Launch slide;
	LagrangianPoint pt;
	struct _ReleaseLocation *next;
} ReleaseLocation;

typedef struct _CETE {
	double CETsum;
	int Encounters;
} CETE;

#endif