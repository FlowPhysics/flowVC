/*
 *  tracers.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */



#ifndef INC_TRACERS_H
#define INC_TRACERS_H

void GenerateStaggeredRelease(void);
void LoadReleasePoints(double *voxdim);
void ReadInTraceMesh(void);
void GenerateTracerMesh(void);
void CreateNewTraceLaunch(int ss);
void ReadInTraceLaunch(int ss);
void WriteOutTraceLaunch(int ss);
void OutputTracers(void);
void FreeTracerData(void);

#endif