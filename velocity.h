/*
 *  velocity.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#ifndef INC_VELOCITY_H
#define INC_VELOCITY_H

void AllocateVelFieldData(void);
void FreeVelFieldData(void);
void LoadCartVelDataFrame(int frame);
void LoadUnstructVelDataFrame(int frame);
void GetVelocity(double t,  LagrangianPoint *pt, double *dXdt);
void GetVelocity_Cartesian(double t,  LagrangianPoint *pt, double *dXdt);
void GetVelocity_Unstructured(const double time,  LagrangianPoint *pt, double *dXdt);
int GetLocalIndex_SurfaceMesh(int global_index); 
void GetVelocity_Analytic(double t,  LagrangianPoint *pt, double *dXdt);
int InsideBoundaryElement(LagrangianPoint *pt);
int TestOutsideCartVelDomain(double X[3]);
int  TestOutsideDomain(double point[3]);

#endif
