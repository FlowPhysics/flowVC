/*
 *  vorticity.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#ifndef INC_VORTICITY_H
#define INC_VORTICITY_H

void AllocateVorticityFieldData(void);
void FreeVorticityData(void);
void LoadVorticityDataFrame(int frame);
void LoadCartVorticityDataFrame(int frame);
void LoadUnstructVorticityDataFrame(int frame);
void SetVorticity(double tc, LagrangianPoint *pt);
void SetVorticity_Cartesian(double time, LagrangianPoint *pt);
void SetVorticity_Unstructured(double time, LagrangianPoint *pt);

#endif