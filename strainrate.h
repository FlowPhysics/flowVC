/*
 *  strainrate.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#ifndef INC_STRAINRATE_H
#define INC_STRAINRATE_H

void AllocateStrainRateData(void);
void FreeStrainRateData(void);
void LoadStrainRateDataFrame(int frame);
void LoadCartStrainRateDataFrame(int frame);
void LoadUnstructStrainRateDataFrame(int frame);
double GetStrainRate(double tc, LagrangianPoint *pt);
double GetStrainRate_Cartesian(double time, LagrangianPoint *pt);
double GetStrainRate_Unstructured(double time, LagrangianPoint *pt);
void OutputAP(void);
void OutputStrainOut(double t, int slide);

#endif