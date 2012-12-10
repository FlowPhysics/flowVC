/*
 *  ftle.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */


#ifndef INC_FTLE_H
#define INC_FTLE_H

#include "structs.h"

void InitializeFTLEArray(void);
void ReadInFTLELaunch(int ss);
void WriteOutFTLELaunch(int ss);
void GetFTLEForPointEarly(int ss, int i, int j, int k, double t1, double te);
void GetFTLEForPoint(int i,int j, int k, double IntTime);
LagrangianPoint Advect_FTLEPoint(int i, int j, int k, double t1, double t2);
void UpdateFTLELocations(void);
void OutputFTLE(int ss, double IntT);
void FreeFTLEData(void);

#endif