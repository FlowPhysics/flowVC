/*
 *  exposuretime.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#ifndef INC_EXPOSURETIME_H
#define INC_EXPOSURETIME_H

void LoadCETMeshData(void);
void ComputeExposureTime(const double *X0, const double *X1, const int X1ElementIndex, const double stepsize);
void OutputCET(void);

#endif