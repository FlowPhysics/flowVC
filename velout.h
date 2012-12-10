/*
 *  velout.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#ifndef INC_VELOUT_H
#define INC_VELOUT_H

void GenerateVelOutMesh(void);
void ReadInVelOutMesh(void);
void OutputVelOut(double t, int slide);
void FreeVelOutData(void);

#endif