/*
 *  parameters.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#ifndef INC_PARAMETERS_H
#define INC_PARAMETERS_H

void ReadInParameters(int argc, const char *argv[]);
void ReadInNextValue(FILE *Parameters_InFileID, void *pt, char type);
void SetDerivedParameters(void);
void CheckParameters(void);

#endif