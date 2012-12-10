/*
 *  memory.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#ifndef INC_MEMORY_H
#define INC_MEMORY_H

double *Calloc1D(double *array, size_t ires); 
double **Calloc2D(double **array, size_t ires, size_t jres);
double ***Calloc3D(double ***array, size_t ires, size_t jres, size_t kres);
double ****Calloc4D(double ****array, size_t tres, size_t ires, size_t jres, size_t kres);
void Free2D(double **array, size_t ires);
void Free3D(double ***array, size_t ires, size_t jres);
void Free4D(double ****array, size_t tres, size_t ires, size_t jres);


#endif