/*
 *  memory.c
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include "io.h"
#include "memory.h"

double *Calloc1D(double *array, size_t ires) {
	/* Allocate memeory for 1D array of doubles, entries initialized to zero */ 
	if((array = (double *)calloc(ires, sizeof(double))) == NULL)
		FatalError("Calloc failed in function Calloc1D");
	return array;
}

double **Calloc2D(double **array, size_t ires, size_t jres) {
	/* Allocate memeory for 2D array of doubles, entries initialized to zero */
	int i;
	if((array = (double **)calloc(ires, sizeof(double *))) == NULL)
		FatalError("Calloc failed in function Calloc2D");
	for(i = 0; i < ires; i++)
		if((array[i] = (double *)calloc(jres, sizeof(double))) == NULL)
			FatalError("Calloc failed in function Calloc2D");
	return array;
}

double ***Calloc3D(double ***array, size_t ires, size_t jres, size_t kres) {
	/* Allocate memeory for 3D array of doubles, entries initialized to zero */
	int i, j;
	if((array = (double ***)calloc(ires, sizeof(double **))) == NULL)
		FatalError("Calloc failed in function Calloc3D");
	for(i = 0; i < ires; i++) {
		if((array[i] = (double **)calloc(jres, sizeof(double *))) == NULL)
			FatalError("Calloc failed in function Calloc3D");
		for (j = 0; j < jres; j++) 
			if ((array[i][j] = (double *)calloc(kres, sizeof(double))) == NULL) 
				FatalError("Calloc failed in function Calloc3D");
	}
	return array;	
}

double ****Calloc4D(double ****array, size_t tres, size_t ires, size_t jres, size_t kres) {
	/* Allocate memeory for 4D array of doubles, entries initialized to zero */
	int t, i, j;
	if((array = (double ****)calloc(tres, sizeof(double ***))) == NULL)
		FatalError("Calloc failed in function Calloc4D");
	for(t = 0; t < tres; t++) {
		if((array[t] = (double ***)calloc(ires, sizeof(double **))) == NULL)
			FatalError("Calloc failed in function Calloc4D");
		for (i = 0; i < ires; i++) {
			if ((array[t][i] = (double **)calloc(jres, sizeof(double *))) == NULL) 
				FatalError("Calloc failed in function Calloc4D");
			for (j = 0; j < jres; j++)
				if ((array[t][i][j] = (double *)calloc(kres, sizeof(double))) == NULL) 
					FatalError("Calloc failed in function Calloc4D");
		}
	}	
	return array;	
}

void Free2D(double **array, size_t ires) {
	
	int i;
	
	for(i = 0; i < ires; i++) { 
		free(array[i]);
		array[i] = NULL;
	} 
	
	free(array);
	array = NULL;
	
}

void Free3D(double ***array, size_t ires, size_t jres) {
	
	int i, j;
	
	for(i = 0; i < ires; i++) { 
		for(j = 0; j < jres; j++) {
			free(array[i][j]);
			array[i][j] = NULL;
		}
		free(array[i]);
		array[i] = NULL;
	} 
	
	free(array);
	array = NULL;
	
}	

void Free4D(double ****array, size_t tres, size_t ires, size_t jres) {
	
	int tt, i, j;
	
	for(tt = 0; tt < tres; tt++) {
		for(i = 0; i < ires; i++) { 
			for(j = 0; j < jres; j++) {
				free(array[tt][i][j]);
				array[tt][i][j] = NULL;
			}
			free(array[tt][i]);
			array[tt][i] = NULL;
		} 
		free(array[tt]);
		array[tt] = NULL;
	}
	free(array);
	array = NULL;
	
}	




