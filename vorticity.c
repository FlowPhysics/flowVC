//
//  vorticity.c
//  
//
//  Created by Shawn Shadden.
//  Copyright 2010 Illinois Institute of Technology. All rights reserved.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "globals.h"
#include "io.h"
#include "macros.h"
#include "memory.h"
#include "mesh.h"
#include "vorticity.h"
#include "structs.h"
#include "tracers.h"
#include "velocity.h"


void AllocateVorticityFieldData(void) {
	/*
	 *	Allocates memory for vorticity field data for 2 time frames (moving buffer of data)
	 */
	
	printf("\n  Allocating memory for vorticity data...");
    fflush(stdout);
	if (Data_MeshType == CARTESIAN) {
		Vel_CartVorArray_wx = Calloc4D(Vel_CartVorArray_wx, 2, Vel_CartMesh.XRes - 1, Vel_CartMesh.YRes - 1, Vel_CartMesh.ZRes - 1);
		Vel_CartVorArray_wy = Calloc4D(Vel_CartVorArray_wy, 2, Vel_CartMesh.XRes - 1, Vel_CartMesh.YRes - 1, Vel_CartMesh.ZRes - 1);
		Vel_CartVorArray_wz  = Calloc4D(Vel_CartVorArray_wz, 2, Vel_CartMesh.XRes - 1, Vel_CartMesh.YRes - 1, Vel_CartMesh.ZRes - 1);
	}
	else {
		Vel_UnstructVorArray_wx = Calloc2D(Vel_UnstructVorArray_wx, 2, Vel_MeshNumElements);
		Vel_UnstructVorArray_wy = Calloc2D(Vel_UnstructVorArray_wy, 2, Vel_MeshNumElements);
		Vel_UnstructVorArray_wz  = Calloc2D(Vel_UnstructVorArray_wz, 2, Vel_MeshNumElements);
	}
	printf("OK!\n");
    fflush(stdout);
	
}

void FreeVorticityData(void) {
	/*** 
	 Deallocate memory for vorticity field for 2 time frames (moving buffer of data)
	 ***/
	if (Data_MeshType == CARTESIAN) {
		Free4D(Vel_CartVorArray_wx, 2, Vel_CartMesh.XRes - 1, Vel_CartMesh.YRes - 1);
		Free4D(Vel_CartVorArray_wy, 2, Vel_CartMesh.XRes - 1, Vel_CartMesh.YRes - 1);
		Free4D(Vel_CartVorArray_wz, 2, Vel_CartMesh.XRes - 1, Vel_CartMesh.YRes - 1);
	}
	else {
		Free2D(Vel_UnstructVorArray_wx, 2);
		Free2D(Vel_UnstructVorArray_wy, 2);
		Free2D(Vel_UnstructVorArray_wz, 2);		
	}
}


void LoadVorticityDataFrame(int frame) {
    if(Data_MeshType == CARTESIAN)  
        LoadVorticityDataFrame(frame);
    else if(Data_MeshType == UNSTRUCTURED)
        LoadUnstructVorticityDataFrame(frame);
    else
        FatalError("Unknown value for Data_MeshType in LoadVorticityDataFrame()");
}


void LoadCartVorticityDataFrame(int frame) {
	
	int i, j, k, ModVal, slot1, slot2;
	double ***tempptr = NULL, timestamp;
	char Data_BinFilePath[LONGSTRING];
	FILE *Data_BinFileID;
	
	if(Int_TimeDirection > 0) {
		slot1 = 0;
		slot2 = 1;
	}    
	else {
		slot1 = 1; 
		slot2 = 0;
	}
	
	ModVal = fmod(frame, Data_TRes - 1); 
	if(ModVal < 0) {
		ModVal += Data_TRes - 1;
	}
	
	if(Vel_CartVorArray_wx == NULL || Vel_CartVorArray_wy == NULL || Vel_CartVorArray_wz == NULL)
		AllocateVorticityFieldData(); 
	
	/* Allocate memory for tempptr, used as placeholder when swapping pointers to data as window of loaded velocity data changes */
	if((tempptr = (double ***)malloc((Vel_CartMesh.XRes - 1) * sizeof(double **))) == NULL)
		FatalError("Malloc failed for tempptr");
	for(i = 0; i < Vel_CartMesh.XRes - 1; i++) 
		if((tempptr[i] = (double **)malloc((Vel_CartMesh.YRes - 1) * sizeof(double *))) == NULL)
			FatalError("Malloc failed for tempptr[%d]", i);
	
	/*** Swap addresses ***/
	/* wx */
	for(i = 0; i < Vel_CartMesh.XRes - 1; i++)
		for(j = 0; j < Vel_CartMesh.YRes - 1; j++)
			tempptr[i][j] = Vel_CartVorArray_wx[slot1][i][j]; /* Save address of slot1 */
	for(i = 0; i < Vel_CartMesh.XRes - 1; i++)
		for(j = 0; j < Vel_CartMesh.YRes - 1; j++)
			Vel_CartVorArray_wx[slot1][i][j] = Vel_CartVorArray_wx[slot2][i][j]; /* Change address of slot1 to address of slot2 */
	for(i = 0; i < Vel_CartMesh.XRes - 1; i++)
		for(j = 0; j < Vel_CartMesh.YRes - 1; j++)
			Vel_CartVorArray_wx[slot2][i][j] = tempptr[i][j]; /* Change address of slot2 to former address of slot1 */
	/* wy */
	for(i = 0; i < Vel_CartMesh.XRes - 1; i++)
		for(j = 0; j < Vel_CartMesh.YRes - 1; j++)
			tempptr[i][j] = Vel_CartVorArray_wy[slot1][i][j]; /* Save address of slot1 */
	for(i = 0; i < Vel_CartMesh.XRes - 1; i++)
		for(j = 0; j < Vel_CartMesh.YRes - 1; j++)
			Vel_CartVorArray_wy[slot1][i][j] = Vel_CartVorArray_wy[slot2][i][j]; /* Change address of slot1 to address of slot2 */
	for(i = 0; i < Vel_CartMesh.XRes - 1; i++)
		for(j = 0; j < Vel_CartMesh.YRes - 1; j++)
			Vel_CartVorArray_wy[slot2][i][j] = tempptr[i][j]; /* Change address of slot2 to former address of slot1 */
	/* wz */
	for(i = 0; i < Vel_CartMesh.XRes - 1; i++)
		for(j = 0; j < Vel_CartMesh.YRes - 1; j++)
			tempptr[i][j] = Vel_CartVorArray_wz[slot1][i][j]; /* Save address of slot1 */
	for(i = 0; i < Vel_CartMesh.XRes - 1; i++)
		for(j = 0; j < Vel_CartMesh.YRes - 1; j++)
			Vel_CartVorArray_wz[slot1][i][j] = Vel_CartVorArray_wz[slot2][i][j]; /* Change address of slot1 to address of slot2 */
	for(i = 0; i < Vel_CartMesh.XRes - 1; i++)
		for(j = 0; j < Vel_CartMesh.YRes - 1; j++)
			Vel_CartVorArray_wz[slot2][i][j] = tempptr[i][j]; /* Change address of slot2 to former address of slot1 */
	
	/*** Read in new time slice of data ***/
	/* Open file */
	sprintf(Data_BinFilePath, "%s%s_vorticity.%d.bin", Path_Data, Data_InFilePrefix, Data_SuffixTMin + ModVal * Data_SuffixTDelta);
	if((Data_BinFileID = fopen(Data_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", Data_BinFilePath);
	/* Read time stamp */
	if(fread(&timestamp, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read time stamp from file %s", Data_BinFilePath);
	printf("\nLoading vorticity data from %s (time stamp = %g)...", Data_BinFilePath, timestamp);
	fflush(stdout);
	/* Read vorticity data */
	for(k = 0; k < Vel_CartMesh.ZRes - 1; k++)
		for(j = 0; j < Vel_CartMesh.YRes - 1; j++)
			for(i = 0; i < Vel_CartMesh.XRes - 1; i++)
				if((fread(&Vel_CartVorArray_wx[slot2][i][j][k], sizeof(double), 1, Data_BinFileID) < 1) ||
				   (fread(&Vel_CartVorArray_wy[slot2][i][j][k], sizeof(double), 1, Data_BinFileID) < 1) ||
				   (fread(&Vel_CartVorArray_wz[slot2][i][j][k], sizeof(double), 1, Data_BinFileID) < 1)) /* Read data from file to address of slot2 */
					FatalError("Could not load vorticity data for index [%d][%d][%d] from %s", i, j, k, Data_BinFilePath);
	/* Close file */
	fclose(Data_BinFileID);
	
	/* Clean up local arrays */
	for(i = 0; i < Vel_CartMesh.XRes - 1; i++) { 
		free(tempptr[i]);
		tempptr[i] = NULL; 
	} 
	free(tempptr);
	tempptr = NULL;
	
	printf("OK!\n");
	fflush(stdout);
	
}


void LoadUnstructVorticityDataFrame(int frame) {
	
	int i, ModVal, slot1, slot2;
	double *tempptr_U, *tempptr_V, *tempptr_W, timestamp;
	char Data_BinFilePath[LONGSTRING];
	FILE *Data_BinFileID;
	
	if(Int_TimeDirection > 0) {
		slot1 = 0;
		slot2 = 1;
	}    
	else {
		slot1 = 1; 
		slot2 = 0;
	}
	
	ModVal = fmod(frame, Data_TRes - 1); 
	if(ModVal < 0) {
		ModVal += Data_TRes - 1;
	}
	
	/* Allocate memory for velocity data if needed */
	if(Vel_UnstructVorArray_wx == NULL || Vel_UnstructVorArray_wy == NULL || Vel_UnstructVorArray_wz == NULL) 
		AllocateVorticityFieldData();
	
	/* Open binary velocity file for reading */
	sprintf(Data_BinFilePath, "%s%s_vorticity.%d.bin", Path_Data, Data_InFilePrefix, Data_SuffixTMin + ModVal * Data_SuffixTDelta);
	if((Data_BinFileID = fopen(Data_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", Data_BinFilePath);
	
	/* Read time stamp */
	if(fread(&timestamp, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read time stamp from file %s", Data_BinFilePath);
	
	printf("\nLoading vorticity data from %s (time stamp = %g)...", Data_BinFilePath, timestamp);
	fflush(stdout);
    
	/* Virtual data shift */
	tempptr_U = Vel_UnstructVorArray_wx[slot1];
	tempptr_V = Vel_UnstructVorArray_wy[slot1];
	tempptr_W = Vel_UnstructVorArray_wz[slot1];
	Vel_UnstructVorArray_wx[slot1] = Vel_UnstructVorArray_wx[slot2];
	Vel_UnstructVorArray_wy[slot1] = Vel_UnstructVorArray_wy[slot2]; 
	Vel_UnstructVorArray_wz[slot1] = Vel_UnstructVorArray_wz[slot2];
	Vel_UnstructVorArray_wx[slot2] = tempptr_U;
	Vel_UnstructVorArray_wy[slot2] = tempptr_V;
	Vel_UnstructVorArray_wz[slot2] = tempptr_W;
	
	/* Read data */
	for(i = 0; i < Vel_MeshNumElements; i++) 
		if(fread(&Vel_UnstructVorArray_wx[slot2][i], sizeof(double), 1, Data_BinFileID) < 1 || 
		   fread(&Vel_UnstructVorArray_wy[slot2][i], sizeof(double), 1, Data_BinFileID) < 1 ||
		   fread(&Vel_UnstructVorArray_wz[slot2][i], sizeof(double), 1, Data_BinFileID) < 1) 
			FatalError("Could not read vorticity data from file %s", Data_BinFilePath);
        
	fclose(Data_BinFileID);
	printf("OK!\n");
	fflush(stdout);
	
}



void SetVorticity(double tc, LagrangianPoint *pt) {
	
	if(Data_MeshType == CARTESIAN)
		SetVorticity_Cartesian(tc, pt);
	else if(Data_MeshType == UNSTRUCTURED) 
		SetVorticity_Unstructured(tc, pt);
	else if(Data_MeshType == ANALYTIC)
		FatalError("Data_MeshType == ANALYTIC in SetVorticity() unsupported");
	else 
		FatalError("Unrecognized value for Data_MeshType in SetVorticity()");
	
}


void SetVorticity_Cartesian(double time, LagrangianPoint *pt) {
	
	double tloc;
	int    i, j, k;
	
	if(TestOutsideDomain(pt->X))
		FatalError("Attempting to interpolate vorticity at point that has left the domain");
	
	tloc = (time - Data_LoadedTMin) / (Data_LoadedTMax - Data_LoadedTMin);
	assert((tloc < (1 + TINY)) && (tloc > (0 - TINY)));
	
	i = (int)floor((pt->X[0] - Vel_CartMesh.XMin) / Vel_CartMesh.XDelta);
	assert((i >= 0) && (i <= (Vel_CartMesh.XRes - 1)));
	if(i == (Vel_CartMesh.XRes - 1)) 
		i = Vel_CartMesh.XRes - 2;
	j = (int)floor((pt->X[1] - Vel_CartMesh.YMin) / Vel_CartMesh.YDelta);
	assert((j >= 0) && (j <= (Vel_CartMesh.XRes - 1)));
	if(j == (Vel_CartMesh.YRes - 1)) 
		j = Vel_CartMesh.YRes - 2;
	if(Dimensions == 3) {
		k = (int)floor((pt->X[2] - Vel_CartMesh.ZMin) / Vel_CartMesh.ZDelta);
		assert((k >= 0) && (k <= (Vel_CartMesh.ZRes - 1)));
		if(k == (Vel_CartMesh.ZRes - 1)) 
			k = Vel_CartMesh.ZRes - 2;
	}
	else
		k = 0;
    
	pt->vorticity[0] = ((1 - tloc) * Vel_CartVorArray_wx[0][i][j][k] + tloc * Vel_CartVorArray_wx[1][i][j][k]);
    pt->vorticity[1] = ((1 - tloc) * Vel_CartVorArray_wy[0][i][j][k] + tloc * Vel_CartVorArray_wy[1][i][j][k]);
    pt->vorticity[2] = ((1 - tloc) * Vel_CartVorArray_wz[0][i][j][k] + tloc * Vel_CartVorArray_wz[1][i][j][k]);
	
}


void SetVorticity_Unstructured(double time, LagrangianPoint *pt) {
	
	double tloc;
	
	pt->ElementIndex = Get_Element_Local_Search(pt->X, pt->ElementIndex);
	
	if(pt->ElementIndex == -1) 
		FatalError("Attempting to interpolate vorticity at point with pt->ElementIndex = -1");
	
	/* Set tloc defining where in between loaded data to interpolate in time */
	tloc = (time - Data_LoadedTMin) / (Data_LoadedTMax - Data_LoadedTMin);
	if((tloc > (1 + TINY)) || (tloc < (0 - TINY))) 
		FatalError("tloc = %f and must be between 0 and 1", tloc);
	
	/* Interpolate vorticity at element in time */
	pt->vorticity[0] = (1 - tloc) * Vel_UnstructVorArray_wx[0][pt->ElementIndex] + tloc * Vel_UnstructVorArray_wx[1][pt->ElementIndex];
    pt->vorticity[1] = (1 - tloc) * Vel_UnstructVorArray_wy[0][pt->ElementIndex] + tloc * Vel_UnstructVorArray_wy[1][pt->ElementIndex];
    pt->vorticity[2] = (1 - tloc) * Vel_UnstructVorArray_wz[0][pt->ElementIndex] + tloc * Vel_UnstructVorArray_wz[1][pt->ElementIndex];
    
    /* printf("vorticity = (%f, %f, %f)\n", pt->vorticity[0], pt->vorticity[1], pt->vorticity[2]);*/
    
	
}
