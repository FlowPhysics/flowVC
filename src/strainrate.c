/*
 *  strainrate.c
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "globals.h"
#include "io.h"
#include "macros.h"
#include "memory.h"
#include "mesh.h"
#include "strainrate.h"
#include "structs.h"
#include "tracers.h"
#include "velocity.h"

void AllocateStrainRateData(void) {
	/*** 
	 Allocate memory for strain rate field for 2 time frames (moving buffer of data)
	 ***/
	printf("Allocating memory for strain rate data...");
	if(Data_MeshType == CARTESIAN) {
		if (Dimensions == 2) 
			Vel_CartStrainRateArray = Calloc4D(Vel_CartStrainRateArray, 2, Vel_CartMesh.XRes - 1, Vel_CartMesh.YRes - 1, 1);
		else
			Vel_CartStrainRateArray = Calloc4D(Vel_CartStrainRateArray, 2, Vel_CartMesh.XRes - 1, Vel_CartMesh.YRes - 1, Vel_CartMesh.ZRes - 1);
	}
	else 
		Vel_UnstructStrainRateArray = Calloc2D(Vel_UnstructStrainRateArray, 2, Vel_MeshNumElements);
	printf("OK\n");
	fflush(stdout);
}

void FreeStrainRateData(void) {
	/*** 
	 Deallocate memory for strain rate field for 2 time frames (moving buffer of data)
	 ***/
	if (Data_MeshType == CARTESIAN) 
		Free4D(Vel_CartStrainRateArray, 2, Vel_CartMesh.XRes - 1, Vel_CartMesh.YRes - 1);
	else 
		Free2D(Vel_UnstructStrainRateArray, 2);
}


void LoadStrainRateDataFrame(int frame) {
    if(Data_MeshType == CARTESIAN)  
        LoadCartStrainRateDataFrame(Data_FirstFrame);
    else if(Data_MeshType == UNSTRUCTURED)
        LoadUnstructStrainRateDataFrame(Data_FirstFrame);
    else
        FatalError("Unknown value for Data_MeshType in LoadStrainRateDataFrame()");
}

void LoadCartStrainRateDataFrame(int frame) {
	
	int i, j, k, imax, jmax, kmax, ModVal = 0, slot1, slot2;
	double ***tempptr, timestamp;
	char Data_BinFilePath[LONGSTRING];
	FILE *Data_BinFileID;
	
	/* 
	 Following variables help ensure that data loaded in 
	 Vel_CartVelArray_*[0][*][*][*] corresponds to earlier 
	 time than data loaded into Vel_CartVelArray_*[1][*][*][*]
	 */
	if(Int_TimeDirection > 0) {
		slot1 = 0;
		slot2 = 1;
	}    
	else {
		slot1 = 1; 
		slot2 = 0;
	}
	
	if((frame > Data_TRes || frame < 0) && !Data_TPeriodic) 
		FatalError("Attempting to load data outside of available data time frame");
	else {
		ModVal = fmod(frame, Data_TRes - 1); 
		if(ModVal < 0) 
			ModVal += Data_TRes - 1;
	}
	
	if(Vel_CartVelArray_U == NULL || Vel_CartVelArray_V == NULL || Vel_CartVelArray_W == NULL)
		FatalError("Velocity mesh parameters do not appear to be loaded into memory");
	
	/* Allocate memory for strain-rate data if needed */
	if(Vel_CartStrainRateArray == NULL) 
		AllocateStrainRateData();
	
	imax = Vel_CartMesh.XRes - 1;
	jmax = Vel_CartMesh.YRes - 1;
	if ((Dimensions == 2) || (Vel_CartMesh.ZRes == 1)) 
		kmax = 1;
	else
		kmax = Vel_CartMesh.ZRes - 1;
	
	/* Allocate memory for tempptr, used as placeholder when swapping pointers to data as window of loaded strain rate data changes */
	if((tempptr = (double ***)malloc(imax * sizeof(double **))) == NULL)
		FatalError("Malloc failed for tempptr");
	for(i = 0; i < imax; i++) 
		if((tempptr[i] = (double **)malloc(jmax * sizeof(double *))) == NULL)
			FatalError("Malloc failed for tempptr[%d]", i);
	/* Swap memory */
	for(j = 0; j < jmax; j++)
		for(i = 0; i < imax; i++)
			tempptr[i][j] = Vel_CartStrainRateArray[slot1][i][j]; /* Save address of slot1 */
	for(j = 0; j < jmax; j++)
		for(i = 0; i < imax; i++)
			Vel_CartStrainRateArray[slot1][i][j] = Vel_CartStrainRateArray[slot2][i][j]; /* Change address of slot1 to address of slot2 */
	for(j = 0; j < jmax; j++)
		for(i = 0; i < imax; i++)
			Vel_CartStrainRateArray[slot2][i][j] = tempptr[i][j]; /* Change address of slot2 to former address of slot1 */
	
	/* Open strain rate data file for reading */
	sprintf(Data_BinFilePath, "%s%s_strain-rate.%d.bin", Path_Data, Data_InFilePrefix, Data_SuffixTMin + ModVal * Data_SuffixTDelta);
	if((Data_BinFileID = fopen(Data_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", Data_BinFilePath);
	/* Read time stamp */
	if(fread(&timestamp, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read time stamp from file %s", Data_BinFilePath);
	printf("Loading strain rate data from %s (time stamp = %g)...", Data_BinFilePath, timestamp);
	fflush(stdout);
	/* Update strain rate array */
	for(k = 0; k < kmax; k++)
		for(j = 0; j < jmax; j++)
			for(i = 0; i < imax; i++) 
				if(fread(&Vel_CartStrainRateArray[slot2][i][j][k], sizeof(double), 1, Data_BinFileID) < 1) /* Read data from file to address of slot2 */
					FatalError("Could not read strain rate index [%d][%d][%d] from %s", i, j, k, Data_BinFilePath);
	/* Close file */
	fclose(Data_BinFileID);
	
	/* Clean up local arrays */
	for(i = 0; i < imax; i++) { 
		free(tempptr[i]);
		tempptr[i] = NULL; 
	} 
	free(tempptr);
	tempptr = NULL;
	
	printf("OK!\n");
	fflush(stdout);
	
}

void LoadUnstructStrainRateDataFrame(int frame) {
	
	int ModVal, slot1, slot2;
	double *tempptr_S, timestamp;
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
	
	/* Open binary strain rate file for reading */
	sprintf(Data_BinFilePath, "%s%s_strain-rate.%d.bin", Path_Data, Data_InFilePrefix, Data_SuffixTMin + ModVal * Data_SuffixTDelta);
	if((Data_BinFileID = fopen(Data_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", Data_BinFilePath);
	
	/* Read time stamp */
	if(fread(&timestamp, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read time stamp from file %s", Data_BinFilePath);
	
	printf("Loading strain rate data from %s (time = %g)...", Data_BinFilePath, timestamp);
	fflush(stdout);
	
	/* Allocate memory for strain-rate data if needed */
	if(Vel_UnstructStrainRateArray == NULL) 
		AllocateStrainRateData();
	
	/* Virtual data shift */
	tempptr_S = Vel_UnstructStrainRateArray[slot1];
	Vel_UnstructStrainRateArray[slot1] = Vel_UnstructStrainRateArray[slot2];
	Vel_UnstructStrainRateArray[slot2] = tempptr_S;
	
	/* Read data */
	if((int)fread(Vel_UnstructStrainRateArray[slot2], sizeof(double), Vel_MeshNumElements, Data_BinFileID) < Vel_MeshNumElements)
		FatalError("Could not read UnstructStrainRateArray data from file %s", Data_BinFilePath);
	
	fclose(Data_BinFileID);
	printf("OK!\n");
	fflush(stdout);
	
}

double GetStrainRate(double tc, LagrangianPoint *pt) {
	
	double result = 0.0;
	
	if(Data_MeshType == CARTESIAN)
		result = GetStrainRate_Cartesian(tc, pt);
	else if(Data_MeshType == UNSTRUCTURED) 
		result = GetStrainRate_Unstructured(tc, pt);
	else if(Data_MeshType == ANALYTIC)
		FatalError("Data_MeshType == ANALYTIC in GetStrainRate() unsupported");
	else 
		FatalError("Unrecognized value for Data_MeshType in GetStrainRate()");
	
	return result;
}

double GetStrainRate_Cartesian(double time, LagrangianPoint *pt) {
	
	double tloc;
	int    i, j, k;
	
	if(TestOutsideDomain(pt->X))
		FatalError("Attempting to interpolate strain rate at point that has left the domain");
	
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
    
	return ((1 - tloc) * Vel_CartStrainRateArray[0][i][j][k] + tloc * Vel_CartStrainRateArray[1][i][j][k]);
	
}

double GetStrainRate_Unstructured(double time, LagrangianPoint *pt) {
	
	double tloc;
	
	pt->ElementIndex = Get_Element_Local_Search(pt->X, pt->ElementIndex);
	
	if(pt->ElementIndex == -1) 
		FatalError("Attempting to interpolate strain rate at point with Element_Index = -1");
	
	/* Set tloc defining where in between loaded data to interpolate in time */
	tloc = (time - Data_LoadedTMin) / (Data_LoadedTMax - Data_LoadedTMin);
	if((tloc > (1 + TINY)) || (tloc < (0 - TINY))) 
		FatalError("tloc = %f and must be between 0 and 1", tloc);
	
	/* Interpolate strain rate at element in time */
	return (1 - tloc) * Vel_UnstructStrainRateArray[0][pt->ElementIndex] + tloc * Vel_UnstructStrainRateArray[1][pt->ElementIndex];
	
}

void OutputAP(void) { 
	
	FILE *OutFileID, *InFileID, *BinFileID;
	char InFilePath[LONGSTRING], OutFilePath[LONGSTRING], buf[LONGSTRING];
	int ii, ss, NumElements, **ConnectivityArray;
	double time;
	LagrangianPoint *TraceIC_Array;
	
	if(Trace_GenerateMesh) { 
		/* Open file to store Cartesian mesh parameters */
		sprintf(OutFilePath, "%s%s_Cartesian.bin", Path_Output, Trace_OutFilePrefix);
		if((OutFileID = fopen(OutFilePath, "wb")) == NULL) 
			FatalError("Could not open file %s", OutFilePath);
		/* Write out Cartesian mesh parameters */
		if(fwrite(&Trace_CartMesh.XMin, sizeof(double), 1, OutFileID) < 1 ||
		   fwrite(&Trace_CartMesh.XMax, sizeof(double), 1, OutFileID) < 1 ||
		   fwrite(&Trace_CartMesh.XRes, sizeof(int),    1, OutFileID) < 1 ||
		   fwrite(&Trace_CartMesh.YMin, sizeof(double), 1, OutFileID) < 1 ||
		   fwrite(&Trace_CartMesh.YMax, sizeof(double), 1, OutFileID) < 1 ||
		   fwrite(&Trace_CartMesh.YRes, sizeof(int),    1, OutFileID) < 1 ||
		   fwrite(&Trace_CartMesh.ZMin, sizeof(double), 1, OutFileID) < 1 ||
		   fwrite(&Trace_CartMesh.ZMax, sizeof(double), 1, OutFileID) < 1 ||
		   fwrite(&Trace_CartMesh.ZRes, sizeof(int),    1, OutFileID) < 1)
			FatalError("Could not write Trace_CartMesh data to %s", OutFilePath);
		fclose(OutFileID);
	}
	else { /* Not Cartesian mesh */    
		/*** Create Coordinate File ***/
		/* Read in grid point locations */
		if((BinFileID = fopen(TraceIC_BinFilePath, "rb")) == NULL) 
			FatalError("Could not open file %s", TraceIC_BinFilePath);
		/* Read Number of Tracers to file */
		fread(&Trace_NumTracers, sizeof(int), 1, BinFileID);
		/* Allocate memory for TraceIC_Array */
		if((TraceIC_Array = (LagrangianPoint *)malloc(Trace_NumTracers * sizeof(LagrangianPoint))) == NULL) 
			FatalError("Malloc failed for TraceIC_Array in OutputRT()");
		/* Read to TraceIC_Array */
		fread(TraceIC_Array, sizeof(LagrangianPoint), Trace_NumTracers, BinFileID);
		/* Close IC file */
		fclose(BinFileID);
		/* Open output file for coordinate data */
		sprintf(OutFilePath, "%s%s_coordinates.bin", Path_Output, Trace_OutFilePrefix);
		if((OutFileID = fopen(OutFilePath, "wb")) == NULL) 
			FatalError("Could not open file %s", OutFilePath);
		/* Write number of nodes */
		if(fwrite(&Trace_NumTracers, sizeof(int), 1, OutFileID) < 1)
			FatalError("Could not write Trace_NumTracers to %s", OutFilePath);
		/* Write coordinates */
		for(ii = 0; ii < Trace_NumTracers; ii++) 
			if(fwrite(TraceIC_Array[ii].X, sizeof(double), 3, OutFileID) < 3)
				FatalError("Could not completely write TraceIC_Array[%d].X data to %s", ii, OutFilePath);
		/* Close coordinates file */
		fclose(OutFileID);
		/* Free TraceIC_Array */
		free(TraceIC_Array);
		TraceIC_Array = NULL; 
		/*** Create Connectivity File ***/
		if(Trace_InFileFormat == VTK_UNSTRUCTURED) { /* Create a file listing node connectivity */
			if(Dimensions == 2)
				FatalError("Trace_InFileFormat = %d is not currently supported for 2D analysis", VTK_UNSTRUCTURED);
			/* Open input file of interpolation points to read connectivity data */
			sprintf(InFilePath, "%s%s", Path_Data, Trace_InFile);
			if((InFileID = fopen(InFilePath, "r")) == NULL) 
				FatalError("Could not open %s", InFilePath);      
			/* Read until connectivity data found (i.e. look for line that starts with word CELLS */
			while(1) {
				if(fgets(buf, LONGSTRING, InFileID) == NULL)
					FatalError("Connectivity data not found in %s", InFilePath);
				if(!strncmp(buf, "CELLS", 5))
					break;
			}
			/* Read in number of elements */
			sscanf(buf, "%*s %d %*d\n", &NumElements); 
			/* Allocate memory for connectivity data */
			if((ConnectivityArray = (int **)malloc(NumElements * sizeof(int *))) == NULL)
				FatalError("Malloc failed for ConnectivityArray");
			for(ii = 0; ii < NumElements; ii++)   
				if((ConnectivityArray[ii] = (int *)malloc(4 * sizeof(int))) == NULL)
					FatalError("Malloc failed for ConnectivityArray[%d]", ii);
			/* Read in data */
			for(ii = 0; ii < NumElements; ii++)   
				if(fscanf(InFileID, "%*d %d %d %d %d\n", &ConnectivityArray[ii][1], &ConnectivityArray[ii][2], &ConnectivityArray[ii][3], &ConnectivityArray[ii][4]) == EOF)
					FatalError("Could not completely read connectivity data from %s", InFilePath); 
			/* Close input file */
			fclose(InFileID);
			/* Open binary output file for connectivity data */
			sprintf(OutFilePath, "%s%s_connectivity.bin", Path_Output, Trace_RTOutFilePrefix);
			if((OutFileID = fopen(OutFilePath, "wb")) == NULL) 
				FatalError("Could not open %s", OutFilePath);
			/* Write out number of elements */
			if(fwrite(&NumElements, sizeof(int), 1, OutFileID) < 1)
				FatalError("Could not write NumElements to %s", OutFilePath);
			/* Write out connectivity */
			for(ii = 0; ii < NumElements; ii++)   
				if(fwrite(ConnectivityArray[ii], sizeof(int), 4, OutFileID) < 4)
					FatalError("Could not write ConnectivityArray[%d] to %s", ii, OutFilePath);
			/* Close connectivity file */
			fclose(OutFileID);
			/* Free ConnectivityArray */
			for(ii = 0; ii < NumElements; ii++)   
				free(ConnectivityArray[ii]);
			free(ConnectivityArray);
		}
	} /* Not Cartesian mesh */
	
	/* Output activation potenial (integrated rate of strain norm) field */  
	for(ss = 0; ss < Trace_NumLaunchTimes; ss++) {
		/* Load in data for slide ss into Trace_MeshPt */
		ReadInTraceLaunch(ss);  
		/* Open output file */
		sprintf(OutFilePath, "%s%s_AP.%d.bin", Path_Output, Trace_OutFilePrefix, ss);
		if((OutFileID = fopen(OutFilePath, "wb")) == NULL) 
			FatalError("Could not open file %s", OutFilePath);
		/* Print time stamp */
		time = Output_TStart + ss*Trace_LaunchTimeSpacing;
		if(fwrite(&time, sizeof(double), 1, OutFileID) < 1)
			FatalError("Could not write time stamp to %s", OutFilePath);
		/* Write data */
		for(ii = 0; ii < Trace_NumTracers; ii++) 
			if(fwrite(&Trace_MeshPt[ii].Scalar, sizeof(double), 1, OutFileID) < 1)
				FatalError("Could not write activation potential data to %s", OutFilePath);
		/* Close output file */
		fclose(OutFileID); 
	}
	
	printf("OK!\n");
	fflush(stdout);
	
}

void OutputStrainOut(double t, int slide) {
	
	int i, NumElements;
	double S;
	char VelOutInFilePath[LONGSTRING], VelOutFilePath[LONGSTRING],  VelOutMeshFilePath[LONGSTRING], buf[LONGSTRING];
	FILE *VelOutInFileID, *VelOutFileID, *VelOutMeshFileID;
	int **ConnectivityArray;
	
	printf("Outputting strain field at time %f...", t);  
	fflush(stdout);
	
	/* Output mesh information if this is the first time outputing data */
	if(slide == 0) {
		if(VelOut_GenerateMesh) { 
			sprintf(VelOutMeshFilePath, "%s%s_Cartesian.bin", Path_Output, VelOut_FilePrefix);
			if((VelOutMeshFileID = fopen(VelOutMeshFilePath, "wb")) == NULL) 
				FatalError("Could not open file %s", VelOutMeshFilePath);
			/* Write out Cartesian mesh parameters */
			if(fwrite(&VelOut_CartMesh.XMin, sizeof(double), 1, VelOutMeshFileID) < 1 ||
			   fwrite(&VelOut_CartMesh.XMax, sizeof(double), 1, VelOutMeshFileID) < 1 ||
			   fwrite(&VelOut_CartMesh.XRes, sizeof(int),    1, VelOutMeshFileID) < 1 ||
			   fwrite(&VelOut_CartMesh.YMin, sizeof(double), 1, VelOutMeshFileID) < 1 ||
			   fwrite(&VelOut_CartMesh.YMax, sizeof(double), 1, VelOutMeshFileID) < 1 ||
			   fwrite(&VelOut_CartMesh.YRes, sizeof(int),    1, VelOutMeshFileID) < 1 ||
			   fwrite(&VelOut_CartMesh.ZMin, sizeof(double), 1, VelOutMeshFileID) < 1 ||
			   fwrite(&VelOut_CartMesh.ZMax, sizeof(double), 1, VelOutMeshFileID) < 1 ||
			   fwrite(&VelOut_CartMesh.ZRes, sizeof(int),    1, VelOutMeshFileID) < 1)
				FatalError("Could not write VelOut_CartMesh data to %s", VelOutMeshFilePath);
			
		}
		else {
			/* Open output file for coordinate data */
			sprintf(VelOutMeshFilePath, "%s%s_coordinates.bin", Path_Output, VelOut_FilePrefix);
			if((VelOutMeshFileID = fopen(VelOutMeshFilePath, "wb")) == NULL) 
				FatalError("Could not open file %s", VelOutMeshFilePath);
			/* Write number of nodes */
			if(fwrite(&VelOut_NumPts, sizeof(int), 1, VelOutMeshFileID) < 1)
				FatalError("Could not write VelOut_NumPts to %s", VelOutMeshFilePath);
			/* Write coordinates */
			for(i = 0; i < VelOut_NumPts; i++) 
				if(fwrite(VelOut_Mesh[i].X, sizeof(double), 3, VelOutMeshFileID) < 3)
					FatalError("Could not completely write VelOut_Mesh[%d].X data to %s", i, VelOutMeshFilePath);
			fclose(VelOutMeshFileID);
			
			if(VelOut_InFileFormat == VTK_UNSTRUCTURED) {	/* Create a file listing node connectivity */
				if(Dimensions == 2)
					FatalError("The specified VelOut_InFileFormat is not currently supported for 2D analysis");
				
				/* Open input file of interpolation points to read connectivity data */
				sprintf(VelOutInFilePath, "%s%s", Path_Data, VelOut_InFile);
				if((VelOutInFileID = fopen(VelOutInFilePath, "r")) == NULL) 
					FatalError("Could not open %s", VelOutInFilePath);
				/* Read until connectivity data found */
				while(1) {
					if(fgets(buf, LONGSTRING, VelOutInFileID) == NULL)
						FatalError("Connectivity data not found in %s", VelOutInFilePath);
					if(!strncmp(buf, "CELLS", 5))
						break;
				}
				sscanf(buf, "%*s %d %*d\n", &NumElements); 
				/* Allocate memory for connectivity data */
				if((ConnectivityArray = (int **)malloc(NumElements * sizeof(int *))) == NULL)
					FatalError("Malloc failed for ConnectivityArray");
				for(i = 0; i < NumElements; i++)   
					if((ConnectivityArray[i] = (int *)malloc(4 * sizeof(int))) == NULL)
						FatalError("Malloc failed for ConnectivityArray[%d]", i);
				/* Read in data */
				for(i = 0; i < NumElements; i++)   
					if(fscanf(VelOutInFileID, "%*d %d %d %d %d\n", &ConnectivityArray[i][1], &ConnectivityArray[i][2], &ConnectivityArray[i][3], &ConnectivityArray[i][4]) == EOF)
						FatalError("Could not completely read connectivity data from %s", VelOutInFilePath); 
				/* Close input file */
				fclose(VelOutInFileID);
				
				/* Open binary output file for connectivity data */
				sprintf(VelOutMeshFilePath, "%s%s_connectivity.bin", Path_Output, VelOut_FilePrefix);
				if((VelOutMeshFileID = fopen(VelOutMeshFilePath, "wb")) == NULL) 
					FatalError("Could not open %s", VelOutMeshFilePath);
				/* Write out number of elements */
				if(fwrite(&NumElements, sizeof(int), 1, VelOutMeshFileID) < 1)
					FatalError("Could not write NumElements to %s", VelOutMeshFilePath);
				/* Write out connectivity */
				for(i = 0; i < NumElements; i++)   
					if(fwrite(ConnectivityArray[i], sizeof(int), 4, VelOutMeshFileID) < 4)
						FatalError("Could not write ConnectivityArray[%d] to %s", i, VelOutMeshFilePath);
				/* Close connectivity file */
				fclose(VelOutMeshFileID);
				
				/* Free ConnectivityArray */
				for(i = 0; i < NumElements; i++)   
					free(ConnectivityArray[i]);
				free(ConnectivityArray);
				
			}
		}
	}
	
	/* Open output file */
	sprintf(VelOutFilePath, "%s%s_strain.%d.bin", Path_Output, VelOut_FilePrefix, slide);
	if((VelOutFileID = fopen(VelOutFilePath, "wb")) == NULL) 
		FatalError("Could not open file %s", VelOutFilePath);
	
	/* Write time stamp */
	if(fwrite(&t, sizeof(double), 1, VelOutFileID) < 1) 
		FatalError("Could not write time stamp to %s", VelOutFilePath);
	
	/* Write strain data */
	for(i = 0; i < VelOut_NumPts; i++) {
		/* Get value at point */
		if(VelOut_Mesh[i].LeftDomain) 
			S = 0.0;	
		else {
			if(Data_MeshType == CARTESIAN)
				S = GetStrainRate_Cartesian(t, &VelOut_Mesh[i]);
			else if(Data_MeshType == UNSTRUCTURED) 
				S = GetStrainRate_Unstructured(t, &VelOut_Mesh[i]);
		}
		/* Output velocity at point */
		if(fwrite(&S, sizeof(double), 1, VelOutFileID) < 1)
			FatalError("Could not write interpolated strain rate data to %s", VelOutFilePath);
	}    
	
	fclose(VelOutFileID);
	
	printf("OK!\n");  
	fflush(stdout);
	
}


