/*
 *  residencetime.c
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "globals.h"
#include "io.h"
#include "macros.h"
#include "residencetime.h"
#include "structs.h"
#include "tracers.h"

void OutputRT(void) {
	
	FILE *OutFileID, *InFileID, *BinFileID;
	char InFilePath[LONGSTRING], OutFilePath[LONGSTRING], buf[LONGSTRING];
	int ii, ss, NumElements, **ConnectivityArray;
	double time, residence_time;
	LagrangianPoint *TraceIC_Array;
	
	printf("\nOutputting residence time data to file...");
	fflush(stdout);
	
	if(Trace_GenerateMesh) { 
		/* Open file to store Cartesian mesh parameters */
		sprintf(OutFilePath, "%s%s_Cartesian.bin", Path_Output, Trace_RTOutFilePrefix);
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
	else {
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
		sprintf(OutFilePath, "%s%s_coordinates.bin", Path_Output, Trace_RTOutFilePrefix);
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
					FatalError("Malloc failed for ConnectivityArray[ii]");
			
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
	}
	
	/* Output residence time data */  
	for(ss = 0; ss < Trace_NumLaunchTimes; ss++) {
		
		/* Load in data for slide ss into Trace_MeshPt */
		ReadInTraceLaunch(ss);  
		
		/* Open output file */
		sprintf(OutFilePath, "%s%s.%d.bin", Path_Output, Trace_RTOutFilePrefix, ss);
		if((OutFileID = fopen(OutFilePath, "wb")) == NULL) 
			FatalError("Could not open file %s", OutFilePath);
		
		/* Print time stamp */
		time = Output_TStart + ss*Trace_LaunchTimeSpacing;
		if(fwrite(&time, sizeof(double), 1, OutFileID) < 1)
			FatalError("Could not write time stamp to %s", OutFilePath);
		
		for(ii = 0; ii < Trace_NumTracers; ii++) {
			if(fabs(Trace_MeshPt[ii].LeftDomainTime + 1) < TINY) 
				residence_time = -1.0;
			else {
				if(!Trace_MeshPt[ii].LeftDomain)  
					Trace_MeshPt[ii].LeftDomainTime = Trace_Launches[ss].StopTime; 
				residence_time = Trace_MeshPt[ii].LeftDomainTime - (Output_TStart + ss*Trace_LaunchTimeSpacing);
			}
			if(fwrite(&residence_time, sizeof(double), 1, OutFileID) < 1)
				FatalError("Could not write residence time data to %s", OutFilePath);
		}
		
		/* Close output file */
		fclose(OutFileID); 
	}
	
	printf("OK!\n");
	fflush(stdout);
	
}