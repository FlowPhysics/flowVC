/*
 *  velout.c
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2011 Flow Physics Group. All rights reserved.
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "globals.h"
#include "io.h"
#include "macros.h"
#include "mesh.h"
#include "velocity.h"
#include "velout.h"

void GenerateVelOutMesh(void) {
	
	int i, j, k, seed = -1, index = -1, count = 0, found = 0, foundGS = 0, guess = -1, jguess = -1, kguess = -1; 
	double Xseed[3];
	
	printf("\nGenerating velocity output mesh...\n"); 
	fflush(stdout);
	
	if(VelOut_CartMesh.XRes == 1)
		VelOut_CartMesh.XDelta = 0;
	else
		VelOut_CartMesh.XDelta = (VelOut_CartMesh.XMax - VelOut_CartMesh.XMin) / (VelOut_CartMesh.XRes - 1);
	if(VelOut_CartMesh.YRes == 1)
		VelOut_CartMesh.YDelta = 0;
	else
		VelOut_CartMesh.YDelta = (VelOut_CartMesh.YMax - VelOut_CartMesh.YMin) / (VelOut_CartMesh.YRes - 1);
	if(VelOut_CartMesh.ZRes == 1)
		VelOut_CartMesh.ZDelta = 0;
	else {
		if(Dimensions == 3)
			VelOut_CartMesh.ZDelta = (VelOut_CartMesh.ZMax - VelOut_CartMesh.ZMin) / (VelOut_CartMesh.ZRes - 1);
		else {
			fprintf(stderr, "Warning: VelOut_CartMesh.ZRes being set to 1 since velocity data is 2D\n");
			VelOut_CartMesh.ZRes = 1;
			VelOut_CartMesh.ZDelta = 0;
		}
	}
	
	VelOut_NumPts = VelOut_CartMesh.XRes * VelOut_CartMesh.YRes * VelOut_CartMesh.ZRes;
	
	/* Allocate memory for velocity output mesh */
	if((VelOut_Mesh = (LagrangianPoint *)malloc(VelOut_NumPts * sizeof(LagrangianPoint))) == NULL) 
		FatalError("Malloc failed for VelOut_Mesh");
	
	/* Try to "seed" local element search by doing global search for center element */
	if(Data_MeshType == UNSTRUCTURED) {
		Xseed[0] = VelOut_CartMesh.XMin + ((int)(VelOut_CartMesh.XRes/2)) * VelOut_CartMesh.XDelta;
		Xseed[1] = VelOut_CartMesh.YMin + ((int)(VelOut_CartMesh.YRes/2)) * VelOut_CartMesh.YDelta;
		Xseed[2] = VelOut_CartMesh.ZMin + ((int)(VelOut_CartMesh.ZRes/2)) * VelOut_CartMesh.ZDelta;
		seed = Get_Element_Global_Search(Xseed);
		if(seed >= 0) {
			printf("  Using center element seed.\n");
			fflush(stdout);
			found++;
			guess = seed;
		}
	}
	
	for(k = 0; k < VelOut_CartMesh.ZRes; k++) {
		kguess = -1;
		for(j = 0; j < VelOut_CartMesh.YRes; j++) {
			jguess = -1;
			for(i = 0; i < VelOut_CartMesh.XRes; i++) {
				/* Set coordinates of point */
				VelOut_Mesh[count].X[0] = VelOut_CartMesh.XMin + i * VelOut_CartMesh.XDelta;
				VelOut_Mesh[count].X[1] = VelOut_CartMesh.YMin + j * VelOut_CartMesh.YDelta;
				VelOut_Mesh[count].X[2] = VelOut_CartMesh.ZMin + k * VelOut_CartMesh.ZDelta;
				
				/* See if coordinates outside of data bounding box */
				VelOut_Mesh[count].LeftDomain = TestOutsideDomain(VelOut_Mesh[count].X);
				if(!VelOut_Mesh[count].LeftDomain) {
					if(Data_MeshType == UNSTRUCTURED) {
						/* Determine element the point is located in */
						if(seed < 0) {
							seed = Get_Element_Global_Search(VelOut_Mesh[count].X);
							if(seed < 0) {
								VelOut_Mesh[count].LeftDomain = 1;
								VelOut_Mesh[count].ElementIndex = -1;
							}
							else {
								printf("  Using first found element seed.\n");
								fflush(stdout);
								VelOut_Mesh[count].ElementIndex = seed;
								found++;
								guess = seed;
								jguess = seed;
								kguess = seed;
							}
						}
						else { 
							index = Get_Element_Local_Search(VelOut_Mesh[count].X, guess);
							if(index < 0) {
								/* An element not found, try local search from seed location */ 
								index = Get_Element_Local_Search(VelOut_Mesh[count].X, seed);
								if(index < 0) {
									/* If an element still not found, do global search if local search checking requested */
									if(LocalSearchChecking) {
										index = Get_Element_Global_Search(VelOut_Mesh[count].X);
										if(index < 0) {
											/* Point definitely outside domain */
											VelOut_Mesh[count].LeftDomain = 1;
											VelOut_Mesh[count].ElementIndex = -1;
										}
										else {
											/* Element found by global search */
											VelOut_Mesh[count].ElementIndex = index;		     
											found++;
											foundGS++;
											guess = index;
											if(jguess < 0)
												jguess = index;
											if(kguess < 0)
												kguess = index;
										}
									}
									else {
										VelOut_Mesh[count].LeftDomain = 1;
										VelOut_Mesh[count].ElementIndex = -1;
									}
								}
								else {
									/* Element found by local search from seed location */
									VelOut_Mesh[count].ElementIndex = index;
									found++;
									guess = index;
									if(jguess < 0)
										jguess = index;
									if(kguess < 0)
										kguess = index;
								}
							}
							else {
								/* Element found by normal local search */
								VelOut_Mesh[count].ElementIndex = index;
								guess = index;
								if(jguess < 0)
									jguess = index;
								if(kguess < 0)
									kguess = index;
							}
						}
					}
					else /* Data_MeshType == CARTESIAN */
						found++;
				}
				count++;
			}
			if(jguess >= 0)
				guess = jguess;
		}
		if(kguess >= 0)
			guess = kguess;
	}
	
	if(LocalSearchChecking) 
		printf("  %d of %d located in domain (%d caught by local search checking)\n", found, VelOut_CartMesh.XRes * VelOut_CartMesh.YRes * VelOut_CartMesh.ZRes, foundGS);
	else
		printf("  %d of %d located in domain\n", found, VelOut_CartMesh.XRes * VelOut_CartMesh.YRes * VelOut_CartMesh.ZRes);  
	fflush(stdout);
	
	/* Set up interpolation times */
	if((VelOut_OutputTime = (double *)malloc(Output_TRes * sizeof(double))) == NULL) 
		FatalError("Malloc failed for VelOut_OutputTime");
	
	if((VelOut_Complete = (int *)malloc(Output_TRes * sizeof(int))) == NULL) 
		FatalError("Malloc failed for VelOut_Complete");
	
	for(i = 0; i < Output_TRes; i++) {
		VelOut_OutputTime[i] = Output_TStart + Int_TimeDirection * i * Output_TDelta;
		VelOut_Complete[i] = 0;
	}
	
	printf("OK!\n");
	fflush(stdout);
	
}

void ReadInVelOutMesh(void) {
	
	int ii, index = -1, found = 0, foundGS = 0, guess = -1;
	FILE *VelOut_InFileID;
	char VelOut_InFilePath[LONGSTRING], buf[LONGSTRING];
	
	sprintf(VelOut_InFilePath, "%s%s", Path_Data, VelOut_InFile);
	if((VelOut_InFileID = fopen(VelOut_InFilePath, "r")) == NULL) 
		FatalError("Could not open file %s", VelOut_InFilePath);
	
	if(VelOut_InFileFormat == ASCII_LIST) { 
		/* First line should contain number of nodal data points to read in */
		fscanf(VelOut_InFileID, "%d\n", &VelOut_NumPts);
		printf("\nReading and initializing data for %d velocity interpolation nodes...", VelOut_NumPts); 
		fflush(stdout);
		
		/* Allocate memory */
		if((VelOut_Mesh = (LagrangianPoint *)malloc(VelOut_NumPts * sizeof(LagrangianPoint))) == NULL) 
			FatalError("Malloc failed for VelOut_Mesh");
		
		/* Read in nodal coordinates */
		for(ii = 0; ii < VelOut_NumPts; ii++) {
			if(Dimensions == 3)
				fscanf(VelOut_InFileID, "%lf %lf %lf\n", &VelOut_Mesh[ii].X[0], &VelOut_Mesh[ii].X[1], &VelOut_Mesh[ii].X[2]);
			else { /* Dimensions == 2 */
				fscanf(VelOut_InFileID, "%lf %lf\n", &VelOut_Mesh[ii].X[0], &VelOut_Mesh[ii].X[1]);
				VelOut_Mesh[ii].X[2] = 0.0;
			}
		}
		
		fclose(VelOut_InFileID);
	}
	else if(VelOut_InFileFormat == VTK_POLYDATA || VelOut_InFileFormat == VTK_UNSTRUCTURED) {
		/* Skip first 4 header lines */
		for(ii = 0; ii < 4; ii++)
			fgets(buf, LONGSTRING, VelOut_InFileID);
		
		/* Get number of nodal data points to read in */
		fscanf(VelOut_InFileID, "%*s %d %*s\n", &VelOut_NumPts); 
		printf("\nReading and initializing data for %d velocity interpolation nodes...", VelOut_NumPts); 
		fflush(stdout);
		
		/* Allocate memory */
		if((VelOut_Mesh = (LagrangianPoint *)malloc(VelOut_NumPts * sizeof(LagrangianPoint))) == NULL) 
			FatalError("Malloc failed for VelOut_Mesh");    
		
		/* Read in nodal coordinates */
		for(ii = 0; ii < VelOut_NumPts; ii++) 
			if(fscanf(VelOut_InFileID, "%lf %lf %lf", &VelOut_Mesh[ii].X[0], &VelOut_Mesh[ii].X[1], &VelOut_Mesh[ii].X[2]) < 3)
				FatalError("Failed to read coordinates for point %d", ii);
		fclose(VelOut_InFileID);
	}
	
	/* Initialize Data */
	for(ii = 0; ii < VelOut_NumPts; ii++) {
		/* Check if point outside bounding box for data */
		VelOut_Mesh[ii].LeftDomain = TestOutsideDomain(VelOut_Mesh[ii].X);
		/* If point in bounding box and velocity data is unstrutured, find element that point is in */
		if(Data_MeshType == UNSTRUCTURED && !VelOut_Mesh[ii].LeftDomain) { 
			if(!found) {
				/* First point located using global search, and result serves as starting point for further (local) searching */
				index = Get_Element_Global_Search(VelOut_Mesh[ii].X);
				if(index < 0) {
					VelOut_Mesh[ii].LeftDomain = 1;
					VelOut_Mesh[ii].ElementIndex = -1;
				}
				else {
					VelOut_Mesh[ii].ElementIndex = index;
					found++;
					guess = index;
				}
			}
			else { 
				index = Get_Element_Local_Search(VelOut_Mesh[ii].X, guess);
				if(index < 0) { 
					if(LocalSearchChecking) {
						index = Get_Element_Global_Search(VelOut_Mesh[ii].X);
						if(index < 0) {
							VelOut_Mesh[ii].LeftDomain = 1;
							VelOut_Mesh[ii].ElementIndex = -1;
						}
						else {
							VelOut_Mesh[ii].ElementIndex = index;
							found++;
							foundGS++;
							guess = index;		
						}
					}
					else {
						VelOut_Mesh[ii].LeftDomain = 1;
						VelOut_Mesh[ii].ElementIndex = -1;
					}
				}
				else {
					VelOut_Mesh[ii].ElementIndex = index;
					found++;
					guess = index;
				}
			}	
		}
	}
	
	if(LocalSearchChecking) 
		printf("  %d of %d located in domain (%d caught by local search checking)\n", found, VelOut_NumPts, foundGS);
	else
		printf("  %d of %d located in domain\n", found, VelOut_NumPts);    
	
	/* Set up interpolation times */
	if((VelOut_OutputTime = (double *)malloc(Output_TRes * sizeof(double))) == NULL)
		FatalError("Malloc failed for VelOut_OutputTime");
	if((VelOut_Complete = (int *)calloc(Output_TRes, sizeof(int))) == NULL) 
		FatalError("Calloc failed for VelOut_Complete");
	for(ii = 0; ii < Output_TRes; ii++) 
		VelOut_OutputTime[ii] = Output_TStart + Int_TimeDirection * ii * Output_TDelta;
	
}

void OutputVelOut(double t, int slide) {
	
	int i, j, NumElements;
	double U[3];
	char VelOutInFilePath[LONGSTRING], VelOutFilePath[LONGSTRING],  VelOutMeshFilePath[LONGSTRING], buf[LONGSTRING];
	FILE *VelOutInFileID, *VelOutFileID, *VelOutMeshFileID;
	int **ConnectivityArray;
	
	printf("Outputting velocity field at time %f...", t);  
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
	sprintf(VelOutFilePath, "%s%s.%d.bin", Path_Output, VelOut_FilePrefix, slide);
	if((VelOutFileID = fopen(VelOutFilePath, "wb")) == NULL) 
		FatalError("Could not open file %s", VelOutFilePath);
	
	/* Write time stamp */
	if(fwrite(&t, sizeof(double), 1, VelOutFileID) < 1) 
		FatalError("Could not write time stamp to %s", VelOutFilePath);
	
	/* Write velocity data */
	for(i = 0; i < VelOut_NumPts; i++) {
		/* Get value at point */
		if(VelOut_Mesh[i].LeftDomain) {
			for(j = 0; j < 3; j++) 
				U[j] = 0.0;	
		}
		else {
			if(Data_MeshType == CARTESIAN)
				GetVelocity_Cartesian(t, &VelOut_Mesh[i], U);
			else if(Data_MeshType == UNSTRUCTURED) 
				GetVelocity_Unstructured(t, &VelOut_Mesh[i], U);
		}
		/* Output velocity at point */
		if(fwrite(U, sizeof(double), 3, VelOutFileID) < 3)
			FatalError("Could not completely write interpolated velocity data to %s", VelOutFilePath);
	}    
	
	fclose(VelOutFileID);
	
	printf("OK!\n");  
	fflush(stdout);
	
}

void FreeVelOutData(void) {
	
	free(VelOut_OutputTime);
	VelOut_OutputTime = NULL;
	free(VelOut_Complete);
	VelOut_Complete = NULL;
	free(VelOut_Mesh);
	VelOut_Mesh = NULL;
	
}