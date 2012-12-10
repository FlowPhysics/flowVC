/*
 *  exposuretime.c
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "exposuretime.h"
#include "globals.h"
#include "io.h"
#include "macros.h"
#include "mesh.h"
#include "mymath.h"


void LoadCETMeshData(void) {
	/*** Reads in unstructured mesh information into CET_MeshElementArray and CET_MeshNodeArray */
	
	int i, elements;
	char Mesh_BinFilePath[LONGSTRING], FileName[SHORTSTRING];
	FILE *Mesh_BinFileID;
	
	/******************************* NODES *********************************/
	/* Open file to read node coordinates */
	sprintf(FileName, "%s_coordinates.bin", Trace_CETMeshPrefix);
	sprintf(Mesh_BinFilePath, "%s%s", Path_Data, FileName);
	if((Mesh_BinFileID = fopen(Mesh_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", Mesh_BinFilePath);
	
	/* Read number of nodes */
	if(fread(&CET_MeshNumNodes, sizeof(int), 1, Mesh_BinFileID) < 1) 
		FatalError("Could not read number of mesh nodes from %s", Mesh_BinFilePath);
	
	printf("Loading coordinate data for %d nodes...", CET_MeshNumNodes);  
	fflush(stdout);
	
	/* Allocate memory for Vel_MeshNodeArray */
	if((CET_MeshNodeArray = (double **)malloc(CET_MeshNumNodes * sizeof(double *))) == NULL)
		FatalError("Malloc failed for CET_MeshNodeArray");
	for(i = 0; i < CET_MeshNumNodes; i++) 
		if((CET_MeshNodeArray[i] = (double *)malloc(3 * sizeof(double))) == NULL)
			FatalError("Malloc failed for CET_MeshNodeArray[%d]", i);
	
	/* Read CET_MeshNodeArray information from file */
	for(i = 0; i < CET_MeshNumNodes; i++) 
		if(fread(CET_MeshNodeArray[i], sizeof(double), 3, Mesh_BinFileID) < 3)
			FatalError("Could not read nodal coordinates completely from %s", Mesh_BinFilePath);
	
	fclose(Mesh_BinFileID);
	printf("OK!\n"); 
	fflush(stdout);
	
	/****************** CONNECTIVITY AND ADJACENCY ***********************/
	/* Open connectivity file to read binary */
	sprintf(FileName, "%s_connectivity.bin", Trace_CETMeshPrefix);
	sprintf(Mesh_BinFilePath, "%s%s", Path_Data, FileName);
	if((Mesh_BinFileID = fopen(Mesh_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", Mesh_BinFilePath);
	
	/* Read number of elements */
	if(fread(&CET_MeshNumElements, sizeof(int), 1, Mesh_BinFileID) < 1)
		FatalError("Could not read number of mesh elements from %s", Mesh_BinFilePath);
	
	printf("Loading connectivity and adjacency data for %d elements...", CET_MeshNumElements); 
	fflush(stdout);
	
	/* Allocate memory for Vel_MeshElementArray */
	if((CET_MeshElementArray = (Element *)malloc(CET_MeshNumElements * sizeof(Element))) == NULL)
		FatalError("Malloc failed for CET_MeshElementArray");
	
	/* Read connectivity information from file */
	for(i = 0; i < CET_MeshNumElements; i++)
		if(fread(CET_MeshElementArray[i].Nodes, sizeof(int), 4, Mesh_BinFileID) < 4)
			FatalError("Could not read connectivity for element %d from %s", i, Mesh_BinFilePath);
	
	/* Close connectivity file */
	fclose(Mesh_BinFileID);
	
	/* Open adjacency file to read binary */
	sprintf(FileName, "%s_adjacency.bin", Trace_CETMeshPrefix);
	sprintf(Mesh_BinFilePath, "%s%s", Path_Data, FileName);
	if((Mesh_BinFileID = fopen(Mesh_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", Mesh_BinFilePath);
	
	/* Read number of elements */
	if(fread(&elements, sizeof(int), 1, Mesh_BinFileID) < 1)
		FatalError("Could not read number of mesh elements from %s", Mesh_BinFilePath);
	if(elements != CET_MeshNumElements)
		FatalError("Incompatible number of elements listed in connectivity and adjacency files");
	
	/* Read adjacency information from file */
	for(i = 0; i < CET_MeshNumElements; i++)
		if(fread(CET_MeshElementArray[i].NeighborIndex, sizeof(int), 4, Mesh_BinFileID) < 4)
			FatalError("Could not read adjacency for element %d from %s", i, Mesh_BinFilePath);
	
	/* Close adjacency file */
	fclose(Mesh_BinFileID);
	
	printf("OK!\n");  
	fflush(stdout);
	
}

void ComputeExposureTime(const double *X0, const double *X1, const int X1ElementIndex, const double stepsize) {
	/*** 
	 Computes CET for a point moving from X0 to X1 in an amount of time stepsize and adds to Trace_CETArray
	 ***/ 
	
	int i, s, Index0, Index1, Indextmp, Indexlast;
	double a, Xtmp[3];
	
	/* Identify elements containing location of X0 and X1 */
	if(Trace_CETCompute && Trace_CETAuxillaryMesh) 
		Index0 = Get_Element_Local_Search_Aux(X0, X1ElementIndex);
	else
		Index0 = Get_Element_Local_Search(X0, X1ElementIndex);
	if(Index0 < 0) {
		if(Trace_CETCompute && Trace_CETAuxillaryMesh) 
			Index0 = Get_Element_Global_Search_Aux(X0);
		else
			Index0 = Get_Element_Global_Search(X0); 
		if(Index0 < 0) 
			FatalError("Point X0 passed to ComputeExposureTime appears to be outside of velocity domain");   
	}
	Index1 = X1ElementIndex;
	
	if(Index0 == Index1) {
		/* Same element contains both X0 and X1, add entire step size to CET sum */
		Trace_CETArray[Index0].CETsum = Trace_CETArray[Index0].CETsum + stepsize;
	}
	else { 
		/* Path between X0 and X1 intersects multiple elements, discretize straight-line path from X0 to X1 into Trace_CETSubsteps points */
		Indexlast = Index0;
		for(s = 1; s < Trace_CETSubsteps; s++) { /* Loop over points (except first one) */
			/* Compute location of point */
			a = ((double)s) / ((double)(Trace_CETSubsteps - 1));
			/* printf("a = %f\n", a);*/
			for(i = 0; i < 3; i++) 
				Xtmp[i] = (1 - a) * X0[i] + a * X1[i]; 
			/*printf("%.9f %.9f %.9f\n", Xtmp[0], Xtmp[1], Xtmp[2]);*/
			/* Determine element containing point's location */
			if(Trace_CETCompute && Trace_CETAuxillaryMesh)
				Indextmp = Get_Element_Local_Search_Aux(Xtmp, Indexlast); 
			else
				Indextmp = Get_Element_Local_Search(Xtmp, Indexlast); 
			/*printf("index = %d\n", Indextmp);*/
			if(Indextmp < 0) {
				if(Trace_CETCompute && Trace_CETAuxillaryMesh)
					Indextmp = Get_Element_Global_Search_Aux(Xtmp); 
				else
					Indextmp = Get_Element_Global_Search(Xtmp); 
				if(Indextmp < 0) {
					fprintf(stderr, "Warning: Particle path intersected region outside boundary, disregarding contribution to exposure times\n");
					continue;
				}
			}
			/* Save time portion to this element */
			Trace_CETArray[Indextmp].CETsum = Trace_CETArray[Indextmp].CETsum + (stepsize / (Trace_CETSubsteps - 1)); 
			/* Increment number of encounters element Indextmp if if previous loction on path was inside a different element */
			if(Indextmp != Indexlast) {
				Trace_CETArray[Indextmp].Encounters++;
				Indexlast = Indextmp;
			}
		}
	} 
}

void OutputCET(void) {
	
	FILE *Trace_CETOutFileID, *Trace_METOutFileID, *Trace_EncOutFileID;
	char Trace_CETOutFilePath[LONGSTRING], Trace_METOutFilePath[LONGSTRING], Trace_EncOutFilePath[LONGSTRING];
	int i;
	double t, vol, xc, *CET, *MET, X1[3], X2[3], X3[3], X4[3], A[3], B[3], C[3], D[3];
	
	printf("\nOutputting exposure time data to file...");
	fflush(stdout);
	
	/* Allocate memory for CET */
	if((CET = (double *)calloc(CET_MeshNumElements, sizeof(double))) == NULL) 
		FatalError("Calloc failed for CET");
	
	/* Allocate memory for MET */
	if((MET = (double *)calloc(CET_MeshNumElements, sizeof(double))) == NULL) 
		FatalError("Calloc failed for MET");
	
	for(i = 0; i < CET_MeshNumElements; i++) {
		if(Trace_CETArray[i].CETsum > TINY) {
			
			if(Dimensions == 3) {
				/* Set up coordinates of tet nodes */
				X1[0] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[0]][0];
				X1[1] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[0]][1];
				X1[2] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[0]][2];
				X2[0] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[1]][0];
				X2[1] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[1]][1];
				X2[2] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[1]][2];
				X3[0] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[2]][0];
				X3[1] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[2]][1];
				X3[2] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[2]][2];
				X4[0] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[3]][0];
				X4[1] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[3]][1];
				X4[2] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[3]][2];
				
				/* Solve the equation tet volume = |((X1 - X4) . ((X2 - X4) x (X3 - X4))| / 6 */
				A[0] = X1[0] - X4[0];
				A[1] = X1[1] - X4[1];
				A[2] = X1[2] - X4[2];
				B[0] = X2[0] - X4[0];
				B[1] = X2[1] - X4[1];
				B[2] = X2[2] - X4[2];
				C[0] = X3[0] - X4[0];
				C[1] = X3[1] - X4[1];
				C[2] = X3[2] - X4[2];
				cross(B, C, D);
				vol = fabs(vdot(A, D, 3)) / 6;
                xc = pow(vol,(1/3));
			}
			else {
				/* Set up coordinates of tet nodes */
				X1[0] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[0]][0];
				X1[1] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[0]][1];
				X2[0] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[1]][0];
				X2[1] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[1]][1];
				X3[0] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[2]][0];
				X3[1] = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[2]][1];
				
				/* Solve the equation triangle volume = |(X2 - X1) x (X3 - X1)| / 2 */
				A[0] = X2[0] - X1[0];
				A[1] = X2[1] - X1[1];
				B[0] = X3[0] - X1[0];
				B[1] = X3[1] - X1[1];
				vol = fabs(A[0] * B[1] - A[1] * B[0]) / 2;
                xc = sqrt(vol);
			}
			/* Normalize CET by element volume */
			CET[i] = Trace_CETArray[i].CETsum / xc / Trace_NumTracers; 
			MET[i] = Trace_CETArray[i].CETsum / xc / Trace_CETArray[i].Encounters;
		}
		else {
			CET[i] = 0;
			MET[i] = 0;
		}
	}
	
	/* Open output files for writing binary */
	sprintf(Trace_CETOutFilePath, "%s%s_CET.0.bin", Path_Output, Trace_RTOutFilePrefix);
	if((Trace_CETOutFileID = fopen(Trace_CETOutFilePath, "wb")) == NULL) {
		fprintf(stderr, "Error: Could not open %s\n", Trace_CETOutFilePath);
		exit(1);
	}
	sprintf(Trace_METOutFilePath, "%s%s_MET.0.bin", Path_Output, Trace_RTOutFilePrefix);
	if((Trace_METOutFileID = fopen(Trace_METOutFilePath, "wb")) == NULL) {
		fprintf(stderr, "Error: Could not open %s\n", Trace_METOutFilePath);
		exit(1);
	}
	sprintf(Trace_EncOutFilePath, "%s%s_ENC.0.bin", Path_Output, Trace_RTOutFilePrefix);
	if((Trace_EncOutFileID = fopen(Trace_EncOutFilePath, "wb")) == NULL) {
		fprintf(stderr, "Error: Could not open %s\n", Trace_EncOutFilePath);
		exit(1);
	}
	
	/* Output dummy time step */
	t = -1.0;
	if(fwrite(&t, sizeof(double), 1, Trace_CETOutFileID) < 1)
		FatalError("Could not output time step to %s", Trace_CETOutFilePath);
	if(fwrite(&t, sizeof(double), 1, Trace_METOutFileID) < 1)
		FatalError("Could not output time step to %s", Trace_CETOutFilePath);
	if(fwrite(&t, sizeof(double), 1, Trace_EncOutFileID) < 1)
		FatalError("Could not output time step to %s", Trace_EncOutFilePath);
	
	/* Output data */
	if(fwrite(CET, sizeof(double), CET_MeshNumElements, Trace_CETOutFileID) < (unsigned int)CET_MeshNumElements)
		fprintf(stderr, "Warning: Could not output all CET data to bin file\n");
	if(fwrite(MET, sizeof(double), CET_MeshNumElements, Trace_METOutFileID) < (unsigned int)CET_MeshNumElements)
		fprintf(stderr, "Warning: Could not output all MET data to bin file\n");
	for(i = 0; i < CET_MeshNumElements; i++) 
		if(fwrite(&Trace_CETArray[i].Encounters, sizeof(int), 1, Trace_EncOutFileID) < 1)
			fprintf(stderr, "Warning: Could not output all encounters data to bin file\n");
	
	
	/* Close file */
	if(fclose(Trace_CETOutFileID))
		fprintf(stderr, "Warning: Could not close %s\n", Trace_CETOutFilePath);
	if(fclose(Trace_METOutFileID))
		fprintf(stderr, "Warning: Could not close %s\n", Trace_METOutFilePath);
	if(fclose(Trace_EncOutFileID))
		fprintf(stderr, "Warning: Could not close %s\n", Trace_EncOutFilePath);
	
	free(CET);
	CET = NULL;
	free(MET);
	MET = NULL;
	
	printf("OK!");
	fflush(stdout);
	
}