/*
 *  mesh.c
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "exposuretime.h"
#include "globals.h"
#include "io.h"
#include "macros.h"
#include "memory.h"
#include "mesh.h"


void LoadMeshData(void) {
	
	printf("Loading mesh data...");
	if (Data_MeshType == CARTESIAN) 
		LoadCartMeshData();
	else if (Data_MeshType == UNSTRUCTURED)  {
		LoadUnstructMeshData();    
		if(Trace_Compute && Trace_CETCompute) {
			if(Trace_CETAuxillaryMesh)
				LoadCETMeshData(); /* Load CET auxillary mesh data */
			else {
				printf("Auxillary mesh not being used.\n");
				fflush(stdout);
				CET_MeshNumNodes = Vel_MeshNumNodes;
				CET_MeshNumElements = Vel_MeshNumElements;
				CET_MeshNodeArray = Vel_MeshNodeArray;
				CET_MeshElementArray = Vel_MeshElementArray;
			}
		}
	}
	printf("OK!\n");
	
}

void FreeMeshData(void) {
	if(Data_MeshType == UNSTRUCTURED) {
		Free2D(Vel_MeshNodeArray, Vel_MeshNumNodes);
		free(Vel_MeshElementArray);
		Vel_MeshElementArray = NULL;
		if(Int_NormalFlow) {
			Free2D(Vel_SurfaceMeshInwardNormals, Vel_SurfaceMeshNumNodes);
			free(Vel_SurfaceMeshNodeIDs);
			Vel_SurfaceMeshNodeIDs = NULL;
			if(Particle_Radius > TINY) {
				free(Vel_MeshElementFlagArray);
				Vel_MeshElementFlagArray = NULL;
			}
		}
		if(Trace_Compute && Trace_CETCompute && Trace_CETAuxillaryMesh) {
			Free2D(CET_MeshNodeArray, CET_MeshNumNodes);
			free(CET_MeshElementArray);
			CET_MeshElementArray = NULL;
		}
	}
}

void LoadCartMeshData(void) {
	
	char Data_BinFilePath[LONGSTRING];
	FILE *Data_BinFileID;
	
	/* Open binary mesh file for reading */
	sprintf(Data_BinFilePath, "%s%s_Cartesian.bin", Path_Data, Data_InFilePrefix);
	if((Data_BinFileID = fopen(Data_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", Data_BinFilePath);
	/* Read grid parameters */
	if(fread(&Vel_CartMesh.XMin, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.XMin from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.XMax, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.XMax from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.XRes, sizeof(int), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.XRes from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.YMin, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.YMin from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.YMax, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.YMax from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.YRes, sizeof(int), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.YRes from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.ZMin, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.ZMin from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.ZMax, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.ZMax from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.ZRes, sizeof(int), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.ZRes from file %s", Data_BinFilePath);
	/* Close file */
	fclose(Data_BinFileID);
	/* Derived parameters */
	if(Vel_CartMesh.XRes < 2 || Vel_CartMesh.XMin >= Vel_CartMesh.XMax)
		FatalError("Insufficient resolution of velocity data in x-direction");
	else
		Vel_CartMesh.XDelta = (Vel_CartMesh.XMax - Vel_CartMesh.XMin) / (Vel_CartMesh.XRes - 1);
	if(Vel_CartMesh.YRes < 2 || Vel_CartMesh.YMin >= Vel_CartMesh.YMax)
		FatalError("Insufficient resolution of velocity data in y-direction");
	else
		Vel_CartMesh.YDelta = (Vel_CartMesh.YMax - Vel_CartMesh.YMin) / (Vel_CartMesh.YRes - 1);
	if(Dimensions == 3) {
		if(Vel_CartMesh.ZRes < 2 || Vel_CartMesh.ZMin >= Vel_CartMesh.ZMax)
			FatalError("Insufficient resolution of velocity data in z-direction");
		else
			Vel_CartMesh.ZDelta = (Vel_CartMesh.ZMax - Vel_CartMesh.ZMin) / (Vel_CartMesh.ZRes - 1);
	}
	else { /* Dimensions == 2 */
		Vel_CartMesh.ZMin = 0;
		Vel_CartMesh.ZMax = 0;
		Vel_CartMesh.ZRes = 1;
		Vel_CartMesh.ZDelta = 0;
	}
}

void LoadUnstructMeshData(void) {
	/*
	 *	Reads in unstructured mesh information into Vel_MeshElementArray and Vel_MeshNodeArray 
	 */
	
	int i, elements;
	char Mesh_BinFilePath[LONGSTRING], FileName[SHORTSTRING];
	FILE *Mesh_BinFileID;
	
	/******************************* NODES *********************************/
	/* Open file to read node coordinates */
	sprintf(FileName, "%s_coordinates.bin", Data_InFilePrefix);
	sprintf(Mesh_BinFilePath, "%s%s", Path_Data, FileName);
	if((Mesh_BinFileID = fopen(Mesh_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", Mesh_BinFilePath);
	
	/* Read number of nodes */
	if(fread(&Vel_MeshNumNodes, sizeof(int), 1, Mesh_BinFileID) < 1) 
		FatalError("Could not read number of mesh nodes from %s", Mesh_BinFilePath);
	
	printf("Loading coordinate data for %d nodes...", Vel_MeshNumNodes);  
	fflush(stdout);
	
	/* Allocate memory for Vel_MeshNodeArray */
	if((Vel_MeshNodeArray = (double **)malloc(Vel_MeshNumNodes * sizeof(double *))) == NULL)
		FatalError("Malloc failed for Vel_MeshNodeArray");
	for(i = 0; i < Vel_MeshNumNodes; i++) 
		if((Vel_MeshNodeArray[i] = (double *)malloc(3 * sizeof(double))) == NULL)
			FatalError("Malloc failed for Vel_MeshNodeArray[%d]", i);
	
	/* Read Vel_MeshNodeArray information from file */
	for(i = 0; i < Vel_MeshNumNodes; i++) 
		if(fread(Vel_MeshNodeArray[i], sizeof(double), 3, Mesh_BinFileID) < 3)
			FatalError("Could not read nodal coordinates completely from %s", Mesh_BinFilePath);
	
	
	fclose(Mesh_BinFileID);
	printf("OK!\n"); 
	fflush(stdout);
	
	/****************** CONNECTIVITY AND ADJACENCY ***********************/
	/* Open connectivity file to read binary */
	sprintf(FileName, "%s_connectivity.bin", Data_InFilePrefix);
	sprintf(Mesh_BinFilePath, "%s%s", Path_Data, FileName);
	if((Mesh_BinFileID = fopen(Mesh_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", Mesh_BinFilePath);
	
	/* Read number of elements */
	if(fread(&Vel_MeshNumElements, sizeof(int), 1, Mesh_BinFileID) < 1)
		FatalError("Could not read number of mesh elements from %s", Mesh_BinFilePath);
	
	printf("Loading connectivity and adjacency data for %d elements...", Vel_MeshNumElements); 
	fflush(stdout);
	
	/* Allocate memory for Vel_MeshElementArray */
	if((Vel_MeshElementArray = (Element *)malloc(Vel_MeshNumElements * sizeof(Element))) == NULL)
		FatalError("Malloc failed for Vel_MeshElementArray");
	
	/* Read connectivity information from file */
	for(i = 0; i < Vel_MeshNumElements; i++)
		if(fread(Vel_MeshElementArray[i].Nodes, sizeof(int), 4, Mesh_BinFileID) < 4)
			FatalError("Could not read connectivity for element %d from %s", i, Mesh_BinFilePath);
	
	/* Close connectivity file */
	fclose(Mesh_BinFileID);
	
	/* Open adjacency file to read binary */
	sprintf(FileName, "%s_adjacency.bin", Data_InFilePrefix);
	sprintf(Mesh_BinFilePath, "%s%s", Path_Data, FileName);
	if((Mesh_BinFileID = fopen(Mesh_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", Mesh_BinFilePath);
	
	/* Read number of elements */
	if(fread(&elements, sizeof(int), 1, Mesh_BinFileID) < 1)
		FatalError("Could not read number of mesh elements from %s", Mesh_BinFilePath);
	if(elements != Vel_MeshNumElements)
		FatalError("Incompatible number of elements listed in connectivity and adjacency files");
	
	/* Read adjacency information from file */
	for(i = 0; i < Vel_MeshNumElements; i++)
		if(fread(Vel_MeshElementArray[i].NeighborIndex, sizeof(int), 4, Mesh_BinFileID) < 4)
			FatalError("Could not read adjacency for element %d from %s", i, Mesh_BinFilePath);
	
	/* Close adjacency file */
	fclose(Mesh_BinFileID);
	
	printf("OK!\n");  
	fflush(stdout);
	
}

void LoadBoundaryElementFlags(void) {
	
	int elements;
	char BinFilePath[LONGSTRING];
	char BinFileName[SHORTSTRING];
	FILE *BinFileID;
	
	/* Open binary file for reading */
	sprintf(BinFileName, "%s_bflags.bin", Data_InFilePrefix);
	sprintf(BinFilePath, "%s%s", Path_Data, BinFileName);
	if((BinFileID = fopen(BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", BinFilePath);
	
	printf("\nReading boundary element flags from %s...", BinFileName);
	fflush(stdout);
	
	/* Read in the number of elements */
	if(fread(&elements, sizeof(int), 1, BinFileID) < 1) 
		FatalError("Could not read value for Vel_SurfaceMeshNumNodes from %s", BinFileName);
	if(elements != Vel_MeshNumElements)
		FatalError("Incompatable number of elements in %s", BinFileName);
	
	/* Allocate memory */
	if((Vel_MeshElementFlagArray = (int *)malloc(Vel_MeshNumElements * sizeof(int))) == NULL)
		FatalError("Malloc failed for Vel_MeshElementFlagArray");
	
	/* Read flags */
	if((int)fread(Vel_MeshElementFlagArray, sizeof(int), Vel_MeshNumElements, BinFileID) < Vel_MeshNumElements)
		FatalError("Could not completely read Vel_MeshNumElements from %s", BinFileName);
	
		
	fclose(BinFileID);
	
	printf("OK!\n");
	fflush(stdout);
	
}

	

void LoadMeshSurfaceNormals(void) {
	
	char BinFilePath[LONGSTRING];
	char BinFileName[SHORTSTRING];
	FILE *BinFileID;
	int i;
	
	/* Open binary file for reading */
	sprintf(BinFileName, "%s_normals.bin", Data_InFilePrefix);
	sprintf(BinFilePath, "%s%s", Path_Data, BinFileName);
	if((BinFileID = fopen(BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", BinFilePath);
	
	/* Read in the number of surface nodes */
	if(fread(&Vel_SurfaceMeshNumNodes, sizeof(int), 1, BinFileID) < 1) 
		FatalError("Could not read value for Vel_SurfaceMeshNumNodes from %s", BinFileName);
	
	printf("\nReading normal vectors for %d surface nodes from %s...", Vel_SurfaceMeshNumNodes, BinFileName);
	fflush(stdout);
	
	/* Allocate memory to store surface node ID numbers */
	if((Vel_SurfaceMeshNodeIDs = (int *)malloc(Vel_SurfaceMeshNumNodes * sizeof(int))) == NULL)
		FatalError("Malloc failed for Vel_SurfaceMeshNodeIDs");
	
	/* Read node ID's */
	if((int)fread(Vel_SurfaceMeshNodeIDs, sizeof(int), Vel_SurfaceMeshNumNodes, BinFileID) < Vel_SurfaceMeshNumNodes)
		FatalError("Could not completely read Vel_SurfaceMeshNodeIDs array from %s", BinFileName);
	
	/* Allocate memory to store normal vectors at surface nodes */
	if((Vel_SurfaceMeshInwardNormals = (double **)malloc(Vel_SurfaceMeshNumNodes * sizeof(double *))) == NULL) 
		FatalError("Malloc failed for Vel_SurfaceMeshInwardNormals");
	for(i = 0; i < Vel_SurfaceMeshNumNodes; i++) 
		if((Vel_SurfaceMeshInwardNormals[i] = (double *)malloc(3 * sizeof(double))) == NULL) 
			FatalError("Malloc failed for Vel_SurfaceMeshInwardNormals[%d]", i);
	
	/* Read in data */
	for(i = 0; i < Vel_SurfaceMeshNumNodes; i++) 
		if(fread(&Vel_SurfaceMeshInwardNormals[i][0], sizeof(double), 1, BinFileID) < 1 || 
		   fread(&Vel_SurfaceMeshInwardNormals[i][1], sizeof(double), 1, BinFileID) < 1 || 
		   fread(&Vel_SurfaceMeshInwardNormals[i][2], sizeof(double), 1, BinFileID) < 1) 
			FatalError("Could not completely read inward normal vector data at surface node %d from %s", i, BinFileName);
	
	fclose(BinFileID);
	printf("OK!\n");
	fflush(stdout);
	
} 


int Get_Element_Global_Search(const double *X) {
	/*** 
	 Returns index of Vel_MeshElementArray for element that contains point X
	 Returns -1 if point does not belong to any element 
	 ***/
	
	int i;
	double x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3; 
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;
	double r, s, t, d;
	double V;
	
	x = X[0];
	y = X[1];
	z = X[2];
	
	i = 0;
	while(i < Vel_MeshNumElements) { /* Search over all elements */
		if(Dimensions == 3) {
			/* Physical coordinates of nodes of the test element */
			x0 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[0]][0];
			y0 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[0]][1];
			z0 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[0]][2];
			x1 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[1]][0];
			y1 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[1]][1];
			z1 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[1]][2];
			x2 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[2]][0];
			y2 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[2]][1];
			z2 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[2]][2];
			x3 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[3]][0];
			y3 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[3]][1];
			z3 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[3]][2];
			
			/* Entries for mapping of physical to natural coordinates of test element */
			a11 = (z3 - z0) * (y2 - y3) - (z2 - z3) * (y3 - y0);
			a21 = (z3 - z0) * (y0 - y1) - (z0 - z1) * (y3 - y0);
			a31 = (z1 - z2) * (y0 - y1) - (z0 - z1) * (y1 - y2);
			a12 = (x3 - x0) * (z2 - z3) - (x2 - x3) * (z3 - z0);
			a22 = (x3 - x0) * (z0 - z1) - (x0 - x1) * (z3 - z0);
			a32 = (x1 - x2) * (z0 - z1) - (x0 - x1) * (z1 - z2);
			a13 = (y3 - y0) * (x2 - x3) - (y2 - y3) * (x3 - x0);
			a23 = (y3 - y0) * (x0 - x1) - (y0 - y1) * (x3 - x0);
			a33 = (y1 - y2) * (x0 - x1) - (y0 - y1) * (x1 - x2);			 
			
			/* Determinant of mapping from natural to physical coordinates of test element */
			V = (x1 - x0) * ((y2 - y0) * (z3 - z0) - (z2 - z0) * (y3 - y0)) + 
			(x2 - x0) * ((y0 - y1) * (z3 - z0) - (z0 - z1) * (y3 - y0)) +
			(x3 - x0) * ((y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0));
			
			/* Natural coordinates of point to be interpolated */
			r = (a11 * (x - x0) + a12 * (y - y0) + a13 * (z - z0)) / V;
			s = (a21 * (x - x0) + a22 * (y - y0) + a23 * (z - z0)) / V;
			t = (a31 * (x - x0) + a32 * (y - y0) + a33 * (z - z0)) / V;
			
			d = fmin(r, fmin(s, fmin(t, 1 - r - s - t)));
			
			if(d > -TINY) /* Point inside test element */
				return(i);
			else /* Check next element */
				i++;
		} 
		else { /* Dimensions == 2 */
			/* Physical coordinates of nodes of the test element */
			x0 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[0]][0];
			y0 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[0]][1];
			x1 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[1]][0];
			y1 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[1]][1];
			x2 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[2]][0];
			y2 = Vel_MeshNodeArray[Vel_MeshElementArray[i].Nodes[2]][1];
			
			/* Entries for mapping of physical to natural coordinates of test element */
			a11 = y2 - y0;
			a12 = x0 - x2;
			a21 = y0 - y1;
			a22 = x1 - x0;
			
			/* Determinant of mapping from natural to physical coordinates of test element */
			V = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
			
			/* Natural coordinates of point to be interpolated */
			r = (a11 * (x - x0) + a12 * (y - y0)) / V;
			s = (a21 * (x - x0) + a22 * (y - y0)) / V;
			
			d = fmin(r, fmin(s, 1 - r - s));
			
			if(d > -TINY) /* Point inside test element */
				return(i);
			else /* Check next element */
				i++;
		}
	}
	
	return -1;
	
}

int Get_Element_Global_Search_Aux(const double *X) {
	/*** 
	 Returns index of CET_MeshElementArray for element that contains point X
	 Returns -1 if point does not belong to any element 
	 ***/
	
	int i;
	double x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3; 
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;
	double r, s, t, d;
	double V;
	
	x = X[0];
	y = X[1];
	z = X[2];
	
	i = 0;
	while(i < CET_MeshNumElements) { /* Search over all elements */
		if(Dimensions == 3) {
			/* Physical coordinates of nodes of the test element */
			x0 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[0]][0];
			y0 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[0]][1];
			z0 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[0]][2];
			x1 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[1]][0];
			y1 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[1]][1];
			z1 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[1]][2];
			x2 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[2]][0];
			y2 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[2]][1];
			z2 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[2]][2];
			x3 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[3]][0];
			y3 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[3]][1];
			z3 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[3]][2];
			
			/* Entries for mapping of physical to natural coordinates of test element */
			a11 = (z3 - z0) * (y2 - y3) - (z2 - z3) * (y3 - y0);
			a21 = (z3 - z0) * (y0 - y1) - (z0 - z1) * (y3 - y0);
			a31 = (z1 - z2) * (y0 - y1) - (z0 - z1) * (y1 - y2);
			a12 = (x3 - x0) * (z2 - z3) - (x2 - x3) * (z3 - z0);
			a22 = (x3 - x0) * (z0 - z1) - (x0 - x1) * (z3 - z0);
			a32 = (x1 - x2) * (z0 - z1) - (x0 - x1) * (z1 - z2);
			a13 = (y3 - y0) * (x2 - x3) - (y2 - y3) * (x3 - x0);
			a23 = (y3 - y0) * (x0 - x1) - (y0 - y1) * (x3 - x0);
			a33 = (y1 - y2) * (x0 - x1) - (y0 - y1) * (x1 - x2);			 
			
			/* Determinant of mapping from natural to physical coordinates of test element */
			V = (x1 - x0) * ((y2 - y0) * (z3 - z0) - (z2 - z0) * (y3 - y0)) + 
			(x2 - x0) * ((y0 - y1) * (z3 - z0) - (z0 - z1) * (y3 - y0)) +
			(x3 - x0) * ((y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0));
			
			/* Natural coordinates of point to be interpolated */
			r = (a11 * (x - x0) + a12 * (y - y0) + a13 * (z - z0)) / V;
			s = (a21 * (x - x0) + a22 * (y - y0) + a23 * (z - z0)) / V;
			t = (a31 * (x - x0) + a32 * (y - y0) + a33 * (z - z0)) / V;
			
			d = fmin(r, fmin(s, fmin(t, 1 - r - s - t)));
			
			if(d > -TINY) /* Point inside test element */
				return(i);
			else /* Check next element */
				i++;
		} 
		else { /* Dimensions == 2 */
			/* Physical coordinates of nodes of the test element */
			x0 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[0]][0];
			y0 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[0]][1];
			x1 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[1]][0];
			y1 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[1]][1];
			x2 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[2]][0];
			y2 = CET_MeshNodeArray[CET_MeshElementArray[i].Nodes[2]][1];
			
			/* Entries for mapping of physical to natural coordinates of test element */
			a11 = y2 - y0;
			a12 = x0 - x2;
			a21 = y0 - y1;
			a22 = x1 - x0;
			
			/* Determinant of mapping from natural to physical coordinates of test element */
			V = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
			
			/* Natural coordinates of point to be interpolated */
			r = (a11 * (x - x0) + a12 * (y - y0)) / V;
			s = (a21 * (x - x0) + a22 * (y - y0)) / V;
			
			d = fmin(r, fmin(s, 1 - r - s));
			
			if(d > -TINY) /* Point inside test element */
				return(i);
			else /* Check next element */
				i++;
		}
	}
	
	return -1;
	
}

int Get_Element_Local_Search(const double *X, int guess) {
	/*** 
	 Returns index to Vel_MeshElementArray for element that contains point X
	 Returns -1 if element not located
	 ***/
	double x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;
	double r, s, t, d;
	double V;
	
	x = X[0];
	y = X[1];
	z = X[2];
	
	while(1) {
		if(Dimensions == 3) {
			/* Physical coordinates of nodes of the test element */
			x0 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[0]][0];
			y0 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[0]][1];
			z0 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[0]][2];
			x1 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[1]][0];
			y1 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[1]][1];
			z1 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[1]][2];
			x2 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[2]][0];
			y2 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[2]][1];
			z2 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[2]][2];
			x3 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[3]][0];
			y3 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[3]][1];
			z3 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[3]][2];
			
			/* Entries for mapping of physical to natural coordinates of test element */
			a11 = (z3 - z0) * (y2 - y3) - (z2 - z3) * (y3 - y0);
			a21 = (z3 - z0) * (y0 - y1) - (z0 - z1) * (y3 - y0);
			a31 = (z1 - z2) * (y0 - y1) - (z0 - z1) * (y1 - y2);
			a12 = (x3 - x0) * (z2 - z3) - (x2 - x3) * (z3 - z0);
			a22 = (x3 - x0) * (z0 - z1) - (x0 - x1) * (z3 - z0);
			a32 = (x1 - x2) * (z0 - z1) - (x0 - x1) * (z1 - z2);
			a13 = (y3 - y0) * (x2 - x3) - (y2 - y3) * (x3 - x0);
			a23 = (y3 - y0) * (x0 - x1) - (y0 - y1) * (x3 - x0);
			a33 = (y1 - y2) * (x0 - x1) - (y0 - y1) * (x1 - x2);			 
			
			/* Determinant of mapping from natural to physical coordinates of test element */
			V = (x1 - x0) * ((y2 - y0) * (z3 - z0) - (z2 - z0) * (y3 - y0)) + 
			(x2 - x0) * ((y0 - y1) * (z3 - z0) - (z0 - z1) * (y3 - y0)) +
			(x3 - x0) * ((y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0));
            
			/* Natural coordinates of point to be interpolated */
			r = (a11 * (x - x0) + a12 * (y - y0) + a13 * (z - z0)) / V;
			s = (a21 * (x - x0) + a22 * (y - y0) + a23 * (z - z0)) / V;
			t = (a31 * (x - x0) + a32 * (y - y0) + a33 * (z - z0)) / V;
			
			d = fmin(r, fmin(s, fmin(t, 1 - r - s - t)));
            
			if(d > -TINY) /* Point inside test element */
				return(guess);
			else { /* Reset test element to neighbor */	
				if(fabs(r - d) < TINY) 
					guess = Vel_MeshElementArray[guess].NeighborIndex[0];
				else if(fabs(s - d) < TINY) 
					guess = Vel_MeshElementArray[guess].NeighborIndex[1];
				else if(fabs(t - d) < TINY) 
					guess = Vel_MeshElementArray[guess].NeighborIndex[2];
				else if(fabs(d - (1 - r - s - t)) < TINY) 
					guess = Vel_MeshElementArray[guess].NeighborIndex[3];
				else {
					FatalError("Indeterminate neighbor in function Get_Element_Local_Search");
				}
				
				if(guess == -1) /* Neighbor not present, could not find element from local search (likely point left domain) */ 
					return(-1);
			}
		}
		else { /* Dimensions == 2 */
			/* Physical coordinates of nodes of the test element */
			x0 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[0]][0];
			y0 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[0]][1];
			x1 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[1]][0];
			y1 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[1]][1];
			x2 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[2]][0];
			y2 = Vel_MeshNodeArray[Vel_MeshElementArray[guess].Nodes[2]][1];
			
			/* printf("x = %f; y =  %f; in element %d x0 = %f; y0 = %f;  x1 = %f;  y1 = %f;  x2 = %f y2 = %f\n", x, y, guess, x0, y0, x1, y1, x2, y2);*/
			
			/* Entries for mapping of physical to natural coordinates of test element */
			a11 = y2 - y0;
			a12 = x0 - x2;
			a21 = y0 - y1;
			a22 = x1 - x0;
			
			/* Determinant of mapping from natural to physical coordinates of test element */
			V = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
			
			/* Natural coordinates of point to be interpolated */
			r = (a11 * (x - x0) + a12 * (y - y0)) / V;
			s = (a21 * (x - x0) + a22 * (y - y0)) / V;
			
			d = fmin(r, fmin(s, 1 - r - s));
			
			if(d > -TINY) /* Point inside test element */
				return(guess);
			else { /* Reset test element to neighbor */	
				if(fabs(r - d) < TINY) 
					guess = Vel_MeshElementArray[guess].NeighborIndex[0];
				else if(fabs(s - d) < TINY) 
					guess = Vel_MeshElementArray[guess].NeighborIndex[1];
				else if(fabs(d - (1 - r - s)) < TINY) 
					guess = Vel_MeshElementArray[guess].NeighborIndex[2];
				else 
					FatalError("Indeterminate neighbor in function Get_Element_Local_Search");
				
				if(guess == -1) /* Neighbor not present, could not find element from local search (point may have exited domain) */ 
					return(-1);
			}
		}
	}
	
}

int Get_Element_Local_Search_Aux(const double *X, int guess) {
	/*** 
	 Returns index to CET_MeshElementArray for element that contains point X
	 Returns -1 if element not located
	 ***/
	double x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;
	double r, s, t, d;
	double V;
	
	x = X[0];
	y = X[1];
	z = X[2];
	
	while(1) {
		if(Dimensions == 3) {
			/* Physical coordinates of nodes of the test element */
			x0 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[0]][0];
			y0 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[0]][1];
			z0 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[0]][2];
			x1 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[1]][0];
			y1 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[1]][1];
			z1 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[1]][2];
			x2 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[2]][0];
			y2 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[2]][1];
			z2 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[2]][2];
			x3 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[3]][0];
			y3 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[3]][1];
			z3 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[3]][2];
			
			/* Entries for mapping of physical to natural coordinates of test element */
			a11 = (z3 - z0) * (y2 - y3) - (z2 - z3) * (y3 - y0);
			a21 = (z3 - z0) * (y0 - y1) - (z0 - z1) * (y3 - y0);
			a31 = (z1 - z2) * (y0 - y1) - (z0 - z1) * (y1 - y2);
			a12 = (x3 - x0) * (z2 - z3) - (x2 - x3) * (z3 - z0);
			a22 = (x3 - x0) * (z0 - z1) - (x0 - x1) * (z3 - z0);
			a32 = (x1 - x2) * (z0 - z1) - (x0 - x1) * (z1 - z2);
			a13 = (y3 - y0) * (x2 - x3) - (y2 - y3) * (x3 - x0);
			a23 = (y3 - y0) * (x0 - x1) - (y0 - y1) * (x3 - x0);
			a33 = (y1 - y2) * (x0 - x1) - (y0 - y1) * (x1 - x2);			 
			
			/* Determinant of mapping from natural to physical coordinates of test element */
			V = (x1 - x0) * ((y2 - y0) * (z3 - z0) - (z2 - z0) * (y3 - y0)) + 
			(x2 - x0) * ((y0 - y1) * (z3 - z0) - (z0 - z1) * (y3 - y0)) +
			(x3 - x0) * ((y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0));
			
			/* Natural coordinates of point to be interpolated */
			r = (a11 * (x - x0) + a12 * (y - y0) + a13 * (z - z0)) / V;
			s = (a21 * (x - x0) + a22 * (y - y0) + a23 * (z - z0)) / V;
			t = (a31 * (x - x0) + a32 * (y - y0) + a33 * (z - z0)) / V;
			
			d = fmin(r, fmin(s, fmin(t, 1 - r - s - t)));
			
			if(d > -TINY) /* Point inside test element */
				return(guess);
			else { /* Reset test element to neighbor */	
				if(fabs(r - d) < TINY) 
					guess = CET_MeshElementArray[guess].NeighborIndex[0];
				else if(fabs(s - d) < TINY) 
					guess = CET_MeshElementArray[guess].NeighborIndex[1];
				else if(fabs(t - d) < TINY) 
					guess = CET_MeshElementArray[guess].NeighborIndex[2];
				else if(fabs(d - (1 - r - s - t)) < TINY) 
					guess = CET_MeshElementArray[guess].NeighborIndex[3];
				else 
					FatalError("Indeterminate neighbor in function Get_Element_Local_Search");
				
				if(guess == -1) /* Neighbor not present, could not find element from local search (likely point left domain) */ 
					return(-1);
			}
		}
		else { /* Dimensions == 2 */
			/* Physical coordinates of nodes of the test element */
			x0 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[0]][0];
			y0 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[0]][1];
			x1 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[1]][0];
			y1 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[1]][1];
			x2 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[2]][0];
			y2 = CET_MeshNodeArray[CET_MeshElementArray[guess].Nodes[2]][1];
			
			/* Entries for mapping of physical to natural coordinates of test element */
			a11 = y2 - y0;
			a12 = x0 - x2;
			a21 = y0 - y1;
			a22 = x1 - x0;
			
			/* Determinant of mapping from natural to physical coordinates of test element */
			V = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
			
			/* Natural coordinates of point to be interpolated */
			r = (a11 * (x - x0) + a12 * (y - y0)) / V;
			s = (a21 * (x - x0) + a22 * (y - y0)) / V;
			
			d = fmin(r, fmin(s, 1 - r - s));
			
			if(d > -TINY) /* Point inside test element */
				return(guess);
			else { /* Reset test element to neighbor */	
				if(fabs(r - d) < TINY) 
					guess = CET_MeshElementArray[guess].NeighborIndex[0];
				else if(fabs(s - d) < TINY) 
					guess = CET_MeshElementArray[guess].NeighborIndex[1];
				else if(fabs(d - (1 - r - s)) < TINY) 
					guess = CET_MeshElementArray[guess].NeighborIndex[2];
				else 
					FatalError("Indeterminate neighbor in function Get_Element_Local_Search");
				
				if(guess == -1) /* Neighbor not present, could not find element from local search (point may have exited domain) */ 
					return(-1);
			}
		}
	}
	
}

