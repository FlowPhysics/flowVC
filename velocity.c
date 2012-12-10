/*
 *  velocity.c
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "globals.h"
#include "io.h"
#include "macros.h"
#include "memory.h"
#include "structs.h"
#include "velocity.h"


void AllocateVelFieldData(void) {
  /*
   *	Allocates memory for velocity field data for 2 time frames (moving buffer of data)
   */
	
  printf("\n  Allocating memory for velocity data...");
  if (Data_MeshType == CARTESIAN) {
    Vel_CartVelArray_U = Calloc4D(Vel_CartVelArray_U, 2, Vel_CartMesh.XRes, Vel_CartMesh.YRes, Vel_CartMesh.ZRes);
    Vel_CartVelArray_V = Calloc4D(Vel_CartVelArray_V, 2, Vel_CartMesh.XRes, Vel_CartMesh.YRes, Vel_CartMesh.ZRes);
    Vel_CartVelArray_W = Calloc4D(Vel_CartVelArray_W, 2, Vel_CartMesh.XRes, Vel_CartMesh.YRes, Vel_CartMesh.ZRes);
  }
  else {
    Vel_UnstructVelArray_U = Calloc2D(Vel_UnstructVelArray_U, 2, Vel_MeshNumNodes);
    Vel_UnstructVelArray_V = Calloc2D(Vel_UnstructVelArray_V, 2, Vel_MeshNumNodes);
    Vel_UnstructVelArray_W = Calloc2D(Vel_UnstructVelArray_W, 2, Vel_MeshNumNodes);
  }
    /*** UPDATE HERE ***/
  printf("OK!\n");
	
}


void FreeVelFieldData(void) {
  /*
   *	Deallocates memory for velocity field data for 2 time frames (moving buffer of data)
   */ 
	
  if (Data_MeshType == CARTESIAN) {
    Free4D(Vel_CartVelArray_U, 2, Vel_CartMesh.XRes, Vel_CartMesh.YRes);
    Free4D(Vel_CartVelArray_V, 2, Vel_CartMesh.XRes, Vel_CartMesh.YRes);
    Free4D(Vel_CartVelArray_W, 2, Vel_CartMesh.XRes, Vel_CartMesh.YRes);
  }
  else {
    Free2D(Vel_UnstructVelArray_U, 2);
    Free2D(Vel_UnstructVelArray_V, 2);
    Free2D(Vel_UnstructVelArray_W, 2);		
  }
}



void LoadCartVelDataFrame(int frame) {
	
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
	
  if(Vel_CartVelArray_U == NULL || Vel_CartVelArray_V == NULL || Vel_CartVelArray_W == NULL)
    AllocateVelFieldData(); 
	
  /* Allocate memory for tempptr, used as placeholder when swapping pointers to data as window of loaded velocity data changes */
  if((tempptr = (double ***)malloc(Vel_CartMesh.XRes * sizeof(double **))) == NULL)
    FatalError("Malloc failed for tempptr");
  for(i = 0; i < Vel_CartMesh.XRes; i++) 
    if((tempptr[i] = (double **)malloc(Vel_CartMesh.YRes * sizeof(double *))) == NULL)
      FatalError("Malloc failed for tempptr[%d]", i);
	
  /*** Swap addresses ***/
  /* U */
  for(i = 0; i < Vel_CartMesh.XRes; i++)
    for(j = 0; j < Vel_CartMesh.YRes; j++)
      tempptr[i][j] = Vel_CartVelArray_U[slot1][i][j]; /* Save address of slot1 */
  for(i = 0; i < Vel_CartMesh.XRes; i++)
    for(j = 0; j < Vel_CartMesh.YRes; j++)
      Vel_CartVelArray_U[slot1][i][j] = Vel_CartVelArray_U[slot2][i][j]; /* Change address of slot1 to address of slot2 */
  for(i = 0; i < Vel_CartMesh.XRes; i++)
    for(j = 0; j < Vel_CartMesh.YRes; j++)
      Vel_CartVelArray_U[slot2][i][j] = tempptr[i][j]; /* Change address of slot2 to former address of slot1 */
  /* V */
  for(i = 0; i < Vel_CartMesh.XRes; i++)
    for(j = 0; j < Vel_CartMesh.YRes; j++)
      tempptr[i][j] = Vel_CartVelArray_V[slot1][i][j]; /* Save address of slot1 */
  for(i = 0; i < Vel_CartMesh.XRes; i++)
    for(j = 0; j < Vel_CartMesh.YRes; j++)
      Vel_CartVelArray_V[slot1][i][j] = Vel_CartVelArray_V[slot2][i][j]; /* Change address of slot1 to address of slot2 */
  for(i = 0; i < Vel_CartMesh.XRes; i++)
    for(j = 0; j < Vel_CartMesh.YRes; j++)
      Vel_CartVelArray_V[slot2][i][j] = tempptr[i][j]; /* Change address of slot2 to former address of slot1 */
  /* W */
  for(i = 0; i < Vel_CartMesh.XRes; i++)
    for(j = 0; j < Vel_CartMesh.YRes; j++)
      tempptr[i][j] = Vel_CartVelArray_W[slot1][i][j]; /* Save address of slot1 */
  for(i = 0; i < Vel_CartMesh.XRes; i++)
    for(j = 0; j < Vel_CartMesh.YRes; j++)
      Vel_CartVelArray_W[slot1][i][j] = Vel_CartVelArray_W[slot2][i][j]; /* Change address of slot1 to address of slot2 */
  for(i = 0; i < Vel_CartMesh.XRes; i++)
    for(j = 0; j < Vel_CartMesh.YRes; j++)
      Vel_CartVelArray_W[slot2][i][j] = tempptr[i][j]; /* Change address of slot2 to former address of slot1 */
	
  /*** Read in new time slice of data ***/
  /* Open file */
  sprintf(Data_BinFilePath, "%s%s_vel.%d.bin", Path_Data, Data_InFilePrefix, Data_SuffixTMin + ModVal * Data_SuffixTDelta);
  if((Data_BinFileID = fopen(Data_BinFilePath, "rb")) == NULL) 
    FatalError("Could not open file %s", Data_BinFilePath);
  /* Read time stamp */
  if(fread(&timestamp, sizeof(double), 1, Data_BinFileID) < 1) 
    FatalError("Could not read time stamp from file %s", Data_BinFilePath);
  printf("\nLoading velocity data from %s (time stamp = %g)...", Data_BinFilePath, timestamp);
  fflush(stdout);
  /* Read velocity data */
  for(k = 0; k < Vel_CartMesh.ZRes; k++)
    for(j = 0; j < Vel_CartMesh.YRes; j++)
      for(i = 0; i < Vel_CartMesh.XRes; i++)
	if((fread(&Vel_CartVelArray_U[slot2][i][j][k], sizeof(double), 1, Data_BinFileID) < 1) ||
	   (fread(&Vel_CartVelArray_V[slot2][i][j][k], sizeof(double), 1, Data_BinFileID) < 1) ||
	   (fread(&Vel_CartVelArray_W[slot2][i][j][k], sizeof(double), 1, Data_BinFileID) < 1)) /* Read data from file to address of slot2 */
	  FatalError("Could not load velocity data for index [%d][%d][%d] from %s", i, j, k, Data_BinFilePath);
  /* Close file */
  fclose(Data_BinFileID);
	
  /* Clean up local arrays */
  for(i = 0; i < Vel_CartMesh.XRes; i++) { 
    free(tempptr[i]);
    tempptr[i] = NULL; 
  } 
  free(tempptr);
  tempptr = NULL;
	
  printf("OK!\n");
  fflush(stdout);
	
}

/*** UPDATE HERE, create analogous function to one below for cell data */

void LoadUnstructVelDataFrame(int frame) {
	
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
  if(Vel_UnstructVelArray_U == NULL || Vel_UnstructVelArray_V == NULL || Vel_UnstructVelArray_W == NULL) 
    AllocateVelFieldData();
	
  /* Open binary velocity file for reading */
  sprintf(Data_BinFilePath, "%s%s_vel.%d.bin", Path_Data, Data_InFilePrefix, Data_SuffixTMin + ModVal * Data_SuffixTDelta);
  if((Data_BinFileID = fopen(Data_BinFilePath, "rb")) == NULL) 
    FatalError("Could not open file %s", Data_BinFilePath);
	
  /* Read time stamp */
  if(fread(&timestamp, sizeof(double), 1, Data_BinFileID) < 1) 
    FatalError("Could not read time stamp from file %s", Data_BinFilePath);
	
  printf("\nLoading velocity data from %s_vel.%d.bin (time stamp = %g)...", Data_InFilePrefix, Data_SuffixTMin + ModVal * Data_SuffixTDelta, timestamp);
  fflush(stdout);

  /* Virtual data shift */
  tempptr_U = Vel_UnstructVelArray_U[slot1];
  tempptr_V = Vel_UnstructVelArray_V[slot1];
  tempptr_W = Vel_UnstructVelArray_W[slot1];
  Vel_UnstructVelArray_U[slot1] = Vel_UnstructVelArray_U[slot2];
  Vel_UnstructVelArray_V[slot1] = Vel_UnstructVelArray_V[slot2]; 
  Vel_UnstructVelArray_W[slot1] = Vel_UnstructVelArray_W[slot2];
  Vel_UnstructVelArray_U[slot2] = tempptr_U;
  Vel_UnstructVelArray_V[slot2] = tempptr_V;
  Vel_UnstructVelArray_W[slot2] = tempptr_W;
	
  /* Read data */
  for(i = 0; i < Vel_MeshNumNodes; i++) 
    if(fread(&Vel_UnstructVelArray_U[slot2][i], sizeof(double), 1, Data_BinFileID) < 1 || 
       fread(&Vel_UnstructVelArray_V[slot2][i], sizeof(double), 1, Data_BinFileID) < 1 ||
       fread(&Vel_UnstructVelArray_W[slot2][i], sizeof(double), 1, Data_BinFileID) < 1) 
      FatalError("Could not read UnstructVelArray data from file %s", Data_BinFilePath);
	
  /* Replace no slip condition on surface with inward velocity if Int_NormalFlow is true */
  if(Int_NormalFlow && Particle_Radius < TINY) {
    for(i = 0; i < Vel_SurfaceMeshNumNodes; i++) {
      Vel_UnstructVelArray_U[slot2][Vel_SurfaceMeshNodeIDs[i]] = Int_TimeDirection * Int_NormalFlowScaling * Vel_SurfaceMeshInwardNormals[i][0];
      Vel_UnstructVelArray_V[slot2][Vel_SurfaceMeshNodeIDs[i]] = Int_TimeDirection * Int_NormalFlowScaling * Vel_SurfaceMeshInwardNormals[i][1];
      Vel_UnstructVelArray_W[slot2][Vel_SurfaceMeshNodeIDs[i]] = Int_TimeDirection * Int_NormalFlowScaling * Vel_SurfaceMeshInwardNormals[i][2];
    }
  }
	
  fclose(Data_BinFileID);
  printf("OK!\n");
  fflush(stdout);
	
}


int GetLocalIndex_SurfaceMesh(int global_index) {
  int i;
  for(i = 0; i < Vel_SurfaceMeshNumNodes; i++) 
    if(Vel_SurfaceMeshNodeIDs[i] == global_index)
      return i;
  return -1;
}

void GetVelocity(double tq, LagrangianPoint *pt, double *dXdt) {
	
  /* UPDATE HERE */  
  if(Data_MeshType == CARTESIAN) 
    GetVelocity_Cartesian(tq, pt, dXdt);
  else if(Data_MeshType == UNSTRUCTURED)
    GetVelocity_Unstructured(tq, pt, dXdt);
  else if(Data_MeshType == ANALYTIC)
    GetVelocity_Analytic(tq, pt, dXdt);	
  else
    FatalError("Unrecognized value for Data_MeshType in GetVelocity()");
	
}

void GetVelocity_Analytic(double tq, LagrangianPoint *pt, double *dXdt) {
	
  dXdt[0] = sin(tq)*pt->X[0];
  dXdt[1] = sin(tq)*pt->X[1];
  dXdt[2] = sin(tq)*pt->X[2];
	
}

void GetVelocity_Cartesian(double tq, LagrangianPoint *pt, double *dXdt) {
	
  double xloc, yloc, zloc, tloc;
  int    i, j, k;
  double V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16;
  double f0000, f0001, f0010, f0011, f0100, f0101, f0110, f0111, f1000, f1001, f1010, f1011, f1100, f1101, f1110, f1111;
	

  /* check outside domain */
  if(TestOutsideCartVelDomain(pt->X)) {
    dXdt[0] = 0.0;
    dXdt[1] = 0.0;
    dXdt[2] = 0.0;
    pt->LeftDomain = 1;
  }
  else {  
    /* Set "local" coordinates (relative to space-time voxel) */
    /* t */
    tloc = (tq - Data_LoadedTMin) / (Data_LoadedTMax - Data_LoadedTMin);
    if((tloc > 1 + TINY) || (tloc < 0 - TINY))
      FatalError("tloc must be between 0 and 1");
		
    /* x */
    i = (int)floor((pt->X[0] - Vel_CartMesh.XMin) / Vel_CartMesh.XDelta);
    if(i == (Vel_CartMesh.XRes - 1)) {
      i = Vel_CartMesh.XRes - 2;
      xloc = 1.0;
    }
    else
      xloc = (pt->X[0] - Vel_CartMesh.XMin - i * Vel_CartMesh.XDelta) / Vel_CartMesh.XDelta;
		
    /* y */
    j = (int)floor((pt->X[1] - Vel_CartMesh.YMin) / Vel_CartMesh.YDelta);
    if(j == (Vel_CartMesh.YRes - 1)) {
      j = Vel_CartMesh.YRes - 2;
      yloc = 1.0;
    }
    else
      yloc = (pt->X[1] - Vel_CartMesh.YMin - j * Vel_CartMesh.YDelta) / Vel_CartMesh.YDelta;
			
    /* z */
    if(Dimensions == 3) {
      k = (int)floor((pt->X[2] - Vel_CartMesh.ZMin) / Vel_CartMesh.ZDelta);
      if(k == (Vel_CartMesh.ZRes - 1)) {
	k = Vel_CartMesh.ZRes - 2;
	zloc = 1.0;
      }
      else 
	zloc = (pt->X[2] - Vel_CartMesh.ZMin - k * Vel_CartMesh.ZDelta) / Vel_CartMesh.ZDelta;
    }
    else {
      k = 0;
      zloc = 0.0;
    }
		
    /* Linear Interpolation coefficients */
    V1 = (tloc)*(xloc)*(yloc)*(zloc);
    V2 = (tloc)*(xloc)*(yloc)*(1-zloc);
    V3 = (tloc)*(xloc)*(1-yloc)*(zloc);
    V4 = (tloc)*(xloc)*(1-yloc)*(1-zloc);
    V5 = (tloc)*(1-xloc)*(yloc)*(zloc);
    V6 = (tloc)*(1-xloc)*(yloc)*(1-zloc);
    V7 = (tloc)*(1-xloc)*(1-yloc)*(zloc);
    V8 = (tloc)*(1-xloc)*(1-yloc)*(1-zloc);
    V9 = (1-tloc)*(xloc)*(yloc)*(zloc);
    V10 = (1-tloc)*(xloc)*(yloc)*(1-zloc);
    V11 = (1-tloc)*(xloc)*(1-yloc)*(zloc);
    V12 = (1-tloc)*(xloc)*(1-yloc)*(1-zloc);
    V13 = (1-tloc)*(1-xloc)*(yloc)*(zloc);
    V14 = (1-tloc)*(1-xloc)*(yloc)*(1-zloc);
    V15 = (1-tloc)*(1-xloc)*(1-yloc)*(zloc);
    V16 = (1-tloc)*(1-xloc)*(1-yloc)*(1-zloc);
		
    /* Vertices of space-time voxel for U */
    if(Dimensions == 3) {
      f0000 = Vel_CartVelArray_U[0][i][j][k];
      f0001 = Vel_CartVelArray_U[0][i][j][k+1];
      f0010 = Vel_CartVelArray_U[0][i][j+1][k];
      f0011 = Vel_CartVelArray_U[0][i][j+1][k+1];
      f0100 = Vel_CartVelArray_U[0][i+1][j][k];
      f0101 = Vel_CartVelArray_U[0][i+1][j][k+1];
      f0110 = Vel_CartVelArray_U[0][i+1][j+1][k];
      f0111 = Vel_CartVelArray_U[0][i+1][j+1][k+1];
      f1000 = Vel_CartVelArray_U[1][i][j][k];
      f1001 = Vel_CartVelArray_U[1][i][j][k+1];
      f1010 = Vel_CartVelArray_U[1][i][j+1][k];
      f1011 = Vel_CartVelArray_U[1][i][j+1][k+1];
      f1100 = Vel_CartVelArray_U[1][i+1][j][k];
      f1101 = Vel_CartVelArray_U[1][i+1][j][k+1];
      f1110 = Vel_CartVelArray_U[1][i+1][j+1][k];
      f1111 = Vel_CartVelArray_U[1][i+1][j+1][k+1];
    }
    else {
      f0000 = Vel_CartVelArray_U[0][i][j][k];
      f0001 = Vel_CartVelArray_U[0][i][j][k];
      f0010 = Vel_CartVelArray_U[0][i][j+1][k];
      f0011 = Vel_CartVelArray_U[0][i][j+1][k];
      f0100 = Vel_CartVelArray_U[0][i+1][j][k];
      f0101 = Vel_CartVelArray_U[0][i+1][j][k];
      f0110 = Vel_CartVelArray_U[0][i+1][j+1][k];
      f0111 = Vel_CartVelArray_U[0][i+1][j+1][k];
      f1000 = Vel_CartVelArray_U[1][i][j][k];
      f1001 = Vel_CartVelArray_U[1][i][j][k];
      f1010 = Vel_CartVelArray_U[1][i][j+1][k];
      f1011 = Vel_CartVelArray_U[1][i][j+1][k];
      f1100 = Vel_CartVelArray_U[1][i+1][j][k];
      f1101 = Vel_CartVelArray_U[1][i+1][j][k];
      f1110 = Vel_CartVelArray_U[1][i+1][j+1][k];
      f1111 = Vel_CartVelArray_U[1][i+1][j+1][k];
    }
		
    /* Linear interpolation for U */
    dXdt[0] = f0000*V16 + f0001*V15 + f0010*V14 + f0011*V13 + f0100*V12 + f0101*V11 + f0110*V10
      + f0111*V9 + f1000*V8 + f1001*V7 + f1010*V6 + f1011*V5 + f1100*V4 + f1101*V3 + f1110*V2 + f1111*V1; 
		
    /* Vertices of space-time voxel for V */
    if(Dimensions == 3) {
      f0000 = Vel_CartVelArray_V[0][i][j][k];
      f0001 = Vel_CartVelArray_V[0][i][j][k+1];
      f0010 = Vel_CartVelArray_V[0][i][j+1][k];
      f0011 = Vel_CartVelArray_V[0][i][j+1][k+1];
      f0100 = Vel_CartVelArray_V[0][i+1][j][k];
      f0101 = Vel_CartVelArray_V[0][i+1][j][k+1];
      f0110 = Vel_CartVelArray_V[0][i+1][j+1][k];
      f0111 = Vel_CartVelArray_V[0][i+1][j+1][k+1];
      f1000 = Vel_CartVelArray_V[1][i][j][k];
      f1001 = Vel_CartVelArray_V[1][i][j][k+1];
      f1010 = Vel_CartVelArray_V[1][i][j+1][k];
      f1011 = Vel_CartVelArray_V[1][i][j+1][k+1];
      f1100 = Vel_CartVelArray_V[1][i+1][j][k];
      f1101 = Vel_CartVelArray_V[1][i+1][j][k+1];
      f1110 = Vel_CartVelArray_V[1][i+1][j+1][k];
      f1111 = Vel_CartVelArray_V[1][i+1][j+1][k+1];
    }
    else {
      f0000 = Vel_CartVelArray_V[0][i][j][k];
      f0001 = Vel_CartVelArray_V[0][i][j][k];
      f0010 = Vel_CartVelArray_V[0][i][j+1][k];
      f0011 = Vel_CartVelArray_V[0][i][j+1][k];
      f0100 = Vel_CartVelArray_V[0][i+1][j][k];
      f0101 = Vel_CartVelArray_V[0][i+1][j][k];
      f0110 = Vel_CartVelArray_V[0][i+1][j+1][k];
      f0111 = Vel_CartVelArray_V[0][i+1][j+1][k];
      f1000 = Vel_CartVelArray_V[1][i][j][k];
      f1001 = Vel_CartVelArray_V[1][i][j][k];
      f1010 = Vel_CartVelArray_V[1][i][j+1][k];
      f1011 = Vel_CartVelArray_V[1][i][j+1][k];
      f1100 = Vel_CartVelArray_V[1][i+1][j][k];
      f1101 = Vel_CartVelArray_V[1][i+1][j][k];
      f1110 = Vel_CartVelArray_V[1][i+1][j+1][k];
      f1111 = Vel_CartVelArray_V[1][i+1][j+1][k];
    }
		
    /* Linear interpolation for V */
    dXdt[1] = f0000*V16 + f0001*V15 + f0010*V14 + f0011*V13 + f0100*V12 + f0101*V11 + f0110*V10
      + f0111*V9 + f1000*V8 + f1001*V7 + f1010*V6 + f1011*V5 + f1100*V4 + f1101*V3 + f1110*V2 + f1111*V1; 
		
    /* Vertices of space-time voxel for W */
    if(Dimensions == 3) {
      f0000 = Vel_CartVelArray_W[0][i][j][k];
      f0001 = Vel_CartVelArray_W[0][i][j][k+1];
      f0010 = Vel_CartVelArray_W[0][i][j+1][k];
      f0011 = Vel_CartVelArray_W[0][i][j+1][k+1];
      f0100 = Vel_CartVelArray_W[0][i+1][j][k];
      f0101 = Vel_CartVelArray_W[0][i+1][j][k+1];
      f0110 = Vel_CartVelArray_W[0][i+1][j+1][k];
      f0111 = Vel_CartVelArray_W[0][i+1][j+1][k+1];
      f1000 = Vel_CartVelArray_W[1][i][j][k];
      f1001 = Vel_CartVelArray_W[1][i][j][k+1];
      f1010 = Vel_CartVelArray_W[1][i][j+1][k];
      f1011 = Vel_CartVelArray_W[1][i][j+1][k+1];
      f1100 = Vel_CartVelArray_W[1][i+1][j][k];
      f1101 = Vel_CartVelArray_W[1][i+1][j][k+1];
      f1110 = Vel_CartVelArray_W[1][i+1][j+1][k];
      f1111 = Vel_CartVelArray_W[1][i+1][j+1][k+1];
    }
    else {
      f0000 = Vel_CartVelArray_W[0][i][j][k];
      f0001 = Vel_CartVelArray_W[0][i][j][k];
      f0010 = Vel_CartVelArray_W[0][i][j+1][k];
      f0011 = Vel_CartVelArray_W[0][i][j+1][k];
      f0100 = Vel_CartVelArray_W[0][i+1][j][k];
      f0101 = Vel_CartVelArray_W[0][i+1][j][k];
      f0110 = Vel_CartVelArray_W[0][i+1][j+1][k];
      f0111 = Vel_CartVelArray_W[0][i+1][j+1][k];
      f1000 = Vel_CartVelArray_W[1][i][j][k];
      f1001 = Vel_CartVelArray_W[1][i][j][k];
      f1010 = Vel_CartVelArray_W[1][i][j+1][k];
      f1011 = Vel_CartVelArray_W[1][i][j+1][k];
      f1100 = Vel_CartVelArray_W[1][i+1][j][k];
      f1101 = Vel_CartVelArray_W[1][i+1][j][k];
      f1110 = Vel_CartVelArray_W[1][i+1][j+1][k];
      f1111 = Vel_CartVelArray_W[1][i+1][j+1][k];
    }
		
    /* Linear interpolation for W */
    dXdt[2] = f0000*V16 + f0001*V15 + f0010*V14 + f0011*V13 + f0100*V12 + f0101*V11 + f0110*V10
      + f0111*V9 + f1000*V8 + f1001*V7 + f1010*V6 + f1011*V5 + f1100*V4 + f1101*V3 + f1110*V2 + f1111*V1; 
  }
	
}

/* UPDATE HERE, create analogous function to one below for cell data */

void GetVelocity_Unstructured(const double tq, LagrangianPoint *pt, double *dXdt) {
	
  double x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;
  double a11, a12, a13, a21, a22, a23, a31, a32, a33;
  double U1t, U2t, U3t, U4t, V1t, V2t, V3t, V4t, W1t, W2t, W3t, W4t;
  double r, s, t, d;
  double V;
  double tloc;
	
  if(pt->ElementIndex == -1) 
    FatalError("Attempting to interpolate velocity at point with Element_Index = -1");
	
  x = pt->X[0];
  y = pt->X[1];
  z = pt->X[2];
	
  /* Set tloc defining where in between loaded data to interpolate in time */
  tloc = (tq - Data_LoadedTMin) / (Data_LoadedTMax - Data_LoadedTMin);
  if((tloc > 1 + TINY) || (tloc < 0 - TINY)) 
    FatalError("tloc = %f and must be between 0 and 1", tloc);
	
  while(1) { 
    /* Locate the point in the velocity mesh */
    if(Dimensions == 3) {
      /* Physical coordinates of nodes of the test element */
      x0 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[0]][0];
      y0 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[0]][1];
      z0 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[0]][2];
      x1 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[1]][0];
      y1 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[1]][1];
      z1 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[1]][2];
      x2 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[2]][0];
      y2 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[2]][1];
      z2 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[2]][2];
      x3 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[3]][0];
      y3 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[3]][1];
      z3 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[3]][2];
			
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
    }
    else { /* Dimensions == 2 */
      /* Physical coordinates of nodes of the test element */
      x0 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[0]][0];
      y0 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[0]][1];
      x1 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[1]][0];
      y1 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[1]][1];
      x2 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[2]][0];
      y2 = Vel_MeshNodeArray[Vel_MeshElementArray[pt->ElementIndex].Nodes[2]][1];
			
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
      t = 0;
			
      d = fmin(r, fmin(s, 1 - r - s));
    }
		
    if((d + TINY) >= 0) {   /* Point inside test element */
      /* Interpolate velocity at nodes in time */
      U1t = (1 - tloc) * Vel_UnstructVelArray_U[0][Vel_MeshElementArray[pt->ElementIndex].Nodes[0]] + 
	tloc * Vel_UnstructVelArray_U[1][Vel_MeshElementArray[pt->ElementIndex].Nodes[0]];
      V1t = (1 - tloc) * Vel_UnstructVelArray_V[0][Vel_MeshElementArray[pt->ElementIndex].Nodes[0]] + 
	tloc * Vel_UnstructVelArray_V[1][Vel_MeshElementArray[pt->ElementIndex].Nodes[0]];
      U2t = (1 - tloc) * Vel_UnstructVelArray_U[0][Vel_MeshElementArray[pt->ElementIndex].Nodes[1]] + 
	tloc * Vel_UnstructVelArray_U[1][Vel_MeshElementArray[pt->ElementIndex].Nodes[1]];
      V2t = (1 - tloc) * Vel_UnstructVelArray_V[0][Vel_MeshElementArray[pt->ElementIndex].Nodes[1]] + 
	tloc * Vel_UnstructVelArray_V[1][Vel_MeshElementArray[pt->ElementIndex].Nodes[1]];
      U3t = (1 - tloc) * Vel_UnstructVelArray_U[0][Vel_MeshElementArray[pt->ElementIndex].Nodes[2]] + 
	tloc * Vel_UnstructVelArray_U[1][Vel_MeshElementArray[pt->ElementIndex].Nodes[2]];
      V3t = (1 - tloc) * Vel_UnstructVelArray_V[0][Vel_MeshElementArray[pt->ElementIndex].Nodes[2]] + 
	tloc * Vel_UnstructVelArray_V[1][Vel_MeshElementArray[pt->ElementIndex].Nodes[2]];
      if(Dimensions == 3) {
	W1t = (1 - tloc) * Vel_UnstructVelArray_W[0][Vel_MeshElementArray[pt->ElementIndex].Nodes[0]] + 
	  tloc * Vel_UnstructVelArray_W[1][Vel_MeshElementArray[pt->ElementIndex].Nodes[0]];
	W2t = (1 - tloc) * Vel_UnstructVelArray_W[0][Vel_MeshElementArray[pt->ElementIndex].Nodes[1]] + 
	  tloc * Vel_UnstructVelArray_W[1][Vel_MeshElementArray[pt->ElementIndex].Nodes[1]];       
	W3t = (1 - tloc) * Vel_UnstructVelArray_W[0][Vel_MeshElementArray[pt->ElementIndex].Nodes[2]] + 
	  tloc * Vel_UnstructVelArray_W[1][Vel_MeshElementArray[pt->ElementIndex].Nodes[2]]; 
	U4t = (1 - tloc) * Vel_UnstructVelArray_U[0][Vel_MeshElementArray[pt->ElementIndex].Nodes[3]] + 
	  tloc * Vel_UnstructVelArray_U[1][Vel_MeshElementArray[pt->ElementIndex].Nodes[3]];
	V4t = (1 - tloc) * Vel_UnstructVelArray_V[0][Vel_MeshElementArray[pt->ElementIndex].Nodes[3]] + 
	  tloc * Vel_UnstructVelArray_V[1][Vel_MeshElementArray[pt->ElementIndex].Nodes[3]];
	W4t = (1 - tloc) * Vel_UnstructVelArray_W[0][Vel_MeshElementArray[pt->ElementIndex].Nodes[3]] + 
	  tloc * Vel_UnstructVelArray_W[1][Vel_MeshElementArray[pt->ElementIndex].Nodes[3]];
				
	/* Get interpolate velocity at point in space using time-interpolated velocity at nodes*/
	dXdt[0] = U1t + (U2t - U1t) * fabs(r) + (U3t - U1t) * fabs(s) + (U4t - U1t) * fabs(t);
	dXdt[1] = V1t + (V2t - V1t) * fabs(r) + (V3t - V1t) * fabs(s) + (V4t - V1t) * fabs(t);
	dXdt[2] = W1t + (W2t - W1t) * fabs(r) + (W3t - W1t) * fabs(s) + (W4t - W1t) * fabs(t);
      }
      else { /* Dimensions = 2 */
	/* Get interpolate velocity at point in space using time-interpolated velocity at nodes*/
	dXdt[0] = U1t + (U2t - U1t) * fabs(r) + (U3t - U1t) * fabs(s);
	dXdt[1] = V1t + (V2t - V1t) * fabs(r) + (V3t - V1t) * fabs(s);
	dXdt[2] = 0.0;
      }
      return;
    }
    else { /* Point not inside test element */
      /* Reset test element to one of neighbors */
      if(Dimensions == 3) {
	if(fabs(r - d) < TINY) 
	  pt->ElementIndex = Vel_MeshElementArray[pt->ElementIndex].NeighborIndex[0];
	else if(fabs(s - d)  < TINY) 
	  pt->ElementIndex = Vel_MeshElementArray[pt->ElementIndex].NeighborIndex[1];
	else if(fabs(t - d) < TINY) 
	  pt->ElementIndex = Vel_MeshElementArray[pt->ElementIndex].NeighborIndex[2];
	else if(fabs((1 - r - s - t) - d) < TINY) 
	  pt->ElementIndex = Vel_MeshElementArray[pt->ElementIndex].NeighborIndex[3];
	else 
	  FatalError("Indeterminate neighbor in function GetVelocity_Unstructured()");
      }
      else {
	if(fabs(r - d) < TINY) 
	  pt->ElementIndex = Vel_MeshElementArray[pt->ElementIndex].NeighborIndex[0];
	else if(fabs(s - d) < TINY) 
	  pt->ElementIndex = Vel_MeshElementArray[pt->ElementIndex].NeighborIndex[1];
	else if(fabs(d - (1 - r - s)) < TINY) 
	  pt->ElementIndex = Vel_MeshElementArray[pt->ElementIndex].NeighborIndex[2];
	else 
	  FatalError("Indeterminate neighbor in function GetVelocity_Unstructured()");
      }
      /* If point has left domain, flag it and return with velocity set to zero */
      if(pt->ElementIndex == -1) {
	pt->LeftDomain = 1;
	dXdt[0] = 0.0;
	dXdt[1] = 0.0;
	dXdt[2] = 0.0;
	return;
      }
    }
  }
	
}


int TestOutsideCartVelDomain(double X[3]) {
  /***  
	inside data bounding box --> return 0
	outside data bounding box --> return 1
  ***/ 
	
  int imin, jmin, kmin, imax, jmax, kmax;
	
  imin = (int)floor((X[0] - Vel_CartMesh.XMin) / Vel_CartMesh.XDelta);
  jmin = (int)floor((X[1] - Vel_CartMesh.YMin) / Vel_CartMesh.YDelta);
  if(Dimensions == 3) 
    kmin = (int)floor((X[2] - Vel_CartMesh.ZMin) / Vel_CartMesh.ZDelta);
  else 
    kmin = 0;
  imax = (int)ceil((X[0] - Vel_CartMesh.XMin) / Vel_CartMesh.XDelta);
  jmax = (int)ceil((X[1] - Vel_CartMesh.YMin) / Vel_CartMesh.YDelta);
  if(Dimensions == 3) 
    kmax = (int)ceil((X[2] - Vel_CartMesh.ZMin) / Vel_CartMesh.ZDelta);
  else 
    kmax = 0;
	
  if(imin < 0 || imax > (Vel_CartMesh.XRes - 1) || jmin < 0 || jmax > (Vel_CartMesh.YRes - 1) 
     || kmin < 0 || kmax > (Vel_CartMesh.ZRes - 1)) {
    return 1;
  }
  else 
    return 0;
}
	

int TestOutsideDomain(double point[3]) {
  /***  
	inside data bounding box --> return 0
	outside data bounding box --> return 1
  ***/ 
	
  if(point[0] < (Data_MeshBounds.XMin - TINY) || point[0] > (Data_MeshBounds.XMax + TINY)) 
    return(1);
  else if(point[1] < (Data_MeshBounds.YMin - TINY) || point[1] > (Data_MeshBounds.YMax + TINY)) 
    return(1);
  else if(point[2] < (Data_MeshBounds.ZMin - TINY) || point[2] > (Data_MeshBounds.ZMax + TINY)) 
    return(1);
  else
    return(0);
	
}
