/*
 *  globals.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#ifndef INC_GLOBALS_H
#define INC_GLOBALS_H

#include "structs.h"
#include "macros.h"

/* Misc globals */
extern int Dimensions;
extern int LocalSearchChecking;

/* Globals related to velocity data */
extern double		****Vel_CartVelArray_U;
extern double		****Vel_CartVelArray_V;
extern double		****Vel_CartVelArray_W;
extern double		****Vel_CartStrainRateArray;
extern double		****Vel_CartVorArray_wx;
extern double		****Vel_CartVorArray_wy;
extern double		****Vel_CartVorArray_wz;
extern double		**Vel_UnstructVelArray_U;
extern double		**Vel_UnstructVelArray_V;
extern double		**Vel_UnstructVelArray_W;
extern double		**Vel_UnstructStrainRateArray;
extern double		**Vel_UnstructVorArray_wx;
extern double		**Vel_UnstructVorArray_wy;
extern double		**Vel_UnstructVorArray_wz;
extern CartMesh		Vel_CartMesh;
extern Element		*Vel_MeshElementArray;
extern Element		*CET_MeshElementArray;
extern double		**Vel_MeshNodeArray;
extern double		**CET_MeshNodeArray;
extern int			Vel_MeshNumElements;
extern int			CET_MeshNumElements;
extern int			Vel_MeshNumNodes;
extern int			CET_MeshNumNodes;
extern int			*Vel_MeshElementFlagArray;
extern int			*Vel_SurfaceMeshNodeIDs;
extern int			Vel_SurfaceMeshNumNodes;
extern double		**Vel_SurfaceMeshInwardNormals;
extern int			Data_MeshType; 
extern double		Data_TMin;
extern double		Data_TMax;
extern int			Data_TRes;
extern double		Data_TDelta;
extern int			Data_TPeriodic;
extern int			Data_XPeriodic;
extern int			Data_SuffixTMin;
extern int			Data_SuffixTDelta;
extern int			Data_FirstFrame;
extern int			Data_LastFrame;
extern double		Data_LoadedTMin;
extern double		Data_LoadedTMax;
extern CartMesh		Data_MeshBounds;

/* Globals related to how much output to generate */
extern double Output_TStart;
extern double Output_TEnd; 
extern double Output_TDelta;
extern int    Output_TRes;

/* Globals related to FTLE computation */
extern int       FTLE_Compute;
extern double    FTLE_IntTLength;
extern int       FTLE_GenerateMesh;
extern int       FTLE_ComputeVariation;
extern int       FTLE_VariationOutFreq;
extern int       FTLE_NumOutput;
extern CartMesh  FTLE_CartMesh;
extern FTLEPoint ***FTLE_MeshPt;
extern Launch    *FTLE_Launches;
extern double    ****FTLE_NewArray;

/* Globals related to integrator settings */
extern int    Int_Type;
extern double Int_TimeStep;
extern double Int_Accuracy;
extern double Int_MinTimeStep;
extern double Int_MaxTimeStep;
extern int    Int_TimeDirection;
extern int    Int_NormalFlow;
extern double Int_NormalFlowScaling;
extern int    Int_Extrapolate;

/* Globals related to particle dynamics */
extern double          Particle_Density;
extern double          Particle_Radius;
extern int             Particle_ICType;
extern double          Fluid_Viscosity;
extern double          Fluid_Density;
extern double          K, R;
extern double          Gravity[3];

/* Globals related to tracer computation */
extern int             Trace_Compute;
extern int             Trace_VorticityCompute;
extern int             Trace_APCompute;
extern double          Trace_IntTLength;
extern int             Trace_NumTracers;
extern int             Trace_GenerateMesh;
extern int             Trace_RTCompute;
extern int             Trace_CETCompute;
extern int             *Trace_CETIndexList;
extern int             Trace_CETAuxillaryMesh;
extern int             Trace_CETSubsteps;
extern double          Trace_LaunchTimeSpacing;
extern int             Trace_NumLaunchTimes;
extern int			   Trace_MultipleInFiles;
extern int             *Trace_NumOutput;
extern int             Trace_AlwaysOutput;
extern int             Trace_InFileFormat;
extern int             Trace_ReleaseStrategy;
extern int             Trace_NumReleaseLocations;
extern double          Trace_ReleaseTMax;
extern CartMesh        Trace_CartMesh;
extern Launch          *Trace_Launches;
extern LagrangianPoint *Trace_MeshPt;
extern ReleaseLocation **Trace_ReleaseList;
extern LagrangianPoint *Trace_ReleasePoints;
extern CETE            *Trace_CETArray;
extern double          *Trace_StrainField;

/* Globals related to outputting interpolated field data */
extern int             VelOut_Compute;
extern int             VelOut_NumPts;
extern int             VelOut_GenerateMesh;
extern int             VelOut_InFileFormat;
extern double          *VelOut_OutputTime;
extern int             *VelOut_Complete;
extern CartMesh        VelOut_CartMesh;
extern LagrangianPoint *VelOut_Mesh;

/* Globally defined file names */
extern char Data_InFilePrefix[LONGSTRING];
extern char FTLE_ICFile[LONGSTRING];
extern char FTLE_OutFilePrefix[LONGSTRING];
extern char Trace_InFile[LONGSTRING];
extern char Trace_OutFilePrefix[LONGSTRING];
extern char TraceOUT_BinFilePath[LONGSTRING];
extern char TraceIC_BinFilePath[LONGSTRING];
extern char TraceN_BinFilePathPrefix[LONGSTRING];
extern char Trace_CETMeshPrefix[LONGSTRING];
extern char Trace_RTOutFilePrefix[LONGSTRING];
extern char VelOut_FilePrefix[LONGSTRING];
extern char VelOut_InFile[LONGSTRING];
extern char Path_Data[LONGSTRING];
extern char Path_Output[LONGSTRING];

/* Globally defined arrays */
double UnitSphere[32][3];

#endif
