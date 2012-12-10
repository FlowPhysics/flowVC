/*
 *  globals.c
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#include "structs.h"
#include "macros.h"

/* Misc globals */
int			Dimensions = 0;
int			LocalSearchChecking = 0;

/* Globals related to velocity data */
double		****Vel_CartVelArray_U = NULL;
double		****Vel_CartVelArray_V = NULL;
double		****Vel_CartVelArray_W = NULL;
double		****Vel_CartStrainRateArray = NULL;
double      ****Vel_CartVorArray_wx = NULL;
double      ****Vel_CartVorArray_wy = NULL;
double      ****Vel_CartVorArray_wz = NULL;
double		**Vel_UnstructVelArray_U = NULL;
double		**Vel_UnstructVelArray_V = NULL;
double		**Vel_UnstructVelArray_W = NULL;
double		**Vel_UnstructStrainRateArray = NULL;
double		**Vel_UnstructVorArray_wx = NULL;
double		**Vel_UnstructVorArray_wy = NULL;
double		**Vel_UnstructVorArray_wz = NULL;
CartMesh	Vel_CartMesh;
Element		*Vel_MeshElementArray = NULL;
Element		*CET_MeshElementArray = NULL;
double		**Vel_MeshNodeArray = NULL;
double		**CET_MeshNodeArray = NULL;
int			Vel_MeshNumElements = 0;
int			CET_MeshNumElements = 0;
int			Vel_MeshNumNodes = 0;
int			CET_MeshNumNodes = 0;
int			*Vel_MeshElementFlagArray = NULL;
int			*Vel_SurfaceMeshNodeIDs = NULL;
int			Vel_SurfaceMeshNumNodes = 0;
double		**Vel_SurfaceMeshInwardNormals = NULL;
int			Data_MeshType = 0; 
double		Data_TMin = 0;
double		Data_TMax = 0;
int			Data_TRes = 0;
double		Data_TDelta = 0;
int			Data_TPeriodic = 0;
int			Data_SuffixTMin = 0;
int			Data_SuffixTDelta = 0;
int			Data_FirstFrame = 0;
int			Data_LastFrame = 0;
double		Data_LoadedTMin = 0;
double		Data_LoadedTMax = 0;
CartMesh	Data_MeshBounds;

/* Globals related to how much output to generate */
double		Output_TStart = 0;
double		Output_TEnd = 0; 
double		Output_TDelta = 0;
int			Output_TRes = 0;

/* Globals related to FTLE computation */
int			FTLE_Compute = 0;
double		FTLE_IntTLength = 0;
int			FTLE_GenerateMesh = 0;
int			FTLE_ComputeVariation = 0;
int			FTLE_VariationOutFreq = 0;
int			FTLE_NumOutput = 0;
CartMesh	FTLE_CartMesh;
FTLEPoint	***FTLE_MeshPt = NULL;
Launch		*FTLE_Launches = NULL;
double		****FTLE_NewArray = NULL;


/* Globals related to integrator settings */
int			Int_Type = 0;
double		Int_TimeStep = 0;
double		Int_Accuracy = 0;
double		Int_MinTimeStep = 0;
double		Int_MaxTimeStep = 0;
int			Int_TimeDirection = 0;
int			Int_NormalFlow = 0;
double		Int_NormalFlowScaling = 0;
int			Int_Extrapolate = 0;

/* Globals related to particle dynamics */
double		Particle_Density = 0;
double      Particle_Radius = 0;
int         Particle_ICType = 0;
double      Fluid_Viscosity = 0;
double      Fluid_Density = 0;
double      K = 0;
double		R = 0;
double      Gravity[3] = {0,0,0};

/* Globals related to tracer computation */
int             Trace_Compute = 0;
int             Trace_VorticityCompute = 0;
int             Trace_APCompute = 0;
double          Trace_IntTLength = 0;
int             Trace_NumTracers = 0;
int             Trace_GenerateMesh = 0;
int             Trace_RTCompute = 0;
int             Trace_CETCompute = 0;
int             *Trace_CETIndexList = NULL;
int             Trace_CETAuxillaryMesh = 0;
int             Trace_CETSubsteps = 0;
double          Trace_LaunchTimeSpacing = 0;
int             Trace_NumLaunchTimes = 0;
int				Trace_MultipleInFiles = 0;
int             *Trace_NumOutput = NULL;
int             Trace_AlwaysOutput = 0;
int             Trace_InFileFormat = 0;
int             Trace_ReleaseStrategy = 0;
int             Trace_NumReleaseLocations = 0;
double          Trace_ReleaseTMax = 0;
CartMesh        Trace_CartMesh;
Launch          *Trace_Launches = NULL;
LagrangianPoint *Trace_MeshPt = NULL;
ReleaseLocation **Trace_ReleaseList = NULL;
LagrangianPoint *Trace_ReleasePoints = NULL;
CETE            *Trace_CETArray = NULL;
double          *Trace_StrainField = NULL;

/* Globals related to outputting interpolated field data */
int             VelOut_Compute = 0;
int             VelOut_NumPts = 0;
int             VelOut_GenerateMesh = 0;
int             VelOut_InFileFormat = 0;
double          *VelOut_OutputTime = NULL;
int             *VelOut_Complete = NULL;
CartMesh        VelOut_CartMesh;
LagrangianPoint *VelOut_Mesh = NULL;

/* Globally defined file names */
char Data_InFilePrefix[LONGSTRING];
char FTLE_ICFile[LONGSTRING];
char FTLE_OutFilePrefix[LONGSTRING];
char Trace_InFile[LONGSTRING];
char Trace_OutFilePrefix[LONGSTRING];
char TraceOUT_BinFilePath[LONGSTRING];
char TraceIC_BinFilePath[LONGSTRING];
char TraceN_BinFilePathPrefix[LONGSTRING];
char Trace_CETMeshPrefix[LONGSTRING];
char Trace_RTOutFilePrefix[LONGSTRING];
char VelOut_FilePrefix[LONGSTRING];
char VelOut_InFile[LONGSTRING];
char Path_Data[LONGSTRING];
char Path_Output[LONGSTRING];

/* Globally defined arrays */
double UnitSphere[32][3] = { 
    {-1, 0, 0},
    {1, 0, 0},
    {0, -1, 0},
    {0, 1, 0}, 
    {0, 0, -1},  
    {0, 0, 1},
    {0, -0.707106781186547, -0.707106781186547},  
    {0, -0.707106781186547, 0.707106781186547}, 
    {0, 0.707106781186547, -0.707106781186547}, 
    {0, 0.707106781186547, 0.707106781186547}, 
    {-0.577350269189626, -0.577350269189626, -0.577350269189626}, 
    {-0.577350269189626, -0.577350269189626, 0.577350269189626}, 
    {-0.577350269189626, 0.577350269189626, -0.577350269189626},
    {-0.577350269189626, 0.577350269189626, 0.577350269189626}, 
    {0.577350269189626, -0.577350269189626, -0.577350269189626},
    {0.577350269189626, -0.577350269189626, 0.577350269189626},
    {0.577350269189626, 0.577350269189626, -0.577350269189626},
    {0.577350269189626, 0.577350269189626, 0.577350269189626}, 
    {-0.707106781186547, -0.707106781186547, 0}, 
    {-0.707106781186547, 0.707106781186547, 0},
    {-0.707106781186547, 0, -0.707106781186547},  
    {-0.707106781186547, 0, 0.707106781186547}, 
    {0.707106781186547, -0.707106781186547, 0}, 
    {0.707106781186547, 0.707106781186547, 0},
    {0.707106781186547, 0, -0.707106781186547},  
    {0.707106781186547, 0, 0.707106781186547}
};
