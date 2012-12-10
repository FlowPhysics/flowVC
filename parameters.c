/*
 *  parameters.c
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
#include "macros.h"
#include "mymath.h"
#include "parameters.h"
#include "io.h"

void ReadInParameters(int argc, const char *argv[]) {
	
    FILE *Parameters_InFileID;
	
    if(argc == 2) {
        if((Parameters_InFileID = fopen(argv[1], "r")) == NULL)
            FatalError("Could not open input file %s", argv[1]);
		
        fprintf(stderr, "Reading in parameters...");
		
        ReadInNextValue(Parameters_InFileID, Path_Data, 's');
        ReadInNextValue(Parameters_InFileID, Path_Output, 's');
		
        ReadInNextValue(Parameters_InFileID, &Dimensions, 'd');
		
        ReadInNextValue(Parameters_InFileID, &Data_MeshType, 'd');
        ReadInNextValue(Parameters_InFileID, Data_InFilePrefix, 's');
        ReadInNextValue(Parameters_InFileID, &Data_SuffixTMin, 'd');
        ReadInNextValue(Parameters_InFileID, &Data_SuffixTDelta, 'd');
        ReadInNextValue(Parameters_InFileID, &Data_TRes, 'd');
        ReadInNextValue(Parameters_InFileID, &Data_TDelta, 'f');
        ReadInNextValue(Parameters_InFileID, &Data_TMin, 'f');
        ReadInNextValue(Parameters_InFileID, &Data_TPeriodic, 'd');
        ReadInNextValue(Parameters_InFileID, &Data_MeshBounds.XMin, 'f');
        ReadInNextValue(Parameters_InFileID, &Data_MeshBounds.XMax, 'f');
        ReadInNextValue(Parameters_InFileID, &Data_MeshBounds.YMin, 'f');
        ReadInNextValue(Parameters_InFileID, &Data_MeshBounds.YMax, 'f');
        ReadInNextValue(Parameters_InFileID, &Data_MeshBounds.ZMin, 'f');
        ReadInNextValue(Parameters_InFileID, &Data_MeshBounds.ZMax, 'f');
		
        ReadInNextValue(Parameters_InFileID, &Fluid_Density, 'f');
        ReadInNextValue(Parameters_InFileID, &Fluid_Viscosity, 'f');
		
        ReadInNextValue(Parameters_InFileID, &Output_TStart, 'f');
        ReadInNextValue(Parameters_InFileID, &Output_TRes, 'd');
        ReadInNextValue(Parameters_InFileID, &Output_TDelta, 'f');
		
        ReadInNextValue(Parameters_InFileID, &Int_Type, 'd');
        ReadInNextValue(Parameters_InFileID, &Int_TimeStep, 'f');
        ReadInNextValue(Parameters_InFileID, &Int_Accuracy, 'f');
        ReadInNextValue(Parameters_InFileID, &Int_MinTimeStep, 'f');
        ReadInNextValue(Parameters_InFileID, &Int_MaxTimeStep, 'f');
        ReadInNextValue(Parameters_InFileID, &Int_TimeDirection, 'd');
        ReadInNextValue(Parameters_InFileID, &Int_NormalFlow, 'd');
        ReadInNextValue(Parameters_InFileID, &Int_NormalFlowScaling, 'f');
        ReadInNextValue(Parameters_InFileID, &Int_Extrapolate, 'd');
		
        ReadInNextValue(Parameters_InFileID, &Particle_Radius, 'f');
        ReadInNextValue(Parameters_InFileID, &Particle_Density, 'f');
        ReadInNextValue(Parameters_InFileID, &Particle_ICType, 'd');
		
        ReadInNextValue(Parameters_InFileID, &Gravity[0], 'f');
        ReadInNextValue(Parameters_InFileID, &Gravity[1], 'f');
        ReadInNextValue(Parameters_InFileID, &Gravity[2], 'f');
		
        ReadInNextValue(Parameters_InFileID, &LocalSearchChecking, 'd');
		
        ReadInNextValue(Parameters_InFileID, &FTLE_Compute, 'd');
        ReadInNextValue(Parameters_InFileID, &FTLE_GenerateMesh, 'd');
        ReadInNextValue(Parameters_InFileID, &FTLE_ICFile, 's');
        ReadInNextValue(Parameters_InFileID, &FTLE_CartMesh.XMin, 'f');
        ReadInNextValue(Parameters_InFileID, &FTLE_CartMesh.XMax, 'f');
        ReadInNextValue(Parameters_InFileID, &FTLE_CartMesh.YMin, 'f');
        ReadInNextValue(Parameters_InFileID, &FTLE_CartMesh.YMax, 'f');
        ReadInNextValue(Parameters_InFileID, &FTLE_CartMesh.ZMin, 'f');
        ReadInNextValue(Parameters_InFileID, &FTLE_CartMesh.ZMax, 'f');
        ReadInNextValue(Parameters_InFileID, &FTLE_CartMesh.XRes, 'd');
        ReadInNextValue(Parameters_InFileID, &FTLE_CartMesh.YRes, 'd');
        ReadInNextValue(Parameters_InFileID, &FTLE_CartMesh.ZRes, 'd');
        ReadInNextValue(Parameters_InFileID, &FTLE_IntTLength, 'f');
        ReadInNextValue(Parameters_InFileID, &FTLE_ComputeVariation, 'd');
        ReadInNextValue(Parameters_InFileID, &FTLE_VariationOutFreq, 'd');
        ReadInNextValue(Parameters_InFileID, FTLE_OutFilePrefix, 's');
		
        ReadInNextValue(Parameters_InFileID, &Trace_Compute, 'd');
        ReadInNextValue(Parameters_InFileID, &Trace_ReleaseStrategy, 'd');
        ReadInNextValue(Parameters_InFileID, &Trace_ReleaseTMax, 'f');
        ReadInNextValue(Parameters_InFileID, &Trace_GenerateMesh, 'd');
        ReadInNextValue(Parameters_InFileID, Trace_InFile, 's');
        ReadInNextValue(Parameters_InFileID, &Trace_MultipleInFiles, 'd');
        ReadInNextValue(Parameters_InFileID, &Trace_InFileFormat, 'd');
        ReadInNextValue(Parameters_InFileID, Trace_OutFilePrefix, 's');
        ReadInNextValue(Parameters_InFileID, &Trace_NumLaunchTimes, 'd');
        ReadInNextValue(Parameters_InFileID, &Trace_LaunchTimeSpacing, 'f');
        ReadInNextValue(Parameters_InFileID, &Trace_IntTLength, 'f');
        ReadInNextValue(Parameters_InFileID, &Trace_AlwaysOutput, 'd');
        ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.XMin, 'f');
        ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.XMax, 'f');
        ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.YMin, 'f');
        ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.YMax, 'f');
        ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.ZMin, 'f');
        ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.ZMax, 'f');
        ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.XRes, 'd');
        ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.YRes, 'd');
        ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.ZRes, 'd');
        ReadInNextValue(Parameters_InFileID, &Trace_VorticityCompute, 'd');
        ReadInNextValue(Parameters_InFileID, &Trace_APCompute, 'd');
        ReadInNextValue(Parameters_InFileID, &Trace_CETCompute, 'd');
        ReadInNextValue(Parameters_InFileID, &Trace_CETAuxillaryMesh, 'd');
        ReadInNextValue(Parameters_InFileID, &Trace_CETMeshPrefix, 's');
        ReadInNextValue(Parameters_InFileID, &Trace_CETSubsteps, 'd');
        ReadInNextValue(Parameters_InFileID, &Trace_RTCompute, 'd');
        ReadInNextValue(Parameters_InFileID, Trace_RTOutFilePrefix, 's');
		
        ReadInNextValue(Parameters_InFileID, &VelOut_Compute, 'd');
        ReadInNextValue(Parameters_InFileID, &VelOut_GenerateMesh, 'd');
        ReadInNextValue(Parameters_InFileID, VelOut_InFile, 's');
        ReadInNextValue(Parameters_InFileID, &VelOut_InFileFormat, 'd');
        ReadInNextValue(Parameters_InFileID, VelOut_FilePrefix, 's');
        ReadInNextValue(Parameters_InFileID, &VelOut_CartMesh.XMin, 'f');
        ReadInNextValue(Parameters_InFileID, &VelOut_CartMesh.XMax, 'f');
        ReadInNextValue(Parameters_InFileID, &VelOut_CartMesh.YMin, 'f');
        ReadInNextValue(Parameters_InFileID, &VelOut_CartMesh.YMax, 'f');
        ReadInNextValue(Parameters_InFileID, &VelOut_CartMesh.ZMin, 'f');
        ReadInNextValue(Parameters_InFileID, &VelOut_CartMesh.ZMax, 'f');
        ReadInNextValue(Parameters_InFileID, &VelOut_CartMesh.XRes, 'd');
        ReadInNextValue(Parameters_InFileID, &VelOut_CartMesh.YRes, 'd');
        ReadInNextValue(Parameters_InFileID, &VelOut_CartMesh.ZRes, 'd');
		
        fclose(Parameters_InFileID);
		
        fprintf(stderr,"Path_Data\t\t=\t%10s\n",Path_Data);
        fprintf(stderr,"Path_Output\t\t=\t%10s\n\n",Path_Output);
		
        fprintf(stderr,"Dimensions\t\t=\t%10d\n", Dimensions);
		
        fprintf(stderr,"Data_MeshType\t\t=\t%10d\n",Data_MeshType);
        fprintf(stderr,"Data_InFilePrefix\t=\t%10s\n",Data_InFilePrefix);
        fprintf(stderr,"Data_SuffixTMin\t\t=\t%10d\n",Data_SuffixTMin);
        fprintf(stderr,"Data_SuffixTDelta\t=\t%10d\n",Data_SuffixTDelta);
        fprintf(stderr,"Data_TRes\t\t=\t%10d\n",Data_TRes);
        fprintf(stderr,"Data_TDelta\t\t=\t%10g\n",Data_TDelta);
        fprintf(stderr,"Data_TMin\t\t=\t%10g\n",Data_TMin);
        fprintf(stderr,"Data_TPeriodic\t\t=\t%10d\n",Data_TPeriodic);
        fprintf(stderr,"Data_MeshBounds.XMin\t=\t%10g\n",Data_MeshBounds.XMin);
        fprintf(stderr,"Data_MeshBounds.XMax\t=\t%10g\n",Data_MeshBounds.XMax);
        fprintf(stderr,"Data_MeshBounds.YMin\t=\t%10g\n",Data_MeshBounds.YMin);
        fprintf(stderr,"Data_MeshBounds.YMax\t=\t%10g\n",Data_MeshBounds.YMax);
        fprintf(stderr,"Data_MeshBounds.ZMin\t=\t%10g\n",Data_MeshBounds.ZMin);
        fprintf(stderr,"Data_MeshBounds.ZMax\t=\t%10g\n\n",Data_MeshBounds.ZMax);
		
        fprintf(stderr,"Fluid_Density\t\t=\t%10g\n",Fluid_Density);
        fprintf(stderr,"Fluid_Viscosity\t\t=\t%10g\n",Fluid_Viscosity);
		
        fprintf(stderr,"Output_TStart\t\t=\t%10g\n",Output_TStart);
        fprintf(stderr,"Output_TRes\t\t=\t%10d\n",Output_TRes);
        fprintf(stderr,"Output_TDelta\t\t=\t%10g\n\n",Output_TDelta);
		
        fprintf(stderr,"Int_Type\t\t=\t%10d\n",Int_Type);
        fprintf(stderr,"Int_TimeStep\t\t=\t%10g\n",Int_TimeStep);
        fprintf(stderr,"Int_Accuracy\t\t=\t%10g\n",Int_Accuracy);
        fprintf(stderr,"Int_MinTimeStep\t\t=\t%10g\n",Int_MinTimeStep);
        fprintf(stderr,"Int_MaxTimeStep\t\t=\t%10g\n",Int_MaxTimeStep);
        fprintf(stderr,"Int_TimeDirection\t=\t%10d\n",Int_TimeDirection);
        fprintf(stderr,"Int_NormalFlow\t\t=\t%10d\n",Int_NormalFlow);
        fprintf(stderr,"Int_NormalFlowScaling\t=\t%10g\n",Int_NormalFlowScaling);
        fprintf(stderr,"Int_Extrapolate\t\t=\t%10d\n\n", Int_Extrapolate);
		
        fprintf(stderr,"Particle_Radius\t\t=\t%10f\n",Particle_Radius);
        fprintf(stderr,"Particle_Density\t=\t%10f\n",Particle_Density);
        fprintf(stderr,"Particle_ICflag\t\t=\t%10d\n",Particle_ICType);
		
        fprintf(stderr,"Gravity[0]\t\t=\t%10f\n",Gravity[0]);
        fprintf(stderr,"Gravity[1]\t\t=\t%10f\n",Gravity[1]);
        fprintf(stderr,"Gravity[2]\t\t=\t%10f\n",Gravity[2]);
		
        fprintf(stderr,"LocalSearchChecking\t=\t%10d\n\n",LocalSearchChecking);
		
        fprintf(stderr,"FTLE_Compute\t\t=\t%10d\n",FTLE_Compute);
        fprintf(stderr,"FTLE_GenerateMesh\t=\t%10d\n",FTLE_GenerateMesh);
        fprintf(stderr,"FTLE_ICFile\t\t=\t%s\n",FTLE_ICFile);
        fprintf(stderr,"FTLE_CartMesh.XMin\t=\t%10g\n",FTLE_CartMesh.XMin);
        fprintf(stderr,"FTLE_CartMesh.XMax\t=\t%10g\n",FTLE_CartMesh.XMax);
        fprintf(stderr,"FTLE_CartMesh.YMin\t=\t%10g\n",FTLE_CartMesh.YMin);
        fprintf(stderr,"FTLE_CartMesh.YMax\t=\t%10g\n",FTLE_CartMesh.YMax);
        fprintf(stderr,"FTLE_CartMesh.ZMin\t=\t%10g\n",FTLE_CartMesh.ZMin);
        fprintf(stderr,"FTLE_CartMesh.ZMax\t=\t%10g\n",FTLE_CartMesh.ZMax);
        fprintf(stderr,"FTLE_CartMesh.XRes\t=\t%10d\n",FTLE_CartMesh.XRes);
        fprintf(stderr,"FTLE_CartMesh.YRes\t=\t%10d\n",FTLE_CartMesh.YRes);
        fprintf(stderr,"FTLE_CartMesh.ZRes\t=\t%10d\n",FTLE_CartMesh.ZRes);
        fprintf(stderr,"FTLE_IntTLength\t\t=\t%10g\n",FTLE_IntTLength);
        fprintf(stderr,"FTLE_ComputeVariation\t=\t%10d\n",FTLE_ComputeVariation);
        fprintf(stderr,"FTLE_VariationOutFreq\t=\t%10d\n",FTLE_VariationOutFreq);
        fprintf(stderr,"FTLE_OutFilePrefix\t=\t%10s\n\n",FTLE_OutFilePrefix);
		
        fprintf(stderr,"Trace_Compute\t\t=\t%10d\n",Trace_Compute);
        fprintf(stderr,"Trace_ReleaseStrategy\t=\t%10d\n",Trace_ReleaseStrategy);
        fprintf(stderr,"Trace_ReleaseTMax\t=\t%10g\n",Trace_ReleaseTMax);
        fprintf(stderr,"Trace_GenerateMesh\t=\t%10d\n",Trace_GenerateMesh);
        fprintf(stderr,"Trace_InFile\t\t=\t%10s\n",Trace_InFile);
        fprintf(stderr,"Trace_MultipleInFiles\t=\t%10d\n",Trace_MultipleInFiles);
        fprintf(stderr,"Trace_InFileFormat\t=\t%10d\n",Trace_InFileFormat);
        fprintf(stderr,"Trace_OutFilePrefix\t=\t%10s\n",Trace_OutFilePrefix);
        fprintf(stderr,"Trace_NumLaunchTimes\t=\t%10d\n",Trace_NumLaunchTimes);
        fprintf(stderr,"Trace_LaunchTimeSpacing\t=\t%10g\n",Trace_LaunchTimeSpacing);
        fprintf(stderr,"Trace_IntTLength\t=\t%10g\n",Trace_IntTLength);
        fprintf(stderr,"Trace_AlwaysOutput\t=\t%10d\n",Trace_AlwaysOutput);
        fprintf(stderr,"Trace_CartMesh.XMin\t=\t%10g\n",Trace_CartMesh.XMin);
        fprintf(stderr,"Trace_CartMesh.XMax\t=\t%10g\n",Trace_CartMesh.XMax);
        fprintf(stderr,"Trace_CartMesh.YMin\t=\t%10g\n",Trace_CartMesh.YMin);
        fprintf(stderr,"Trace_CartMesh.YMax\t=\t%10g\n",Trace_CartMesh.YMax);
        fprintf(stderr,"Trace_CartMesh.ZMin\t=\t%10g\n",Trace_CartMesh.ZMin);
        fprintf(stderr,"Trace_CartMesh.ZMax\t=\t%10g\n",Trace_CartMesh.ZMax);
        fprintf(stderr,"Trace_CartMesh.XRes\t=\t%10d\n",Trace_CartMesh.XRes);
        fprintf(stderr,"Trace_CartMesh.YRes\t=\t%10d\n",Trace_CartMesh.YRes);
        fprintf(stderr,"Trace_CartMesh.ZRes\t=\t%10d\n",Trace_CartMesh.ZRes);
        fprintf(stderr,"Trace_VorticityCompute\t=\t%10d\n",Trace_VorticityCompute);
        fprintf(stderr,"Trace_APCompute\t\t=\t%10d\n",Trace_APCompute);
        fprintf(stderr,"Trace_CETCompute\t=\t%10d\n",Trace_CETCompute);
        fprintf(stderr,"Trace_CETAuxillaryMesh\t=\t%10d\n",Trace_CETAuxillaryMesh);
        fprintf(stderr,"Trace_CETMeshPrefix\t=\t%10s\n",Trace_CETMeshPrefix);
        fprintf(stderr,"Trace_CETSubsteps\t=\t%10d\n",Trace_CETSubsteps);
        fprintf(stderr,"Trace_RTCompute\t\t=\t%10d\n",Trace_RTCompute);
        fprintf(stderr,"Trace_RTOutFilePrefix\t=\t%10s\n\n",Trace_RTOutFilePrefix);
		
        fprintf(stderr,"VelOut_Compute\t\t=\t%10d\n",VelOut_Compute);
        fprintf(stderr,"VelOut_GenerateMesh\t=\t%10d\n",VelOut_GenerateMesh);
        fprintf(stderr,"VelOut_InFile\t\t=\t%10s\n",VelOut_InFile);
        fprintf(stderr,"VelOut_InFileFormat\t=\t%10d\n",VelOut_InFileFormat);
        fprintf(stderr,"VelOut_FilePrefix\t=\t%10s\n",VelOut_FilePrefix);
        fprintf(stderr,"VelOut_CartMesh.XMin\t=\t%10g\n",VelOut_CartMesh.XMin);
        fprintf(stderr,"VelOut_CartMesh.XMax\t=\t%10g\n",VelOut_CartMesh.XMax);
        fprintf(stderr,"VelOut_CartMesh.YMin\t=\t%10g\n",VelOut_CartMesh.YMin);
        fprintf(stderr,"VelOut_CartMesh.YMax\t=\t%10g\n",VelOut_CartMesh.YMax);
        fprintf(stderr,"VelOut_CartMesh.ZMin\t=\t%10g\n",VelOut_CartMesh.ZMin);
        fprintf(stderr,"VelOut_CartMesh.ZMax\t=\t%10g\n",VelOut_CartMesh.ZMax);
        fprintf(stderr,"VelOut_CartMesh.XRes\t=\t%10d\n",VelOut_CartMesh.XRes);
        fprintf(stderr,"VelOut_CartMesh.YRes\t=\t%10d\n",VelOut_CartMesh.YRes);
        fprintf(stderr,"VelOut_CartMesh.ZRes\t=\t%10d\n\n",VelOut_CartMesh.ZRes);
    }
    else {
        fprintf(stderr, "usage: %s inputfile\n", argv[0]);
        exit(1);
    }
}

void ReadInNextValue(FILE *Parameters_InFileID, void *pt, char type) {
	
    int c;
    char buf[LONGSTRING];
	
    while(1) {
        c = fgetc(Parameters_InFileID);
		
        if(c == ' ' || c == '\n' || c == '\r') {
            continue;
        }
		
        if(c == '#') {
            /* Comment line, read until end of line */
            while(c != '\n') {
                c = fgetc(Parameters_InFileID);
            }
        }
        else {
            /* Line should contain variable information of form VARIABLE = VALUE*/
            ungetc(c, Parameters_InFileID);
            fscanf(Parameters_InFileID, "%s %*s", buf);
            if(type == 's')
                fscanf(Parameters_InFileID,"%s\n", (char *) pt);
            else if(type == 'd')
                fscanf(Parameters_InFileID,"%d\n", (int *) pt);
            else if(type == 'f')
                fscanf(Parameters_InFileID,"%lf\n", (double *) pt);
            else
                FatalError("Unsupported data type passed to ReadInNextValue().\n");
            break;
        }
    }
}

void SetDerivedParameters(void) {
	
    double Data_TMinReq;
    double Data_TMaxReq;
	
    printf("Setting derived parameters...\n");
	
    Data_TMax = Data_TMin + (Data_TRes - 1) * Data_TDelta;
    if((Output_TStart < Data_TMin || Output_TStart > Data_TMax) && !Data_TPeriodic)
        FatalError("Output_TStart outside data range");
    Output_TEnd = Output_TStart + (Output_TRes - 1) * Int_TimeDirection * Output_TDelta;
    if((Output_TEnd > Data_TMax || Output_TEnd < Data_TMin) && !Data_TPeriodic)
        FatalError("Output_TEnd outside data range");
    else
        printf("  Output_TEnd\t\t=\t%10g\n", Output_TEnd);
	
    if(FTLE_Compute) {
        FTLE_CartMesh.XDelta = (FTLE_CartMesh.XMax - FTLE_CartMesh.XMin) / (FTLE_CartMesh.XRes - 1);
        FTLE_CartMesh.YDelta = (FTLE_CartMesh.YMax - FTLE_CartMesh.YMin) / (FTLE_CartMesh.YRes - 1);
        if(Dimensions == 3)
            FTLE_CartMesh.ZDelta = (FTLE_CartMesh.ZMax - FTLE_CartMesh.ZMin) / (FTLE_CartMesh.ZRes - 1);
        else {
            FTLE_CartMesh.ZMin = 0;
            FTLE_CartMesh.ZMax = 0;
            FTLE_CartMesh.ZDelta = 0;
            FTLE_CartMesh.ZRes = 1;
        }
		
        Data_TMinReq = fmin(Output_TStart, Output_TEnd + Int_TimeDirection * FTLE_IntTLength);
        if((Data_TMinReq < Data_TMin) && !Data_TPeriodic) {
            fprintf(stderr, "\nWARNING: Not enough data to accomodate requested FTLE integration time for all slides.\n");
            Data_TMinReq = Data_TMin;
        }
        Data_TMaxReq = fmax(Output_TStart, Output_TEnd + Int_TimeDirection * FTLE_IntTLength);
        if((Data_TMaxReq > Data_TMax) && !Data_TPeriodic) {
            fprintf(stderr, "\nWARNING: Not enough data to accomodate requested FTLE integration time for all slides.\n");
            Data_TMaxReq = Data_TMax;
        }
        Data_FirstFrame = (Int_TimeDirection > 0) ? floor((Data_TMinReq - Data_TMin) / Data_TDelta) : ceil((Data_TMaxReq - Data_TMin) / Data_TDelta);
        Data_LastFrame = (Int_TimeDirection > 0) ? ceil((Data_TMaxReq - Data_TMin) / Data_TDelta) : floor((Data_TMinReq - Data_TMin) / Data_TDelta);
    }
    else {
        Data_FirstFrame = (Int_TimeDirection > 0) ? floor((Output_TStart - Data_TMin) / Data_TDelta) : ceil((Output_TStart - Data_TMin) / Data_TDelta);
        Data_LastFrame = (Int_TimeDirection > 0) ? ceil((Output_TEnd - Data_TMin) / Data_TDelta) : floor((Output_TEnd - Data_TMin) / Data_TDelta);
    }
	
    printf("  Data_FirstFrame\t=\t%10d\n", Data_FirstFrame);
    printf("  Data_LastFrame\t=\t%10d\n", Data_LastFrame);
	
    if(Trace_GenerateMesh) {
        if(Trace_CartMesh.XRes == 1)
            Trace_CartMesh.XDelta = 0.0;
        else
            Trace_CartMesh.XDelta = (Trace_CartMesh.XMax - Trace_CartMesh.XMin) / ((double)Trace_CartMesh.XRes - 1.0);
        if(Trace_CartMesh.YRes == 1)
            Trace_CartMesh.YDelta = 0.0;
        else
            Trace_CartMesh.YDelta = (Trace_CartMesh.YMax - Trace_CartMesh.YMin) / ((double)Trace_CartMesh.YRes - 1.0);
        if(Trace_CartMesh.ZRes == 1)
            Trace_CartMesh.ZDelta = 0.0;
        else {
            if(Dimensions == 3)
                Trace_CartMesh.ZDelta = (Trace_CartMesh.ZMax - Trace_CartMesh.ZMin) / ((double)Trace_CartMesh.ZRes - 1.0);
            else {
                Trace_CartMesh.ZMin = 0;
                Trace_CartMesh.ZMax = 0;
                Trace_CartMesh.ZDelta = 0;
                fprintf(stderr, "\nWarning: Trace_CartMesh.ZRes being set to 1\n");
                Trace_CartMesh.ZRes = 1;
            }
        }
    }
	
    if(strcmp(Path_Output, "pwd")) {
        if(Path_Output[strlen(Path_Output) - 1] != '/')
            sprintf(Path_Output, "%s/", Path_Output);
    }
    else
        sprintf(Path_Output, "./");
	
    if(strcmp(Path_Data, "pwd")) {
        if(Path_Data[strlen(Path_Data) - 1] != '/' )
            sprintf(Path_Data, "%s/", Path_Data);
    }
    else
        sprintf(Path_Data, "./");
	
    if(Trace_Compute) {
        sprintf(TraceOUT_BinFilePath, "%s%s_OUT.bin", Path_Data, Trace_OutFilePrefix);
        sprintf(TraceIC_BinFilePath, "%s%s.IC", Path_Data, Trace_OutFilePrefix);
        sprintf(TraceN_BinFilePathPrefix, "%s%s_slide", Path_Data, Trace_OutFilePrefix);
    }
	
    if(Particle_Radius < 0)
        FatalError("Particle_Radius cannot be negative");
    else if(Particle_Radius > TINY) {
        K = (4.5 * Fluid_Viscosity) / (Fluid_Density * Particle_Radius * Particle_Radius);
        R = (2 * Fluid_Density) / (Fluid_Density + 2 * Particle_Density);
        printf("K = %f, R = %f\n", K, R);
    }
	
    printf("OK!\n");
	
}

void CheckParameters(void) {
	
    printf("Checking parameters...\n");
	
    if(Dimensions != 2 && Dimensions != 3)
        FatalError("Dimensions must be equal to 2 or 3");
    if(Data_MeshType != CARTESIAN && Data_MeshType != UNSTRUCTURED)
        FatalError("Unsupported Data_MeshType");
    if(Data_TRes < 2)
        FatalError("Data_TRes must be 2 or greater.");
    if(Data_TDelta <= 0)
        FatalError("Data_TDelta <= 0");
    if(Data_TPeriodic != 0 && Data_TPeriodic != 1)
        FatalError("Data_TPeriodic should be 0 or 1");
    else if(Data_TPeriodic)
        printf("  Data specified as periodic: Ensure first and last data files correspond to same point in cycle.\n");
    if(Data_MeshBounds.XMin >= Data_MeshBounds.XMax)
        FatalError("Data_MeshBounds.XMin >= Data_MeshBounds.XMax");
    if(Data_MeshBounds.YMin >= Data_MeshBounds.YMax)
        FatalError("Data_MeshBounds.YMin >= Data_MeshBounds.YMax");
    if(Data_MeshBounds.ZMin > Data_MeshBounds.ZMax)
        FatalError("Data_MeshBounds.ZMin > Data_MeshBounds.ZMax");
    if(Output_TRes < 1)
        FatalError("Output_TRes < 1");
    if(Output_TDelta < TINY && Output_TRes > 1)
        FatalError("Output_TDelta must be positive");
    if(Int_Type != 0 && Int_Type != 1 && Int_Type != 2)
        fprintf(stderr, "Warning: Unrecognized value for Int_Type, setting to 0\n");
    if(Int_TimeStep < TINY && Int_Type != 2)
        FatalError("Int_TimeStep must be positive");
    if((Int_MinTimeStep < TINY || Int_MaxTimeStep < TINY || Int_Accuracy < TINY) && Int_Type == 2)
        FatalError("RKF integration parameters must be positive");
    if(Int_TimeDirection != 1 && Int_TimeDirection != -1)
        FatalError("Unsupported Int_TimeDirection value, must be -1 or 1");
    if(Int_TimeDirection == -1 && Particle_Radius > TINY)
        FatalError("Cannot integrate finite size particles backward in time");
    if(Particle_Radius > TINY)
        if(fabs(Int_NormalFlowScaling) > 1)
            fprintf(stderr, "Warning: Int_NormalFlowScaling is restitution coefficient for finite sized particles and should be less than 1.\n");
    
    /*
     if(Int_NormalFlow && Int_NormalFlowScaling < 0)
     FatalError("Int_NormalFlowScaling < 0");
     */
    if(FTLE_Compute) {
        if(Trace_Compute)
            FatalError("Both FTLE_Compute = 1 and Trace_Compute = 1 cannot be handled concurrently");
        if(FTLE_CartMesh.XMin >= FTLE_CartMesh.XMax)
            FatalError("FTLE_CartMesh.XMin >= FTLE_CartMesh.XMax");
        if(FTLE_CartMesh.YMin >= FTLE_CartMesh.YMax)
            FatalError("FTLE_CartMesh.YMin >= FTLE_CartMesh.YMax");
        if(FTLE_CartMesh.ZMin > FTLE_CartMesh.ZMax)
            FatalError("FTLE_CartMesh.ZMin > FTLE_CartMesh.ZMax");
        if(FTLE_CartMesh.XRes < 2)
            FatalError("FTLE_CartMesh.XRes < 2");
        if(FTLE_CartMesh.YRes < 2)
            FatalError("FTLE_CartMesh.YRes < 2");
        if(FTLE_CartMesh.ZRes < 2)
            if((Dimensions == 2 && FTLE_CartMesh.ZRes) < 1 || Dimensions == 3)
                FatalError("FTLE_CartMesh.ZRes insufficient");
        if(FTLE_IntTLength < TINY)
            FatalError("FTLE_IntTLength must be positive");
        if(FTLE_ComputeVariation) {
            if(Output_TRes != 1)
                FatalError("Output_TRes has to be 1 if FTLE_ComputeVariation = 1");
            if(FTLE_VariationOutFreq < 1)
                FatalError("FTLE_VariationOutFreq < 1");
        }
    }
    if(Trace_Compute) {
        if(FTLE_Compute)
            FatalError("Both FTLE_Compute = 1 and Trace_Compute = 1 cannot be handled concurrently");
        if(Trace_ReleaseStrategy == STAGGERED) {
            if(Trace_APCompute)
                fprintf(stderr, "\nWarning: Trace_APCompute not appropriate for Trace_ReleaseStrategy = %d\n\n\a", STAGGERED);
            if(Int_TimeDirection < 0)
                FatalError("Staggered tracer release not supported for backward time integration");
            if(Trace_AlwaysOutput && Trace_IntTLength < (Output_TEnd - Output_TStart)) {
                fprintf(stderr, "\nWarning: Trace_IntTLength should be set larger than Output_TEnd - Output_TStart for Trace_AlwaysOutput = 1\n");
                fprintf(stderr, "\t\t Resetting Trace_IntTLength = %f\n\n\a", Output_TEnd - Output_TStart);
                Trace_IntTLength = Output_TEnd - Output_TStart + TINY;
            }
            if(Trace_ReleaseTMax < 0)
                FatalError("Trace_ReleaseTMax < 0");
            if(Dimensions == 2)
                FatalError("Staggered release not supported for 2D analysis");
            if(Trace_InFileFormat != BINARY_LIST)
                FatalError("Trace_InFileFormat must be %d for Trace_ReleaseStrategy = %d", BINARY_LIST, STAGGERED);
        }
        if(Trace_GenerateMesh) {
            if(Trace_CartMesh.XMin > Trace_CartMesh.XMax)
                FatalError("Trace_CartMesh.XMin > Trace_CartMesh.XMax");
            if(Trace_CartMesh.YMin > Trace_CartMesh.YMax)
                FatalError("Trace_CartMesh.YMin > Trace_CartMesh.YMax");
            if(Trace_CartMesh.ZMin > Trace_CartMesh.ZMax)
                FatalError("Trace_CartMesh.ZMin > Trace_CartMesh.ZMax");
            if(Trace_CartMesh.XRes < 1)
                FatalError("Trace_CartMesh.XRes < 1");
            if(Trace_CartMesh.YRes < 1)
                FatalError("Trace_CartMesh.YRes < 1");
            if(Trace_CartMesh.ZRes < 1)
                FatalError("Trace_CartMesh.ZRes < 1");
        }
        else if(Trace_InFileFormat != INITIALIZED && Trace_InFileFormat != ASCII_LIST &&
                Trace_InFileFormat != VTK_POLYDATA && Trace_InFileFormat != VTK_UNSTRUCTURED && Trace_InFileFormat != BINARY_LIST)
            FatalError("Unrecognized value for Trace_InFileFormat");
        if(Trace_NumLaunchTimes < 1)
            FatalError("Trace_NumLaunchTimes < 1");
        if(Trace_LaunchTimeSpacing < TINY && Trace_NumLaunchTimes > 1)
            FatalError("Trace_LaunchTimeSpacing must be positive");
        if(Trace_IntTLength < 0)
            FatalError("Trace_IntTLength < 0");
        if(Trace_CETCompute) {
            if(Trace_RTCompute)
                FatalError("Cannot compute exposure times and residence times concurrently");
            if(Trace_CETSubsteps < 2)
                FatalError("Trace_CETSubsteps must be greater than 1");
        }
        if((Trace_RTCompute || Trace_CETCompute) && Int_TimeDirection < 0)
            FatalError("Residence/Exposure time computations not supported for backward time integrations");
        if(Trace_InFileFormat == INITIALIZED && Trace_RTCompute && !Trace_GenerateMesh) 
            fprintf(stderr, "  Warning: Trace_InFileFormat = %d not recommended when Trace_RTCompute = 1 as connectivity data will be lost.\n", INITIALIZED);
    }   
    if(VelOut_Compute) {
        if(VelOut_GenerateMesh) {
            if(VelOut_CartMesh.XMin > VelOut_CartMesh.XMax)
                FatalError("VelOut_CartMesh.XMin > VelOut_CartMesh.XMax");
            if(VelOut_CartMesh.YMin > VelOut_CartMesh.YMax)
                FatalError("VelOut_CartMesh.YMin > VelOut_CartMesh.YMax");
            if(VelOut_CartMesh.ZMin > VelOut_CartMesh.ZMax)
                FatalError("VelOut_CartMesh.ZMin > VelOut_CartMesh.ZMax");
            if(VelOut_CartMesh.XRes < 1)
                FatalError("VelOut_CartMesh.XRes < 1");
            if(VelOut_CartMesh.YRes < 1)
                FatalError("VelOut_CartMesh.YRes < 1");
            if(VelOut_CartMesh.ZRes < 1)
                FatalError("VelOut_CartMesh.ZRes < 1");
        }
        if(VelOut_InFileFormat != ASCII_LIST && VelOut_InFileFormat != VTK_POLYDATA && VelOut_InFileFormat != VTK_UNSTRUCTURED)
            FatalError("Unsupported value for VelOut_InFileFormat");
    }
	
    printf("OK!\n\n");
    fflush(stdout);
	
}

