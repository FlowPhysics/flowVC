/*
 *  main.c
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Illinois Institute of Technology. All rights reserved.
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "exposuretime.h"
#include "ftle.h"
#include "globals.h"
#include "integration.h"
#include "io.h"
#include "macros.h"
#include "memory.h"
#include "mesh.h"
#include "parameters.h"
#include "residencetime.h"
#include "structs.h"
#include "strainrate.h"
#include "tracers.h"
#include "velocity.h"
#include "velout.h"
#include "vorticity.h"

/* main function */
int main (int argc, const char * argv[]) {
    double t1, t2, runtime, nextoutputtime, firstoutputtime;
    int    df, ss, ii, jj, ee, i, j, k, startindex = 0, num_integrated, outputframe = 0, outputcount = 0;
    char   BinFile[LONGSTRING];
    time_t SimulationStartTime, InitializationTime, SimulationStopTime;
    FILE   *Trace_BinFileID;
    FileData *Trace_BinFileArray;
    LagrangianPoint TempPoint;
    ReleaseLocation *rnode;
    
    /* Get simulation start time so we can time how long simulation runs */
    time(&SimulationStartTime);
    
    /* Input paramenters */
    ReadInParameters(argc, argv);
    CheckParameters();
    SetDerivedParameters();
    
    /* Load mesh */
    LoadMeshData();
    
    /* If performing a staggered tracer release, compute release times at each release location based on velocity data */
    if(Trace_Compute && Trace_ReleaseStrategy == STAGGERED)
        GenerateStaggeredRelease();
    
    /* Load velocity for first frame needed to start integration */
    if(Data_MeshType == CARTESIAN)
        LoadCartVelDataFrame(Data_FirstFrame);
    else if(Data_MeshType == UNSTRUCTURED) {
        if(Int_NormalFlow)
            LoadMeshSurfaceNormals(); /* Mesh normals needed to impose inflow along mesh walls */
        LoadUnstructVelDataFrame(Data_FirstFrame);
    }
    
    /* Do initialization */
    if(VelOut_Compute) {
        if(VelOut_GenerateMesh)
            GenerateVelOutMesh();
        else
            ReadInVelOutMesh();
    }
    if(FTLE_Compute) {
        InitializeFTLEArray();
        Trace_BinFileID = NULL;
        Trace_BinFileArray = NULL;
    }
    else if(Trace_Compute) {
        /* If computing exposure time data, allocate memory and initialize exposure time array */
        if(Trace_CETCompute) {
            if((Trace_CETArray = (CETE *)malloc(CET_MeshNumElements * sizeof(CETE))) == NULL)
                FatalError("Malloc failed for Trace_CETArray");
            for(ee = 0; ee < CET_MeshNumElements; ee++) {
                Trace_CETArray[ee].CETsum = 0; /* Cummulative exposure time for element ee */
                Trace_CETArray[ee].Encounters = 0; /* Numer of encounters for element ee */
            }
        }
        if(Trace_ReleaseStrategy == STAGGERED) {
            Trace_BinFileID = NULL;
            /* Generate an output file for each output time */
            if((Trace_BinFileArray = (FileData *)malloc(Output_TRes * sizeof(FileData))) == NULL)
                FatalError("Malloc failed for Trace_BinFileArray");
            for(ss = 0; ss < Output_TRes; ss++) {
                /* Open file */
                sprintf(Trace_BinFileArray[ss].FilePath, "%s%s.%d.bin", Path_Output, Trace_OutFilePrefix, ss);
                if((Trace_BinFileArray[ss].FileID = fopen(Trace_BinFileArray[ss].FilePath, "wb")) == NULL)
                    FatalError("Could not open file %s", Trace_BinFileArray[ss].FilePath);
                /* Write output time as first entry in file */
                nextoutputtime = Output_TStart + ss * Output_TDelta;
                if(fwrite(&nextoutputtime, sizeof(double), 1, Trace_BinFileArray[ss].FileID) < 1)
                    FatalError("Could not write output time to %s", Trace_BinFileArray[ss].FilePath);
            }
        }
        else {
            Trace_BinFileArray = NULL;
            /* Set up mesh of tracers to be integrated */
            if(!Trace_GenerateMesh)
                ReadInTraceMesh();
            else
                GenerateTracerMesh();
            /* Open binary file for tracer output */
            if((Trace_BinFileID = fopen(TraceOUT_BinFilePath, "wb")) == NULL)
                FatalError("Could not open file %s", TraceOUT_BinFilePath);
            /* Initialize output counters */
            outputframe = 0;
            nextoutputtime = Output_TStart;
            for(ss = 0; ss < Trace_NumLaunchTimes; ss++)
                CreateNewTraceLaunch(ss);
            /* If computing activation potential (integrated strain rate), then load first strain rate field */
            if(Trace_APCompute)
                LoadStrainRateDataFrame(Data_FirstFrame);
        }
        if(Int_NormalFlow && Particle_Radius > TINY)
            LoadBoundaryElementFlags();
        
        /* Allocate memory for array to store number of tracers output each output time frame */
        if((Trace_NumOutput = (int *)calloc(Output_TRes, sizeof(int))) == NULL)
            FatalError("Calloc failed for Trace_NumOutput");
        
        /* If outputting vorticity of each particle, then load first vorticity field */
        if(Trace_VorticityCompute)
            LoadVorticityDataFrame(Data_FirstFrame);
    } /* end if(Trace_Compute) */
    else {
        Trace_BinFileID = NULL;
        Trace_BinFileArray = NULL;
    }
    
    /* Output initialization time */
    time(&InitializationTime);
    runtime = difftime(InitializationTime, SimulationStartTime);
    printf("\nInitialization time: %d hr %d min %d s\n", (int)(runtime/3600.0),
           (int)(fmod(runtime, 3600.0)/60.0), (int)(fmod(fmod(runtime, 3600.0), 60.0)));
    fflush(stdout);
    
    /*** MAIN LOOP ***/
    for(df = Data_FirstFrame; df != Data_LastFrame; df = df + Int_TimeDirection) {
        
        /* Load next velocity data frame */
        if(Data_MeshType == CARTESIAN)
            LoadCartVelDataFrame(df + Int_TimeDirection);
        else if(Data_MeshType == UNSTRUCTURED)
            LoadUnstructVelDataFrame(df + Int_TimeDirection);
        
        /* If computing activation potential, load strain rate data */
        if(Trace_Compute && (Trace_ReleaseStrategy != STAGGERED) && Trace_APCompute)
            LoadStrainRateDataFrame(df + Int_TimeDirection);
        /* If outputting vorticity, load vorticity data */
        if(Trace_Compute && Trace_VorticityCompute)
            LoadVorticityDataFrame(df + Int_TimeDirection);
        
        /* Determine time interval of loaded data */
        if(Int_TimeDirection > 0) {
            Data_LoadedTMin = Data_TMin + df * Data_TDelta;
            Data_LoadedTMax = Data_LoadedTMin + Data_TDelta;
        }
        else {
            Data_LoadedTMin = Data_TMin + (df - 1) * Data_TDelta;
            Data_LoadedTMax = Data_LoadedTMin + Data_TDelta;
        }
        printf("Data loaded in memory spans t = %g to %g\n", Data_LoadedTMin, Data_LoadedTMax);
        fflush(stdout);
        
        /**************************************** Perform Computations *******************************************/
        /* If requested, interpolate velocity data to output mesh */
        if(VelOut_Compute)
            for(ss = 0; ss < Output_TRes; ss++)
            /* Check if current output time for ss is in time interval of loaded data and that output hasn't already been completed */
                if((VelOut_OutputTime[ss] - Data_LoadedTMin > -TINY) && (VelOut_OutputTime[ss] - Data_LoadedTMax < TINY) && !VelOut_Complete[ss]) {
                    OutputVelOut(VelOut_OutputTime[ss], ss);
                    /* OutputStrainOut(VelOut_OutputTime[ss], ss); */
                    VelOut_Complete[ss] = 1;
                }
        
        /* If requested, compute FTLE */
        if(FTLE_Compute) {
            for(ss = 0; ss < Output_TRes; ss++) { /* For output time */
                /* Consider only output times that are not complete AND, have been launched OR need to be launched */
                if(FTLE_Launches[ss].Status != COMPLETE &&
                   (FTLE_Launches[ss].Status == LAUNCHED ||
                    (FTLE_Launches[ss].StartTime <= Data_LoadedTMax + TINY && FTLE_Launches[ss].StartTime >= Data_LoadedTMin - TINY))) {
                       ReadInFTLELaunch(ss); /* Load data for release ss from bin file to FTLE_MeshPt */
                       /* Initialize slide if it hasn't been launched yet and set integration bounds */
                       if(FTLE_Launches[ss].Status == UNLAUNCHED) {
                           FTLE_Launches[ss].Status = LAUNCHED;
                           t1 = FTLE_Launches[ss].StartTime;
                           /* Loop over all the points */
                           for(k = 0; k < FTLE_CartMesh.ZRes; k++)
                               for(j = 0; j < FTLE_CartMesh.YRes; j++)
                                   for(i = 0; i < FTLE_CartMesh.XRes; i++)
                                       if(!FTLE_MeshPt[i][j][k].Pt.LeftDomain) {
                                           if(Particle_Radius > TINY && Particle_ICType == 0) {
                                               FTLE_MeshPt[i][j][k].Pt.V[0] = 0.0;
                                               FTLE_MeshPt[i][j][k].Pt.V[1] = 0.0;
                                               FTLE_MeshPt[i][j][k].Pt.V[2] = 0.0;
                                           }
                                           else {
                                               GetVelocity(FTLE_Launches[ss].StartTime, &(FTLE_MeshPt[i][j][k].Pt), FTLE_MeshPt[i][j][k].Pt.V);
                                           }
                                       }
                           printf("Initiated FTLE output %d at time %g\n", ss, t1);
                           fflush(stdout);
                       }
                       else
                           t1 = (Int_TimeDirection > 0) ? Data_LoadedTMin : Data_LoadedTMax;
                       t2 = (Int_TimeDirection > 0) ?
                       fmin(Data_LoadedTMax, FTLE_Launches[ss].StopTime) : fmax(Data_LoadedTMin, FTLE_Launches[ss].StopTime);
                       
                       /* Integrate release ss from time t1 to time t2 */
                       printf("Integrating FTLE grid initiated at %g from %g to %g...\n", FTLE_Launches[ss].StartTime, t1, t2);
                       fflush(stdout);
                       num_integrated = 0;
                       /* Loop over all the points */
                       for(k = 0; k < FTLE_CartMesh.ZRes; k++) {
                           for(j = 0; j < FTLE_CartMesh.YRes; j++) {
                               for(i = 0; i < FTLE_CartMesh.XRes; i++) {
                                   /* Make sure point initially inside of domain AND it hasn't since left, unless extrapolating position */
                                   if((FTLE_MeshPt[i][j][k].Pt.LeftDomainTime > -TINY) && (!FTLE_MeshPt[i][j][k].Pt.LeftDomain || Int_Extrapolate)) {
                                       /* Create TempPoint to integrate in case point leaves the domain */
                                       TempPoint.X[0] = FTLE_MeshPt[i][j][k].Pt.X[0];
                                       TempPoint.X[1] = FTLE_MeshPt[i][j][k].Pt.X[1];
                                       TempPoint.X[2] = FTLE_MeshPt[i][j][k].Pt.X[2];
                                       TempPoint.LeftDomain = FTLE_MeshPt[i][j][k].Pt.LeftDomain;
                                       TempPoint.LeftDomainTime = FTLE_MeshPt[i][j][k].Pt.LeftDomainTime;
                                       TempPoint.V[0] = FTLE_MeshPt[i][j][k].Pt.V[0];
                                       TempPoint.V[1] = FTLE_MeshPt[i][j][k].Pt.V[1];
                                       TempPoint.V[2] = FTLE_MeshPt[i][j][k].Pt.V[2];
                                       if(Data_MeshType == UNSTRUCTURED)
                                           TempPoint.ElementIndex = FTLE_MeshPt[i][j][k].Pt.ElementIndex;
                                       
                                       Advect(&TempPoint, t1, t2);
                                       num_integrated++;
                                       if(TempPoint.LeftDomain && !Int_Extrapolate) {
                                           /* Compute FTLE at point (early) */
                                           if(!FTLE_MeshPt[i][j][k].HaveFTLE)
                                               GetFTLEForPointEarly(ss, i, j, k, t1, TempPoint.LeftDomainTime);
                                           /* Compute FTLE for neighbors (early) */
                                           if(i > 0)  /* Get FTLE for i-1 neighbor */
                                               if(!FTLE_MeshPt[i-1][j][k].HaveFTLE)
                                                   GetFTLEForPointEarly(ss, i-1, j, k, t1, TempPoint.LeftDomainTime);
                                           if(i < FTLE_CartMesh.XRes - 1)  /* Get FTLE for i+1 neighbor */
                                               if(!FTLE_MeshPt[i+1][j][k].HaveFTLE)
                                                   GetFTLEForPointEarly(ss, i+1, j, k, t1, TempPoint.LeftDomainTime);
                                           if(j > 0)  /* Get FTLE for j-1 neighbor */
                                               if(!FTLE_MeshPt[i][j-1][k].HaveFTLE)
                                                   GetFTLEForPointEarly(ss, i, j-1, k, t1, TempPoint.LeftDomainTime);
                                           if(j < FTLE_CartMesh.YRes - 1) /* Get FTLE for j+1 neighbor */
                                               if(!FTLE_MeshPt[i][j+1][k].HaveFTLE)
                                                   GetFTLEForPointEarly(ss, i, j+1, k, t1, TempPoint.LeftDomainTime);
                                           if(k > 0)  /* Get FTLE for k-1 neighbor */
                                               if(!FTLE_MeshPt[i][j][k-1].HaveFTLE)
                                                   GetFTLEForPointEarly(ss, i, j, k-1, t1, TempPoint.LeftDomainTime);
                                           if(k < FTLE_CartMesh.ZRes - 1)  /* Get FTLE for k+1 neighbor */
                                               if(!FTLE_MeshPt[i][j][k+1].HaveFTLE)
                                                   GetFTLEForPointEarly(ss, i, j, k+1, t1, TempPoint.LeftDomainTime);
                                           
                                           /* Set LeftDomain flag last, as GetFTLEForPointEarly calls require integration of FTLE_MeshPt[i][j][k].Pt, which will fail if this flag is true */
                                           FTLE_MeshPt[i][j][k].Pt.LeftDomain = 1;
                                       }
                                       else {
                                           /* Point still in domain, save position for update of FTLE_MeshPt once all points integrated over current time interval */
                                           FTLE_NewArray[i][j][k][0] = TempPoint.X[0];
                                           FTLE_NewArray[i][j][k][1] = TempPoint.X[1];
                                           FTLE_NewArray[i][j][k][2] = TempPoint.X[2];
                                           FTLE_MeshPt[i][j][k].Pt.V[0] = TempPoint.V[0];
                                           FTLE_MeshPt[i][j][k].Pt.V[1] = TempPoint.V[1];
                                           FTLE_MeshPt[i][j][k].Pt.V[2] = TempPoint.V[2];
                                           if(Int_Extrapolate)
                                               FTLE_MeshPt[i][j][k].Pt.LeftDomain = TempPoint.LeftDomain;
                                           if(Data_MeshType == UNSTRUCTURED)
                                               FTLE_MeshPt[i][j][k].Pt.ElementIndex = TempPoint.ElementIndex;
                                       }
                                   }
                               }
                           }
                       }
                       UpdateFTLELocations(); /* Update position of points still in domain, i.e. copy FTLE_NewArray values to FTLE_MeshPt */
                       WriteOutFTLELaunch(ss); /* Store FTLE_MeshPt in a bin file for slide ss */
                       printf(" %d points integrated.\n", num_integrated);
                       fflush(stdout);
                       
                       /* Write FTLE to output file if slide ss is complete  */
                       if(((Int_TimeDirection > 0) ? (FTLE_Launches[ss].StopTime - Data_LoadedTMax) : (Data_LoadedTMin - FTLE_Launches[ss].StopTime)) < TINY) {
                           printf("Stopping FTLE output slide %d and writing to output file.\n", ss);
                           fflush(stdout);
                           OutputFTLE(ss, FTLE_IntTLength);
                           FTLE_Launches[ss].Status = COMPLETE;
                       }
                       /* Or write FTLE to output file if computing variation of FTLE with integration time */
                       else if(FTLE_ComputeVariation && !fmod(df - Data_FirstFrame, FTLE_VariationOutFreq)) {
                           printf("Exporting FTLE values at integration length %f.\n", fabs(t2 - FTLE_Launches[ss].StartTime)); fflush(stdout);
                           OutputFTLE(ss, fabs(t2 - FTLE_Launches[ss].StartTime));
                       }
                   }
            }
        }
        else if(Trace_Compute) {
            if(Trace_ReleaseStrategy == STAGGERED) {
                num_integrated = 0;
                /* Loop over all release locations */
                for(ss = 0; ss < Trace_NumReleaseLocations; ss++) {
                    /* Loop over all release times */
                    for(rnode = Trace_ReleaseList[ss]; rnode != NULL; rnode = rnode->next) {
                        /* Make sure point launched, or needs to be launched and not complete */
                        if(rnode->slide.Status != COMPLETE &&
                           (rnode->slide.Status == LAUNCHED || (rnode->slide.StartTime >= Data_LoadedTMin && rnode->slide.StartTime <= Data_LoadedTMax))) {
                            /* Make sure point still in domain  */
                            if(!rnode->pt.LeftDomain || Int_Extrapolate) {
                                /* Check if point hasn't been launched yet */
                                if(rnode->slide.Status == UNLAUNCHED) {
                                    rnode->slide.Status = LAUNCHED;
                                    /* Increment number of encounters in element where point is launched */
                                    if(Trace_CETCompute)
                                        Trace_CETArray[Trace_CETAuxillaryMesh? rnode->pt.AuxElementIndex : rnode->pt.ElementIndex].Encounters++;
                                    /* If launchtime is Output_TStart, then write out locations */
                                    if(fabs(rnode->slide.StartTime - Output_TStart) < TINY) {
                                        if(fwrite(rnode->pt.X, sizeof(double), 3, Trace_BinFileArray[0].FileID) < 3)
                                            FatalError("Could not write rnode->pt.X to %s", Trace_BinFileArray[0].FilePath);
                                        if(Trace_VorticityCompute)
                                            if(fwrite(rnode->pt.vorticity, sizeof(double), 3, Trace_BinFileArray[0].FileID) < 3)
                                                FatalError("Could not write rnode->pt.vorticity to %s", Trace_BinFileArray[0].FilePath);
                                        /*if(Trace_APCompute)
                                         if(fwrite(&(rnode->pt.Scalar), sizeof(double), 1, Trace_BinFileArray[0].FileID) < 1)
                                         FatalError("Could not write rnode->pt.Scalar to %s", Trace_BinFileArray[0].FilePath);*/
                                        Trace_NumOutput[0]++;
                                    }
                                }
                                /* Integrate point forward */
                                t1 = fmax(rnode->slide.StartTime, Data_LoadedTMin);
                                t2 = fmin(rnode->slide.StopTime, Data_LoadedTMax);
                                /* See how many times we need to output data over time interval from t1 to t2 */
                                outputcount = 0;
                                for(ii = 0; ii < Output_TRes; ii++)
                                    if((Output_TStart + ii*Output_TDelta > t1 + TINY) && (Output_TStart + ii*Output_TDelta <= t2 + TINY))
                                        outputcount++;
                                /* Break up integration if needed */
                                if(outputcount == 0)  /* not needed */
                                    Advect(&rnode->pt, t1, t2);
                                else { /* break up integration */
                                    /* Determine next output time */
                                    for(ii = 0; ii < Output_TRes; ii++) {
                                        nextoutputtime = Output_TStart + ii*Output_TDelta;
                                        if(nextoutputtime > t1 + TINY) {
                                            startindex = ii;
                                            break;
                                        }
                                    }
                                    t2 = nextoutputtime;
                                    Advect(&rnode->pt, t1, t2);
                                    /* Output to file */
                                    if(!rnode->pt.LeftDomain || Int_Extrapolate || Trace_AlwaysOutput) {
                                        if(fwrite(rnode->pt.X, sizeof(double), 3, Trace_BinFileArray[startindex].FileID) < 3)
                                            FatalError("Could not write rnode->pt.X to %s", Trace_BinFileArray[startindex].FilePath);
                                        if(Trace_VorticityCompute)
                                            if(fwrite(&(rnode->pt.vorticity), sizeof(double), 3, Trace_BinFileArray[startindex].FileID) < 3)
                                                FatalError("Could not write rnode->pt.vorticity to %s", Trace_BinFileArray[startindex].FilePath);
                                        /*if(Trace_APCompute)
                                         if(fwrite(&(rnode->pt.Scalar), sizeof(double), 1, Trace_BinFileArray[ii].FileID) < 1)
                                         FatalError("Could not write rnode->pt.Scalar to %s", Trace_BinFileArray[ii].FilePath);*/
                                        Trace_NumOutput[startindex]++;
                                    }
                                    /* Break up integration for further output times (if needed) */
                                    for(ii = 1; ii < outputcount; ii++) {
                                        if(!rnode->pt.LeftDomain || Int_Extrapolate) {
                                            t1 = t2;
                                            t2 = t1 + Output_TDelta;
                                            Advect(&rnode->pt, t1, t2);
                                            /* Output to file */
                                            if(!rnode->pt.LeftDomain || Int_Extrapolate  || Trace_AlwaysOutput) {
                                                if(fwrite(rnode->pt.X, sizeof(double), 3, Trace_BinFileArray[startindex + ii].FileID) < 3)
                                                    FatalError("Could not write rnode->pt.X to %s", Trace_BinFileArray[startindex + ii].FilePath);
                                                if(Trace_VorticityCompute)
                                                    if(fwrite(rnode->pt.vorticity, sizeof(double), 3, Trace_BinFileArray[startindex + ii].FileID) < 3)
                                                        FatalError("Could not write rnode->pt.vorticity to %s", Trace_BinFileArray[startindex + ii].FilePath);
                                                /*if(Trace_APCompute)
                                                 if(fwrite(&(rnode->pt.Scalar), sizeof(double), 1, Trace_BinFileArray[(int)((t2 + TINY)/Output_TDelta)].FileID) < 1)
                                                 FatalError("Could not write rnode->pt.Scalar to %s", Trace_BinFileArray[(int)((t2 + TINY)/Output_TDelta)].FilePath);*/
                                                Trace_NumOutput[startindex + ii]++;
                                            }
                                        }
                                        else if(Trace_AlwaysOutput) {
                                            while(ii < outputcount) {
                                                if(fwrite(rnode->pt.X, sizeof(double), 3, Trace_BinFileArray[startindex + ii].FileID) < 3)
                                                    FatalError("Could not write rnode->pt.X to %s", Trace_BinFileArray[startindex + ii].FilePath);
                                                if(Trace_VorticityCompute)
                                                    if(fwrite(&(rnode->pt.vorticity), sizeof(double), 3, Trace_BinFileArray[startindex+ii].FileID) < 3)
                                                        FatalError("Could not write rnode->pt.vorticity to %s", Trace_BinFileArray[startindex + ii].FilePath);
                                                Trace_NumOutput[startindex + ii]++;
                                            }
                                            break;
                                        }
                                    }
                                    /* Continue until end of interval or StopTime */
                                    if(!rnode->pt.LeftDomain || Int_Extrapolate) {
                                        t1 = t2;
                                        t2 = fmin(rnode->slide.StopTime, Data_LoadedTMax);
                                        if(t2 - t1 > 0) {
                                            Advect(&rnode->pt, t1, t2);
                                        }
                                        else if(t2 - t1 < -TINY)
                                            FatalError("Integration of point beyond requested stop time or range of loaded data");
                                    }
                                }
                                if(rnode->slide.StopTime <= Data_LoadedTMax)
                                    rnode->slide.Status = COMPLETE;
                                num_integrated++;
                            } /* end if particle inside domain */
                            else if(Trace_AlwaysOutput) {
                                t1 = fmax(rnode->slide.StartTime, Data_LoadedTMin);
                                t2 = fmin(rnode->slide.StopTime, Data_LoadedTMax);
                                /* Determine number of outputs */
                                outputcount = 0;
                                for(ii = 0; ii < Output_TRes; ii++)
                                    if((Output_TStart + ii*Output_TDelta > t1 + TINY) && (Output_TStart + ii*Output_TDelta <= t2 + TINY))
                                        outputcount++;
                                /* Determine first output time */
                                for(ii = 0; ii < Output_TRes; ii++) {
                                    firstoutputtime = Output_TStart + ii*Output_TDelta;
                                    if(firstoutputtime > t1 + TINY)
                                        break;
                                }
                                /* Write data to each output file */
                                for(jj=0; jj<outputcount; jj++) {
                                    if(fwrite(rnode->pt.X, sizeof(double), 3, Trace_BinFileArray[ii+jj].FileID) < 3)
                                        FatalError("Could not write rnode->pt.X to %s", Trace_BinFileArray[ii+jj].FilePath);
                                    if(Trace_VorticityCompute)
                                        if(fwrite(&(rnode->pt.vorticity), sizeof(double), 3, Trace_BinFileArray[ii+jj].FileID) < 3)
                                            FatalError("Could not write rnode->pt.vorticity to %s", Trace_BinFileArray[ii+jj].FilePath);
                                    Trace_NumOutput[ii+jj]++;
                                }
                            }
                        } /* If particle should be advected during this time interval */
                    } /* For each particle released from this location */
                } /* For each release location */
                printf(" %d tracers integrated.\n", num_integrated); fflush(stdout);
            } /* Staggered release */
            else { /* Normal release */
                
                /* Set start of integration interval */
                t1 = (Int_TimeDirection > 0) ? Data_LoadedTMin : Data_LoadedTMax;
                
                /* Integrate until end of loaded reached */
                while(((Int_TimeDirection > 0) ? (t1 - Data_LoadedTMax) : (Data_LoadedTMin - t1)) < 0) {
                    
                    /* Set end of integration interval to next output time, or end of loaded data */
                    t2 = (Int_TimeDirection > 0) ? fmin(nextoutputtime, Data_LoadedTMax) : fmax(Data_LoadedTMin, nextoutputtime);
                    
                    /* Loop over each release */
                    for(ss = 0; ss < Trace_NumLaunchTimes; ss++) {
                        
                        /* Consider only releases that have been launched, or need to be launched */
                        if(Trace_Launches[ss].Status != COMPLETE &&
                           (Trace_Launches[ss].Status == LAUNCHED || (Trace_Launches[ss].StartTime <= fmax(t1, t2) && Trace_Launches[ss].StartTime >= fmin(t1, t2)))) {
                            
                            /* Load tracer data from file to Trace_MeshPt */
                            ReadInTraceLaunch(ss);
                            
                            /* Initiate release if needed */
                            if(Trace_Launches[ss].Status == UNLAUNCHED) {
                                Trace_Launches[ss].Status = LAUNCHED;
                                /* Update encounter counter for elements containing each release point */
                                if(Trace_CETCompute)
                                    for(ii = 0; ii < Trace_NumTracers; ii++)
                                        if(!Trace_MeshPt[ii].LeftDomain)
                                            Trace_CETArray[Trace_CETAuxillaryMesh? Trace_MeshPt[ii].AuxElementIndex : Trace_MeshPt[ii].ElementIndex].Encounters++;
                                /* Reset start of integration to release time */
                                t1 = Trace_Launches[ss].StartTime;
                                /* Initialize velocity of each particle */
                                for(ii = 0; ii < Trace_NumTracers; ii++) {
                                    if(!Trace_MeshPt[ii].LeftDomain) {
                                        if(Particle_Radius > TINY && Particle_ICType == 0) {
                                            Trace_MeshPt[ii].V[0] = 0.0;
                                            Trace_MeshPt[ii].V[1] = 0.0;
                                            Trace_MeshPt[ii].V[2] = 0.0;
                                        }
                                        else
                                            GetVelocity(Trace_Launches[ss].StartTime, &Trace_MeshPt[ii], Trace_MeshPt[ii].V);
                                    }
                                }
                                printf("Initiated release %d for %d tracers with a launch time at %g\n", ss, Trace_NumTracers, t1); fflush(stdout);
                            }
                            
                            printf("Integrating release %d from %g to %g...", ss, t1, t2); fflush(stdout);
                            num_integrated = 0;
                            /* Loop over tracers */
                            for(ii = 0; ii < Trace_NumTracers; ii++) {
                                if(Trace_MeshPt[ii].LeftDomainTime > -TINY && (!Trace_MeshPt[ii].LeftDomain || Int_Extrapolate)) {
                                    Advect(&Trace_MeshPt[ii], t1, t2);
                                    num_integrated++;
                                }
                                /* If nextoutputtime reached, output tracer positions to file */
                                if(fabs(nextoutputtime - t2) < TINY) {
                                    if(Trace_MeshPt[ii].LeftDomainTime > -TINY && (!Trace_MeshPt[ii].LeftDomain || Trace_AlwaysOutput)) {
                                        if(fwrite(Trace_MeshPt[ii].X, sizeof(double), 3, Trace_BinFileID) < 3)
                                            FatalError("Could not write Trace_MeshPt[%d].X to %s", ii, TraceOUT_BinFilePath);
                                        if(Trace_VorticityCompute)
                                            if(fwrite(Trace_MeshPt[ii].vorticity, sizeof(double), 3, Trace_BinFileID) < 3)
                                                FatalError("Could not write Trace_MeshPt[%d].vorticity to %s", ii, TraceOUT_BinFilePath);
                                        /*if(Trace_APCompute)
                                         if(fwrite(&Trace_MeshPt[ii].Scalar, sizeof(double), 1, Trace_BinFileID) < 1)
                                         FatalError("Could not write Trace_MeshPt[%d].Scalar to %s", ii, TraceOUT_BinFilePath);*/
                                        Trace_NumOutput[outputframe]++;
                                    }
                                }
                            }
                            printf(" %d tracers integrated.\n", num_integrated);
                            fflush(stdout);
                            WriteOutTraceLaunch(ss); /* Saves Trace_MeshPt array to binary file for current slide */
                        }
                        if(((Int_TimeDirection > 0) ? (Trace_Launches[ss].StopTime - t2) : (t2 - Trace_Launches[ss].StopTime)) < TINY)
                            Trace_Launches[ss].Status = COMPLETE;
                    } /* for(ss = 0; ss < Trace_NumLaunchTimes; ss++) */
                    if(fabs(t2 - nextoutputtime) < TINY) {
                        printf("Number of tracers output: %d\n", Trace_NumOutput[outputframe]);
                        fflush(stdout);
                        nextoutputtime += Int_TimeDirection * Output_TDelta;
                        outputframe++;
                        if(outputframe >= Output_TRes)
                            break;
                    }
                    t1 = t2;
                }
            }
        }
    }
    
    /* Do any additional clean up */
    if(FTLE_Compute) {
        if(FTLE_ComputeVariation) {
            /* Delete temp FTLE files */
            sprintf(BinFile, "%s%s.0.bin", Path_Data, FTLE_OutFilePrefix);
            if(remove(BinFile))
                fprintf(stderr, "Warning: Could not delete file %s\n", BinFile);
        }
    }
    if(Trace_Compute) {
        /* Output RT results */
        if(Trace_RTCompute)
            OutputRT();
        else if(Trace_CETCompute)
            OutputCET();
        
        /* Output AP field */
        if(Trace_APCompute)
            OutputAP();
        
        /* Output tracer positions */
        if(Trace_ReleaseStrategy == STAGGERED) {
            for(ss = 0; ss < Trace_NumLaunchTimes; ss++)
                if(fclose(Trace_BinFileArray[ss].FileID))
                    fprintf(stderr, "Warning: could not close file %s\n", Trace_BinFileArray[ss].FilePath);
            free(Trace_BinFileArray);
        }
        else {
            /* Remove temporary bin files */
            for(ss = 0; ss < Trace_NumLaunchTimes; ss++) {
                sprintf(BinFile, "%s.%d.bin", TraceN_BinFilePathPrefix, ss);
                if(remove(BinFile))
                    fprintf(stderr, "Warning: Could not delete file %s\n", BinFile);
            }
            fclose(Trace_BinFileID);
            /* Convert tracer position data from temp files to output files */
            OutputTracers();
        }
    }
    
    printf("\nCleaning memory\n");
    fflush(stdout);
    
    FreeVelFieldData();
    
    if(FTLE_Compute) {
        FreeFTLEData();
    }
    if(Trace_Compute) {
        FreeTracerData();
        if(Trace_VorticityCompute)
            FreeVorticityData();
        if(Trace_APCompute)
            FreeStrainRateData();
    }
    if(VelOut_Compute) {
        FreeVelOutData();
    }
    
    printf("\nSimulation complete!\n");
    
    time(&SimulationStopTime);
    runtime = difftime(SimulationStopTime, SimulationStartTime);
    printf("Total run time: %d:%d:%d\n", (int)(runtime/3600), (int)(fmod(runtime, 3600)/60), (int)(fmod(fmod(runtime, 3600), 60)));
    printf("Processor time: %d:%d:%d\n", (int)(clock()/CLOCKS_PER_SEC/3600), (int)(fmod(clock()/CLOCKS_PER_SEC, 3600)/60),
           (int)(fmod(fmod(clock()/CLOCKS_PER_SEC, 3600), 60)));
    
    
    return(0);
}
