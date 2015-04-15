/*
 *  integration.c
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
#include "integration.h"
#include "io.h"
#include "macros.h"
#include "mesh.h"
#include "mymath.h"
#include "strainrate.h"
#include "vorticity.h"
#include "velocity.h"

void Advect(LagrangianPoint *pt, double tstart, double tend) {
    
    int ii;
    double stoptime;
    
    if(pt->LeftDomain && Int_Extrapolate)
        for(ii = 0; ii < 3; ii++)
            pt->X[ii] = pt->X[ii] + pt->V[ii] * (tend - tstart);
    else {
        if(Particle_Radius > TINY)
            stoptime = pEuler(pt, tstart, tend);
        else {
            if(Int_Type == 1)
                stoptime = RK4(pt, tstart, tend);
            else if(Int_Type == 2)
                stoptime = RKF(pt, tstart, tend);
            else
                stoptime = Euler(pt, tstart, tend);
        }
        if((fabs(stoptime - tend) > TINY) && Int_Extrapolate)
            for(ii = 0; ii < 3; ii++)
                pt->X[ii] = pt->X[ii] + pt->V[ii] * (tend - stoptime);
    }
}

double pEuler(LagrangianPoint *pt, double tstart, double tend) {
    /***
     psuedo-Euler integration routine of x' = v and v' = f
     where f defined by Maxey-Riley ignoring Fauxen correction and memory terms.
     Note, v'= f solved analytically and then x' = v updated via Euler.
     ***/
    
    double h, tc, dx, X0[3], Xs[3], vec1[3], vec2[3], u0[3], dudt0[3], dotu0[3], gradu0[3][3], wp[3], strain_rate0 = 0, strain_rate1 = 0;
    double u, v, w, d, bx1[3], bx2[3], abn[3], vt[3], vn[3];
    int i, ii, jj, sindex, lbdnodes[3], bdnodes[3], count;
    LagrangianPoint MPA[6];
    
    if(pt->LeftDomain)
        FatalError("Attempting to integrate particle that already left domain");
    if(Data_MeshType == UNSTRUCTURED)
        if(pt->ElementIndex == -1)
            FatalError("Attempting to integrate particle with Element_Index = -1");
    
    /* Initialize local time variables */
    tc = tstart;
    h = Int_TimeStep;
    
    /* Choose (arbitrary) spacing for finite differencing */
    dx = 0.001;
    
    
    /* Intergrate until tc reaches tend */
    while((tc + TINY) < tend) { /* Integration not complete */
        /* If trying to integrate past tend, reset step size, h, to integrate until tend */
        if ((tc + h - tend) > 0.0)
            h = tend - tc;
        
        /* Compute strain rate */
        if(Trace_Compute && (Trace_ReleaseStrategy != STAGGERED) && Trace_APCompute) {
            if(Data_MeshType == CARTESIAN)
                strain_rate0 = GetStrainRate_Cartesian(tc, pt);
            else if(Data_MeshType == UNSTRUCTURED)
                strain_rate0 = GetStrainRate_Unstructured(tc, pt);
        }
        
        /* fluid velocity at time tc and particle position x(tc) */
        GetVelocity(tc, pt, u0);
        
        /*** Collision Model ***/
        /* Check if using collisions */
        if(Int_NormalFlow) {
            /* Loop over sphere points */
            for(ii = 0; ii < 26; ii++) {
                /* Locate element containing sphere point Xs[ii] */
                Xs[0] = pt->X[0] + Particle_Radius * UnitSphere[ii][0];
                Xs[1] = pt->X[1] + Particle_Radius * UnitSphere[ii][1];
                Xs[2] = pt->X[2] + Particle_Radius * UnitSphere[ii][2];
                sindex = Get_Element_Local_Search(Xs, pt->ElementIndex);
                if(sindex >= 0) {
                    /* Check if Xs[ii] in boundary element */
                    if(Vel_MeshElementFlagArray[sindex]) {
                        count = 0;
                        /* Check how many nodes on boundary */
                        if(Vel_MeshElementFlagArray[sindex] == 1) { /* Only one node on boundary */
                            /* Determine which nodes have zero velocity */
                            for(i=0; i<(Dimensions+1); i++) { /* loop over nodes */
                                u = Vel_UnstructVelArray_U[0][Vel_MeshElementArray[sindex].Nodes[i]];
                                v = Vel_UnstructVelArray_V[0][Vel_MeshElementArray[sindex].Nodes[i]];
                                w = Vel_UnstructVelArray_W[0][Vel_MeshElementArray[sindex].Nodes[i]];
                                if(u*u + v*v + w*w < TINY) {
                                    /* Save off index for the no slip node */
                                    bdnodes[0] = Vel_MeshElementArray[sindex].Nodes[i];
                                    count++;
                                    if(count > 1) FatalError("Element %d has %d no slip nodes but expecting 1.", sindex, count);
                                }
                            }
                            if(count < 1) FatalError("Not enough no slip nodes");
                            /* Get local index of surface normal vectors */
                            lbdnodes[0] = GetLocalIndex_SurfaceMesh(bdnodes[0]);
                            /* Get unit normal */
                            for(i=0; i<3; i++)
                                abn[i] = Vel_SurfaceMeshInwardNormals[lbdnodes[0]][i];
                            /* Get coordinates of boundary node */
                            for(i=0; i<3; i++)
                                bx1[i] = Vel_MeshNodeArray[bdnodes[0]][i];
                            /* Get distance from sphere point to boundary point */
                            d = dist(bx1, Xs, Dimensions);
                            /* Check for acutal collision (tolerance set to 15% of particle radius) */
                            if(d <= 0.15*Particle_Radius) {
                                /* Compute normal component and reverse normal component if needed */
                                if(vdot(abn, pt->V, Dimensions) < 0) {
                                    for(i=0; i<3; i++) {
                                        vn[i] = vdot(abn, pt->V, Dimensions)*abn[i];
                                        vt[i] = pt->V[i] - vn[i];
                                    }
                                    for(i=0; i<3; i++)
                                        pt->V[i] = vt[i] - Int_NormalFlowScaling*vn[i]; /* Int_NormalFlowScaling used as restitution coefficient */
                                }
                                break; /* If collision detected, no need to check other points */
                            }
                        }
                        else if(Vel_MeshElementFlagArray[sindex] == 2) { /* Two boundary nodes */
                            /* Determine which nodes have zero velocity */
                            for(i=0; i<(Dimensions+1); i++) {
                                u = Vel_UnstructVelArray_U[0][Vel_MeshElementArray[sindex].Nodes[i]];
                                v = Vel_UnstructVelArray_V[0][Vel_MeshElementArray[sindex].Nodes[i]];
                                w = Vel_UnstructVelArray_W[0][Vel_MeshElementArray[sindex].Nodes[i]];
                                if(u*u + v*v + w*w < TINY) {
                                    /* Save off index for the no slip node */
                                    bdnodes[count] = Vel_MeshElementArray[sindex].Nodes[i];
                                    count++;
                                    if(count > 2) FatalError("Element %d has %d no slip nodes but expecting 2.", sindex, count);
                                }
                            }
                            if(count < 2) FatalError("Not enough no slip nodes");
                            /* Get local index of surface normal vectors */
                            lbdnodes[0] = GetLocalIndex_SurfaceMesh(bdnodes[0]);
                            lbdnodes[1] = GetLocalIndex_SurfaceMesh(bdnodes[1]);
                            /* Get average normal */
                            TwoVectorMean(Vel_SurfaceMeshInwardNormals[lbdnodes[0]], Vel_SurfaceMeshInwardNormals[lbdnodes[1]], abn);
                            
                            /* Get coordinates of boundary nodes */
                            for(i=0; i<3; i++) {
                                bx1[i] = Vel_MeshNodeArray[bdnodes[0]][i];
                                bx2[i] = Vel_MeshNodeArray[bdnodes[1]][i];
                            }
                            /* Get distance from sphere point to boundary edge */
                            d = distline(bx1, bx2, Xs);
                            /* Check for acutal collision (tolerance set to 15% of particle radius) */
                            if(d <= 0.15*Particle_Radius) {
                                /* Compute normal component and reverse normal component if needed */
                                if(vdot(abn, pt->V, Dimensions) < 0) {
                                    for(i=0; i<3; i++) {
                                        vn[i] = vdot(abn, pt->V, Dimensions)*abn[i];
                                        vt[i] = pt->V[i] - vn[i];
                                    }
                                    for(i=0; i<3; i++)
                                        pt->V[i] = vt[i] - Int_NormalFlowScaling*vn[i]; /* Int_NormalFlowScaling used as restitution coefficient */
                                }
                                break; /* If collision detected, no need to check other points */
                            }
                        }
                        else if(Vel_MeshElementFlagArray[sindex] == 3) { /* 3 boundary nodes */
                            /* Determine which nodes have zero velocity */
                            for(i=0; i<(Dimensions+1); i++) {
                                u = Vel_UnstructVelArray_U[0][Vel_MeshElementArray[sindex].Nodes[i]];
                                v = Vel_UnstructVelArray_V[0][Vel_MeshElementArray[sindex].Nodes[i]];
                                w = Vel_UnstructVelArray_W[0][Vel_MeshElementArray[sindex].Nodes[i]];
                                if((u*u + v*v + w*w) < TINY) {
                                    /* Save off index for the no slip node */
                                    bdnodes[count] = Vel_MeshElementArray[sindex].Nodes[i];
                                    count++;
                                    if(count > 3) FatalError("Element %d has %d no slip nodes but expecting %d.", sindex, count, 3);
                                }
                            }
                            if(count < 3) FatalError("Not enough no slip nodes");
                            /* Get local index of surface normal vectors */
                            lbdnodes[0] = GetLocalIndex_SurfaceMesh(bdnodes[0]);
                            lbdnodes[1] = GetLocalIndex_SurfaceMesh(bdnodes[1]);
                            lbdnodes[2] = GetLocalIndex_SurfaceMesh(bdnodes[2]);
                            /* Get average normal */
                            ThreeVectorMean(Vel_SurfaceMeshInwardNormals[lbdnodes[0]], Vel_SurfaceMeshInwardNormals[lbdnodes[1]], Vel_SurfaceMeshInwardNormals[lbdnodes[2]], abn) ;
                            /* Get coordinates of a boundary node */
                            for(i=0; i<3; i++) {
                                bx1[i] = Vel_MeshNodeArray[bdnodes[0]][i];
                                /*bx2[i] = Vel_MeshNodeArray[bdnodes[1]][i];
                                 bx3[i] = Vel_MeshNodeArray[bdnodes[2]][i];*/
                            }
                            /* Get distance from sphere point to boundary face */
                            d = fabs(vdot(abn, Xs, Dimensions) - vdot(abn, bx1, Dimensions));
                            /* Check for actual collision (tolerance set to 15% of particle radius) */
                            if(d <= 0.15*Particle_Radius) {
                                /* Compute normal component and reverse normal component if needed */
                                if(vdot(abn, pt->V, Dimensions) < 0) {
                                    for(i=0; i<3; i++) {
                                        vn[i] = vdot(abn, pt->V, Dimensions)*abn[i];
                                        vt[i] = pt->V[i] - vn[i];
                                    }
                                    for(i=0; i<3; i++)
                                        pt->V[i] = vt[i] - Int_NormalFlowScaling*vn[i]; /* Int_NormalFlowScaling used as restitution coefficient */
                                }
                                break; /* If collision detected, no need to check other points */
                            }
                        }
                        else if(Vel_MeshElementFlagArray[sindex] > 3)
                            FatalError("Element has all no slip nodes");
                    } /* end if(Vel_MeshElementFlagArray[sindex] == 1) */
                } /* end if(sindex >= 0) */
            } /* end loop over sphere boundary */
        } /* end if(Int_NormalFlow) */
        
        /* time derivative of fluid velocity at time tc and particle position x(tc) */
        GetVelocity(tc + h/2, pt, vec2);
        if((tc - h/2) < Data_LoadedTMin) { /* Use forward difference */
            GetVelocity(tc, pt, vec1);
            for(ii = 0; ii < 3; ii++)
                dotu0[ii] = (vec2[ii] - vec1[ii]) / (h / 2);
        }
        else { /* Use central difference */
            GetVelocity(tc - h/2, pt, vec1);
            for(ii = 0; ii < 3; ii++)
                dotu0[ii] = (vec2[ii] - vec1[ii]) / h;
        }
        
        /* spatial derivative of fluid velocity at time tc and particle position x(tc) */
        /* interpolation points */
        for(ii=0;ii<6;ii++) {
            MPA[ii].LeftDomain = pt->LeftDomain;
            MPA[ii].ElementIndex = pt->ElementIndex;
        }
        MPA[0].X[0] = pt->X[0] - dx/2;
        MPA[0].X[1] = pt->X[1];
        MPA[0].X[2] = pt->X[2];
        MPA[1].X[0] = pt->X[0] + dx/2;
        MPA[1].X[1] = pt->X[1];
        MPA[1].X[2] = pt->X[2];
        MPA[2].X[0] = pt->X[0];
        MPA[2].X[1] = pt->X[1] - dx/2;
        MPA[2].X[2] = pt->X[2];
        MPA[3].X[0] = pt->X[0];
        MPA[3].X[1] = pt->X[1] + dx/2;
        MPA[3].X[2] = pt->X[2];
        MPA[4].X[0] = pt->X[0];
        MPA[4].X[1] = pt->X[1];
        MPA[4].X[2] = pt->X[2] - dx/2;
        MPA[5].X[0] = pt->X[0];
        MPA[5].X[1] = pt->X[1];
        MPA[5].X[2] = pt->X[2] + dx/2;
        
        /* fluid velocity at interpolation points */
        for(ii = 0; ii < 6; ii++)
            GetVelocity(tc, &MPA[ii], MPA[ii].V);
        
        /* compute gradient */
        for(jj=0; jj < 6; jj = jj + 2)
            if(MPA[jj].LeftDomain && MPA[jj+1].LeftDomain)
                for(ii = 0; ii < 3; ii++)
                    gradu0[ii][jj/2] = 0;
            else if(MPA[jj].LeftDomain) /* Forward difference */
                for(ii = 0; ii < 3; ii++)
                    gradu0[ii][jj/2] = (MPA[jj+1].V[ii] - u0[ii]) / (dx / 2);
            else if(MPA[jj+1].LeftDomain) /* Backward difference */
                for(ii = 0; ii < 3; ii++)
                    gradu0[ii][jj/2] = (u0[ii] - MPA[jj].V[ii]) / (dx / 2);
            else /* Central differencing */
                for(ii = 0; ii < 3; ii++)
                    gradu0[ii][jj/2] = (MPA[jj+1].V[ii] - MPA[jj].V[ii]) / dx;
        
        /* material derivative of fluid at time tc and particle position x(tc) */
        for(ii = 0; ii < 3; ii++)
            dudt0[ii] = dotu0[ii] + (u0[0] * gradu0[ii][0]) + (u0[1] * gradu0[ii][1]) + (u0[2] * gradu0[ii][2]);
        
        /* Difference in velocity between particle and fluid  */
        for(ii = 0; ii < 3; ii++)
            wp[ii] = exp(-K * R *h) * ( (pt->V[ii] - u0[ii]) - ((1 - 1.5 * R) * (Gravity[ii] - dudt0[ii]) / (K * R)) ) + ( (1 - 1.5 * R) * (Gravity[ii] - dudt0[ii]) / (K * R) );
        
        /* Update velocity and position of particle */
        for(ii = 0; ii < 3; ii++) {
            pt->V[ii] = u0[ii] + wp[ii];
            pt->X[ii] = pt->X[ii] + pt->V[ii] * h;
            if(Data_MeshType == CARTESIAN && Data_XPeriodic == SPATIALPERIODIC) ReMapPt(pt);
        }
        
        /* update time step */
        tc = tc + h;
        
        /* Update which element in velocity mesh contains point */
        if(Data_MeshType == UNSTRUCTURED)
            pt->ElementIndex = Get_Element_Local_Search(pt->X, pt->ElementIndex);
        
        /* Check if point left domain */
        /* if((Data_MeshType == CARTESIAN && TestOutsideDomain(pt->X)) || (Data_MeshType == UNSTRUCTURED && pt->ElementIndex < 0)) { */
        if((Data_MeshType == CARTESIAN && (TestOutsideDomain(pt->X) || TestOutsideCartVelDomain(pt->X)) && !Data_XPeriodic)
           || (Data_MeshType == UNSTRUCTURED && (pt->ElementIndex < 0 || TestOutsideDomain(pt->X)))) {
            /* Point left domain */
            pt->LeftDomain = 1;
            pt->LeftDomainTime = tc - h;
            /* Revert to previous position and break */
            for(ii = 0; ii < 3; ii++)
                pt->X[ii] = pt->X[ii] - pt->V[ii] * h;
            return pt->LeftDomainTime;
        }
        else { /* Point still in domain */
            if(Trace_Compute) {
                if(Trace_CETCompute) { /* If running exposure time computations */
                    /* Compute previous position */
                    for(ii = 0; ii < 3; ii++)
                        X0[0] = pt->X[ii] - pt->V[ii] * h;
                    if(Trace_CETAuxillaryMesh) { /* If using an auxillary mesh for exposure time computations */
                        /* Update which element in auxillary mesh contains point */
                        pt->AuxElementIndex = Get_Element_Local_Search_Aux(pt->X, pt->AuxElementIndex);
                        /* Check if point outside of any element */
                        if(pt->AuxElementIndex < 0) {
                            pt->LeftDomain = 1;
                            pt->LeftDomainTime = tc - h;
                            /* Return previous position and break */
                            for(ii = 0; ii < 3; ii++)
                                pt->X[ii] = X0[ii];
                            return pt->LeftDomainTime;
                        }
                        else /* Compute exposure time using Aux mesh */
                            ComputeExposureTime(X0, pt->X, pt->AuxElementIndex, h);
                    }
                    else /* Compute exposure time using Vel mesh */
                        ComputeExposureTime(X0, pt->X, pt->ElementIndex, h);
                } /* end if running exposure time computations */
                if((Trace_ReleaseStrategy != STAGGERED) && Trace_APCompute) {
                    if(Data_MeshType == CARTESIAN)
                        strain_rate1 = GetStrainRate_Cartesian(tc, pt);
                    else if(Data_MeshType == UNSTRUCTURED)
                        strain_rate1 = GetStrainRate_Unstructured(tc, pt);
                    pt->Scalar += (h * ((strain_rate0 + strain_rate1) / 2));
                }
                if(Trace_VorticityCompute)
                    SetVorticity(tc, pt);
            } /* if(Trace_Compute) */
        } /* Point still in domain */
    } /* Integration not complete */
    
    return tend;
    
} /* End of function pEuler() */

double Euler(LagrangianPoint *pt, double tstart, double tend) {
    /***
     Euler integration routine of x' = v
     ***/
    
    double h, tc, strain_rate0 = 0, strain_rate1 = 0, ve[3], X0[3];
    
    if(pt->LeftDomain)
        FatalError("Attempting to integrate point that already left domain.");
    if(Data_MeshType == UNSTRUCTURED && pt->ElementIndex < 0)
        FatalError("Attempting to integrate point with Element_Index = -1");
    
    /* Initialize local time variable */
    tc = tstart;
    /* Make time step negative if tend < tstart */
    h = SIGN(Int_TimeStep, tend - tstart);
    
    /* Intergrate until tc reaches tend */
    while((tc - tend) * (tend - tstart) < 0.0 ) { /* Integration not complete */
        /* If trying to integrate past tend, reset step size, h, to integrate until tend */
        if ((tc + h - tend) * (tc + h - tstart) > 0.0)
            h = SIGN(tend - tc, tend - tstart);
        
        /* Compute strain rate */
        if(Trace_Compute && (Trace_ReleaseStrategy != STAGGERED) && Trace_APCompute)
            strain_rate0 = GetStrainRate(tc, pt);
        
        /* Get velocity at location pt.X at time tc */
        GetVelocity(tc, pt, ve);
        
        /* Update position, time, and (if needed) ElementIndex */
        pt->X[0] = pt->X[0] + h * ve[0];
        pt->X[1] = pt->X[1] + h * ve[1];
        pt->X[2] = pt->X[2] + h * ve[2];
        if(Data_MeshType == CARTESIAN && Data_XPeriodic == SPATIALPERIODIC) ReMapPt(pt);
        tc = tc + h;
        if(Data_MeshType == UNSTRUCTURED)
        /* Update which element in velocity mesh contains point */
            pt->ElementIndex = Get_Element_Local_Search(pt->X, pt->ElementIndex);
        
        /* Check if point left domain */
        if((Data_MeshType == CARTESIAN && (TestOutsideDomain(pt->X) || TestOutsideCartVelDomain(pt->X)) && !Data_XPeriodic)
           || (Data_MeshType == UNSTRUCTURED && (pt->ElementIndex < 0 || TestOutsideDomain(pt->X)))) {
            /* Point left domain */
            pt->LeftDomain = 1;
            pt->LeftDomainTime = tc - h;
            /* Return previous position and break */
            pt->X[0] = pt->X[0] - h * ve[0];
            pt->X[1] = pt->X[1] - h * ve[1];
            pt->X[2] = pt->X[2] - h * ve[2];
            return pt->LeftDomainTime;
        }
        else { /* Point still in domain */
            /* Update velocity of particle */
            if(Int_Extrapolate)
                GetVelocity(tc, pt, pt->V);
            if(Trace_Compute) {
                if(Trace_CETCompute) { /* If running exposure time computations */
                    /* Compute previous position */
                    X0[0] = pt->X[0] - h * ve[0];
                    X0[1] = pt->X[1] - h * ve[1];
                    X0[2] = pt->X[2] - h * ve[2];
                    if(Trace_CETAuxillaryMesh) { /* If using an auxillary mesh for exposure time computations */
                        /* Update which element in auxillary mesh contains point */
                        pt->AuxElementIndex = Get_Element_Local_Search_Aux(pt->X, pt->AuxElementIndex);
                        /* Check if point outside of any element */
                        if(pt->AuxElementIndex < 0) {
                            pt->LeftDomain = 1;
                            pt->LeftDomainTime = tc - h;
                            /* Return previous position and break */
                            pt->X[0] = X0[0];
                            pt->X[1] = X0[1];
                            pt->X[2] = X0[2];
                            return pt->LeftDomainTime;
                        }
                        else /* Compute exposure time using Aux mesh */
                            ComputeExposureTime(X0, pt->X, pt->AuxElementIndex, h);
                    }
                    else /* Compute exposure time using Vel mesh */
                        ComputeExposureTime(X0, pt->X, pt->ElementIndex, h);
                } /* end if running exposure time computations */
                if((Trace_ReleaseStrategy != STAGGERED) && Trace_APCompute) {
                    strain_rate1 = GetStrainRate(tc, pt);
                    pt->Scalar += (h * ((strain_rate0 + strain_rate1) / 2));
                }
                if(Trace_VorticityCompute)
                    SetVorticity(tc, pt);
            } /* If Trace_Compute */
        } /* Point still in domain */
    } /* Integration not complete */
    
    return tend;
    
} /* End of function Euler() */



double RK4(LagrangianPoint *pt, double tstart, double tend) {
    /***
     Explicit 4th order Runge Kutta integration routine
     ***/
    
    double tc, h, strain_rate0 = 0, strain_rate1 = 0, k1[3], k2[3], k3[3], k4[3], X0[3];
    LagrangianPoint MP1, MP2, MP3, MP4;
    
    if(pt->LeftDomain)
        FatalError("Attempting to integrate point that already left domain");
    else {
        MP1.LeftDomain = pt->LeftDomain;
        MP2.LeftDomain = pt->LeftDomain;
        MP3.LeftDomain = pt->LeftDomain;
        MP4.LeftDomain = pt->LeftDomain;
    }
    if(Data_MeshType == UNSTRUCTURED) {
        if(pt->ElementIndex < 0)
            FatalError("Attempting to integrate point with Element_Index = -1");
        else {
            MP1.ElementIndex = pt->ElementIndex;
            MP2.ElementIndex = pt->ElementIndex;
            MP3.ElementIndex = pt->ElementIndex;
            MP4.ElementIndex = pt->ElementIndex;
        }
    }
    
    /* Initialize time */
    tc = tstart;
    h = SIGN(Int_TimeStep, tend - tstart); /* Make time step negative if tend < tstart */
    
    /* Intergrate until tc reaches tend */
    while((tc - tend) * (tend - tstart) < 0.0 ) {
        /* If trying to integrate past tend, reset step size to integrate until tend */
        if((tc + h - tend) * (tc + h - tstart) > 0.0)
            h = SIGN(tend - tc, tend - tstart);
        
        /* Compute strain rate */
        if(Trace_Compute && (Trace_ReleaseStrategy != STAGGERED) && Trace_APCompute) {
            if(Data_MeshType == CARTESIAN)
                strain_rate0 = GetStrainRate_Cartesian(tc, pt);
            else if(Data_MeshType == UNSTRUCTURED)
                strain_rate0 = GetStrainRate_Unstructured(tc, pt);
        }
        
        /* k1 */
        MP1.X[0] = pt->X[0];
        MP1.X[1] = pt->X[1];
        MP1.X[2] = pt->X[2];
        GetVelocity(tc, &MP1, k1);
        
        /* k2 */
        MP2.X[0] = pt->X[0] + 0.5 * k1[0] * h;
        MP2.X[1] = pt->X[1] + 0.5 * k1[1] * h;
        MP2.X[2] = pt->X[2] + 0.5 * k1[2] * h;
        if(Data_MeshType == CARTESIAN && Data_XPeriodic == SPATIALPERIODIC) ReMapPt(&MP2);
        GetVelocity(tc + h/2, &MP2, k2);
        
        /* k3 */
        MP3.X[0] = pt->X[0] + 0.5 * k2[0] * h;
        MP3.X[1] = pt->X[1] + 0.5 * k2[1] * h;
        MP3.X[2] = pt->X[2] + 0.5 * k2[2] * h;
        if(Data_MeshType == CARTESIAN && Data_XPeriodic == SPATIALPERIODIC) ReMapPt(&MP3);
        GetVelocity(tc + h/2, &MP3, k3);
        
        /* k4 */
        MP4.X[0] = pt->X[0] + k3[0] * h;
        MP4.X[1] = pt->X[1] + k3[1] * h;
        MP4.X[2] = pt->X[2] + k3[2] * h;
        if(Data_MeshType == CARTESIAN && Data_XPeriodic == SPATIALPERIODIC) ReMapPt(&MP4);
        GetVelocity(tc + h, &MP4, k4);
        
        /* Update position, time, and (if needed) element index */
        pt->X[0] = pt->X[0] + h * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6;
        pt->X[1] = pt->X[1] + h * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6;
        pt->X[2] = pt->X[2] + h * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6;
        tc = tc + h;
        if(Data_MeshType == UNSTRUCTURED)
            pt->ElementIndex = Get_Element_Local_Search(pt->X, pt->ElementIndex);
        if(Data_MeshType == CARTESIAN && Data_XPeriodic == SPATIALPERIODIC) ReMapPt(pt);
        
        /* Check if point left domain */
        if(((Data_MeshType == CARTESIAN && (TestOutsideDomain(pt->X) || TestOutsideCartVelDomain(pt->X) || MP1.LeftDomain || MP2.LeftDomain || MP3.LeftDomain || MP4.LeftDomain) && !Data_XPeriodic)) || (Data_MeshType == UNSTRUCTURED && (TestOutsideDomain(pt->X) || pt->ElementIndex < 0 || MP1.LeftDomain || MP2.LeftDomain || MP3.LeftDomain || MP4.LeftDomain))) {
            /* Point left domain, or data used to update location came from outside domain */
            pt->LeftDomain = 1;
            pt->LeftDomainTime = tc - h;
            /* Compute previous position */
            pt->X[0] = pt->X[0] - h * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6;
            pt->X[1] = pt->X[1] - h * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6;
            pt->X[2] = pt->X[2] - h * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6;
            return pt->LeftDomainTime;
        }
        else { /* Point still in domain */
            /* Update velocity of particle */
            if(Int_Extrapolate) {
                GetVelocity(tc, pt, pt->V);
            }
            if(Trace_Compute) {
                if(Trace_CETCompute) {
                    /* Compute previous position */
                    X0[0] = pt->X[0] - h * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6;
                    X0[1] = pt->X[1] - h * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6;
                    X0[2] = pt->X[2] - h * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6;
                    if(Trace_CETCompute && Trace_CETAuxillaryMesh) {
                        pt->AuxElementIndex = Get_Element_Local_Search_Aux(pt->X, pt->AuxElementIndex);
                        if(pt->AuxElementIndex < 0) {
                            pt->LeftDomain = 1;
                            pt->LeftDomainTime = tc - h;
                            /* Keep initial location */
                            pt->X[0] = X0[0];
                            pt->X[1] = X0[1];
                            pt->X[2] = X0[2];
                            return pt->LeftDomainTime;
                        }
                        else
                            ComputeExposureTime(X0, pt->X, pt->AuxElementIndex, h);
                    }
                    else
                        ComputeExposureTime(X0, pt->X, pt->ElementIndex, h);
                }
                if((Trace_ReleaseStrategy != STAGGERED) && Trace_APCompute) {
                    if(Data_MeshType == CARTESIAN)
                        strain_rate1 = GetStrainRate_Cartesian(tc, pt);
                    else if(Data_MeshType == UNSTRUCTURED)
                        strain_rate1 = GetStrainRate_Unstructured(tc, pt); /* Note, tc is already updated */
                    pt->Scalar += (h * ((strain_rate0 + strain_rate1) / 2));
                }
                if(Trace_VorticityCompute)
                    SetVorticity(tc, pt);
            } /* end if(Trace_Compute) */
        }
    }
    
    return tend;
    
}

double RKF(LagrangianPoint *pt, double tstart, double tend) {
    
    double k1[3], k2[3], k3[3], k4[3], k5[3], k6[3], X0[3];
    double a3, b3, a4, b4, c4, a5, b5, c5, d5, a6, b6, c6, d6, e6, a7, b7, c7, d7, e7, a8, b8, c8, d8;
    double h, tol, err, errx, erry, errz, ss, tc, ts, hmax, strain_rate0 = 0, strain_rate1 = 0;
    LagrangianPoint MP1, MP2, MP3, MP4, MP5, MP6;
    
    if(pt->LeftDomain)
        FatalError("Attempting to integrate point that already left domain");
    else {
        MP1.LeftDomain = pt->LeftDomain;
        MP2.LeftDomain = pt->LeftDomain;
        MP3.LeftDomain = pt->LeftDomain;
        MP4.LeftDomain = pt->LeftDomain;
        MP5.LeftDomain = pt->LeftDomain;
        MP6.LeftDomain = pt->LeftDomain;
    }
    
    if(Data_MeshType == UNSTRUCTURED) {
        if(pt->ElementIndex == -1)
            FatalError("Attempting to integrate point with Element_Index = -1");
        else {
            MP1.ElementIndex = pt->ElementIndex;
            MP2.ElementIndex = pt->ElementIndex;
            MP3.ElementIndex = pt->ElementIndex;
            MP4.ElementIndex = pt->ElementIndex;
            MP5.ElementIndex = pt->ElementIndex;
            MP6.ElementIndex = pt->ElementIndex;
        }
    }
    
    /* Set RKF coefficients */
    a3 =     3.0 / 32.0;
    b3 =     9.0 / 32.0;
    a4 =  1932.0 / 2197.0;
    b4 = -7200.0 / 2197.0;
    c4 =  7296.0 / 2197.0;
    a5 =   439.0 / 216.0;
    b5 =    -8.0;
    c5 =  3680.0 / 513.0;
    d5 =  -845.0 / 4104.0;
    a6 =    -8.0 / 27.0;
    b6 =     2.0;
    c6 = -3544.0 / 2565.0;
    d6 =  1859.0 / 4104.0;
    e6 =   -11.0 / 40.0;
    a7 =     1.0 / 360.0;
    b7 =  -128.0 / 4275.0;
    c7 = -2197.0 / 75240.0;
    d7 =     1.0 / 50.0;
    e7 =     2.0 / 55.0;
    a8 =    25.0 / 216.0;
    b8 =  1408.0 / 2565.0;
    c8 =  2197.0 / 4104.0;
    d8 =    -1.0 / 5.0;
    
    /* Timing parameters */
    h    = SIGN(Int_MaxTimeStep, tend - tstart);
    hmax = Int_MaxTimeStep;
    tol  = Int_Accuracy;
    tc   = tstart;
    
    /* Intergrate until tc reaches tend */
    while((tc - tend) * (tend - tstart) < 0.0 ) {
        /* If trying to integrate longer than needed, reset step size to integrate until tend */
        if ((tc + h - tend) * (tc + h - tstart) > 0.0) {
            h = SIGN(tend - tc, tend - tstart);
        }
        
        /* Compute strain rate */
        if(Trace_Compute && (Trace_ReleaseStrategy != STAGGERED) && Trace_APCompute) {
            if(Data_MeshType == CARTESIAN)
                strain_rate0 = GetStrainRate_Cartesian(tc, pt);
            else if(Data_MeshType == UNSTRUCTURED)
                strain_rate0 = GetStrainRate_Unstructured(tc, pt);
        }
        
        /* k1 */
        MP1.X[0] = pt->X[0];
        MP1.X[1] = pt->X[1];
        MP1.X[2] = pt->X[2];
        GetVelocity(tc, &MP1, k1);
        
        /* k2 */
        MP2.X[0] = pt->X[0] + 0.25 * k1[0] * h;
        MP2.X[1] = pt->X[1] + 0.25 * k1[1] * h;
        MP2.X[2] = pt->X[2] + 0.25 * k1[2] * h;
        if(Data_MeshType == CARTESIAN && Data_XPeriodic == SPATIALPERIODIC) ReMapPt(&MP2);
        ts = tc + 0.25 * h;
        GetVelocity(ts, &MP2, k2);
        
        /* k3 */
        MP3.X[0] = pt->X[0] + (a3 * k1[0] + b3 * k2[0]) * h;
        MP3.X[1] = pt->X[1] + (a3 * k1[1] + b3 * k2[1]) * h;
        MP3.X[2] = pt->X[2] + (a3 * k1[2] + b3 * k2[2]) * h;
        if(Data_MeshType == CARTESIAN && Data_XPeriodic == SPATIALPERIODIC) ReMapPt(&MP3);
        ts = tc + 0.375 * h;
        GetVelocity(ts, &MP3, k3);
        
        /* k4 */
        MP4.X[0] = pt->X[0] + (a4 * k1[0] + b4 * k2[0] + c4 * k3[0]) * h;
        MP4.X[1] = pt->X[1] + (a4 * k1[1] + b4 * k2[1] + c4 * k3[1]) * h;
        MP4.X[2] = pt->X[2] + (a4 * k1[2] + b4 * k2[2] + c4 * k3[2]) * h;
        if(Data_MeshType == CARTESIAN && Data_XPeriodic == SPATIALPERIODIC) ReMapPt(&MP4);
        ts = tc + (12.0 * h) / 13.0;
        GetVelocity(ts, &MP4, k4);
        
        /* k5 */
        MP5.X[0] = pt->X[0] + (a5 * k1[0] + b5 * k2[0] + c5 * k3[0] + d5 * k4[0]) * h;
        MP5.X[1] = pt->X[1] + (a5 * k1[1] + b5 * k2[1] + c5 * k3[1] + d5 * k4[1]) * h;
        MP5.X[2] = pt->X[2] + (a5 * k1[2] + b5 * k2[2] + c5 * k3[2] + d5 * k4[2]) * h;
        if(Data_MeshType == CARTESIAN && Data_XPeriodic == SPATIALPERIODIC) ReMapPt(&MP5);
        ts = tc + h;
        GetVelocity(ts, &MP5, k5);
        
        /* k6 */
        MP6.X[0] = pt->X[0] + (a6 * k1[0] + b6 * k2[0] + c6 * k3[0] + d6 * k4[0] + e6 * k5[0]) * h;
        MP6.X[1] = pt->X[1] + (a6 * k1[1] + b6 * k2[1] + c6 * k3[1] + d6 * k4[1] + e6 * k5[1]) * h;
        MP6.X[2] = pt->X[2] + (a6 * k1[2] + b6 * k2[2] + c6 * k3[2] + d6 * k4[2] + e6 * k5[2]) * h;
        if(Data_MeshType == CARTESIAN && Data_XPeriodic == SPATIALPERIODIC) ReMapPt(&MP6);
        ts = tc + 0.5 * h;
        GetVelocity(ts, &MP6, k6);
        
        /* Check if any intermediate interpolation points were outside the domain */
        if(((Data_MeshType == CARTESIAN && (MP1.LeftDomain || MP2.LeftDomain || MP3.LeftDomain || MP4.LeftDomain || MP5.LeftDomain || MP6.LeftDomain) && !Data_XPeriodic)) || (Data_MeshType == UNSTRUCTURED && (MP1.LeftDomain || MP2.LeftDomain || MP3.LeftDomain || MP4.LeftDomain || MP5.LeftDomain || MP6.LeftDomain))) {
            pt->LeftDomain = 1;
            pt->LeftDomainTime = tc;
            return pt->LeftDomainTime;
        }
        
        /* Error between RK4 and RK5 methods */
        errx = fabs(a7 * k1[0] + b7 * k3[0] + c7 * k4[0] + d7 * k5[0] + e7 * k6[0]);
        erry = fabs(a7 * k1[1] + b7 * k3[1] + c7 * k4[1] + d7 * k5[1] + e7 * k6[1]);
        errz = fabs(a7 * k1[2] + b7 * k3[2] + c7 * k4[2] + d7 * k5[2] + e7 * k6[2]);
        err = fmax(fmax(errx, erry), errz); /* Maximum error */
        
        /* Set scaling factor to determine optimal step size (NOTE: choice for updating step size is conventional, but somewhat arbitrary) */
        if (err == 0.0) {
            ss = 0.0;
        }
        else {
            ss = pow(fabs((tol * h) / (2 * err)), 0.25);
        }
        
        if (err < tol || fabs(h) <= Int_MinTimeStep) {  /* Below error tolerance or min time step reached, accept RK4 approximation */
            /* Update position, time, and (if needed) element index */
            pt->X[0] = pt->X[0] + (a8 * k1[0] + b8 * k3[0] + c8 * k4[0] + d8 * k5[0]) * h;
            pt->X[1] = pt->X[1] + (a8 * k1[1] + b8 * k3[1] + c8 * k4[1] + d8 * k5[1]) * h;
            pt->X[2] = pt->X[2] + (a8 * k1[2] + b8 * k3[2] + c8 * k4[2] + d8 * k5[2]) * h;
            if(Data_MeshType == CARTESIAN && Data_XPeriodic == SPATIALPERIODIC) ReMapPt(pt);
            tc = tc + h;
            if(Data_MeshType == UNSTRUCTURED)
                pt->ElementIndex = Get_Element_Local_Search(pt->X, pt->ElementIndex);
            /* Check if point left domain */
            if((Data_MeshType == CARTESIAN && (TestOutsideDomain(pt->X) || TestOutsideCartVelDomain(pt->X)) && !Data_XPeriodic)
               || (Data_MeshType == UNSTRUCTURED && (pt->ElementIndex < 0 || TestOutsideDomain(pt->X)))) {
                /* Point left domain or data used to update location came from outside domain */
                pt->LeftDomain = 1;
                pt->LeftDomainTime = tc - h;
                /* Compute previous position */
                pt->X[0] = pt->X[0] - (a8 * k1[0] + b8 * k3[0] + c8 * k4[0] + d8 * k5[0]) * h;
                pt->X[1] = pt->X[1] - (a8 * k1[1] + b8 * k3[1] + c8 * k4[1] + d8 * k5[1]) * h;
                pt->X[2] = pt->X[2] - (a8 * k1[2] + b8 * k3[2] + c8 * k4[2] + d8 * k5[2]) * h;
                return pt->LeftDomainTime;
            }
            else { /* Point still in domain */
                /* Update velocity of particle */
                if(Int_Extrapolate) {
                    GetVelocity(tc, pt, pt->V);
                }
                if(Trace_Compute) {
                    if(Trace_CETCompute) {
                        /* Compute previous position */
                        X0[0] = pt->X[0] - (a8 * k1[0] + b8 * k3[0] + c8 * k4[0] + d8 * k5[0]) * h;
                        X0[1] = pt->X[1] - (a8 * k1[1] + b8 * k3[1] + c8 * k4[1] + d8 * k5[1]) * h;
                        X0[2] = pt->X[2] - (a8 * k1[2] + b8 * k3[2] + c8 * k4[2] + d8 * k5[2]) * h;
                        if(Trace_CETCompute && Trace_CETAuxillaryMesh) {
                            pt->AuxElementIndex = Get_Element_Local_Search_Aux(pt->X, pt->AuxElementIndex);
                            if(pt->AuxElementIndex < 0) {
                                pt->LeftDomain = 1;
                                pt->LeftDomainTime = tc - h;
                                /* Compute previous position */
                                pt->X[0] = X0[0];
                                pt->X[1] = X0[1];
                                pt->X[2] = X0[2];
                                return pt->LeftDomainTime;
                            }
                            else
                                ComputeExposureTime(X0, pt->X, pt->AuxElementIndex, h);
                        }
                        else
                            ComputeExposureTime(X0, pt->X, pt->ElementIndex, h);
                    }
                    if((Trace_ReleaseStrategy != STAGGERED) &&Trace_APCompute) {
                        if(Data_MeshType == CARTESIAN)
                            strain_rate1 = GetStrainRate_Cartesian(tc, pt);
                        else if(Data_MeshType == UNSTRUCTURED)
                            strain_rate1 = GetStrainRate_Unstructured(tc, pt);
                        pt->Scalar += (h * ((strain_rate0 + strain_rate1) / 2));
                    }
                    if(Trace_VorticityCompute)
                        SetVorticity(tc, pt);
                } /* end if(Trace_Compute) */
                if(fabs(h) > Int_MinTimeStep && ss > 1.0) { /* Update time step size if needed */
                    ss = (ss > 4.0)? 4.0:ss;
                    /* Set optimal time step */
                    h = SIGN((fabs(ss * h) > hmax)? hmax : ss * h, tend - tstart);
                }
            }
        }
        else { /* Do not accept, decrease step size */
            if(ss < 1.0) {
                ss = (ss < 0.1)? 0.1:ss;
                /* Set optimal time step */
                h = SIGN((fabs(ss * h) < Int_MinTimeStep)? Int_MinTimeStep : ss * h, tend - tstart);
            }
            else
                h = SIGN((fabs(0.5 * h) < Int_MinTimeStep)? Int_MinTimeStep : 0.5 * h, tend - tstart);
        }
    }
    
    return tend;
    
} /* End of RKF */


