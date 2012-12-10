/*
 *  integration.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#ifndef INC_INTEGRATION_H
#define INC_INTEGRATION_H

void Advect(LagrangianPoint *pt, double tstart, double tend);
double pEuler(LagrangianPoint *pt, double tstart, double tend);
double Euler(LagrangianPoint *pt, double tstart, double tend);
double RK4(LagrangianPoint *pt, double tstart, double tend);
double RKF(LagrangianPoint *pt, double tstart, double tend);

#endif