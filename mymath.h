/*
 *  mymath.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */


#ifndef INC_MYMATH_H
#define INC_MYMATH_H

double vdot(double *v1, double *v2, int dim);
void cross(const double *u, const double *v, double *result);
double dist(const double *v1, const double *v2, const int dim);
void vdiff(const double *v1, const double *v2, double *result);
double distline(const double *x1, const double *x2, const double *x0);
void TwoVectorMean(const double *v1, const double *v2, double *m);
void ThreeVectorMean(const double *v1, const double *v2, const double *v3, double *m);
double pythag(double a, double b);
double GetMaxEigenvalue(double a[][3]);

#endif

