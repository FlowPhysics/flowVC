/*
 *  mymath.c
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#include <math.h>
#include "io.h"
#include "macros.h"
#include "mymath.h"

double vdot(double *v1, double *v2, int dim) {
  int i;
  double sum = 0;
  for(i = 0; i < dim; i++) 
    sum = sum + v1[i]*v2[i];
  return(sum); 
}

void cross(const double *u, const double *v, double *result) {
  result[0] = u[1] * v[2] - u[2] * v[1];
  result[1] = u[2] * v[0] - u[0] * v[2];
  result[2] = u[0] * v[1] - u[1] * v[0];
}  

double dist(const double *v1, const double *v2, const int dim) {
  int i;
  double dv[3];
  for(i = 0; i < dim; i++) 
    dv[i] = v1[i] - v2[i];
  return(sqrt(vdot(dv, dv, dim)));
}

void vdiff(const double *v1, const double *v2, double *result) {
  result[0] = v1[0] - v2[0];
  result[1] = v1[1] - v2[1];
  result[2] = v1[2] - v2[2];
}


double distline(const double *x1, const double *x2, const double *x0) {
  /* Computes minimum distance between line spanned by x1 and x2 to the point x0 */

  double d21[3], d10[3], c[3], d;
  
  vdiff(x2, x1, d21);
  vdiff(x1, x0, d10);
  cross(d21,d10,c);
  d = sqrt(vdot(c,c,3)) / sqrt(vdot(d21,d21,3));
  
  return d;

}

void TwoVectorMean(const double *v1, const double *v2, double *m) {
  int i;
  for(i=0; i<3; i++)
    m[i] = (v1[i] + v2[i]) / 2;
}
void ThreeVectorMean(const double *v1, const double *v2, const double *v3, double *m) {
  int i;
  for(i=0; i<3; i++)
    m[i] = (v1[i] + v2[i] + v3[i]) / 3;
}



double pythag(double a, double b) {
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

double GetMaxEigenvalue(double a[][3]) {
  /*Householder reduction of a real, symmetric matrix a[0..n-1][0..n-1]. On output, a is replaced
    by the orthogonal matrix Q effecting the transformation. d[0..n-1] returns the diagonal elements
    of the tridiagonal matrix, and e[0..n-1] the off-diagonal elements, with e[0]=0. Several
    statements, as noted in comments, can be omitted if only eigenvalues are to be found, in which
    case a contains no useful information on output. Otherwise they are to be included.*/
	
  int l,k,j,i;
  double scale,hh,h,g,f;
  int m,iter;
  double s,r,p,dd,c,b,e[3],d[3];
	
  for (i=3-1;i>0;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 0) {
      for (k=0;k<=l;k++)
	scale += fabs(a[i][k]);
      if (scale == 0.0)       /*Skip transformation.*/
	e[i]=a[i][l];
      else {
	for (k=0;k<=l;k++) {
	  a[i][k] /= scale;     /*Use scaled a’s for transformation.*/
	  h += a[i][k]*a[i][k];    /*Form sigma in h.*/
	}
	f=a[i][l];
	g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
	e[i]=scale*g;
	h -= f*g;        /*Now h is equation (11.2.4).*/
	a[i][l]=f-g;       /*Store u in the ith row of a.*/
	f=0.0;
	for (j=0;j<=l;j++) {
	  g=0.0;        /*Form an element of A · u in g.*/
	  for (k=0;k<=j;k++)
	    g += a[j][k]*a[i][k];
	  for (k=j+1;k<=l;k++)
	    g += a[k][j]*a[i][k];
	  e[j]=g/h;       /*Form element of p in temporarily unused element of e.*/
	  f += e[j]*a[i][j];
	}
	hh=f/(h+h);        /*Form K, equation (11.2.11).*/
	for (j=0;j<=l;j++) {      /*Form q and store in e overwriting p.*/
	  f=a[i][j];
	  e[j]=g=e[j]-hh*f;
	  for (k=0;k<=j;k++)     /*Reduce a, equation (11.2.13).*/
	    a[j][k] -= (f*e[k]+g*a[i][k]);
	}
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  e[0]=0.0;
  for (i=0;i<3;i++) {  
    d[i]=a[i][i]; 
  }
	
  for (i=1;i<3;i++) 
    e[i-1]=e[i];        /*Convenient to renumber the elements of e.*/
  e[3-1]=0.0; 
  for (l=0;l<3;l++) {
    iter=0;
    do {
      for (m=l;m<3-1;m++) {     /*Look for a single small subdiagonal element to split the matrix.*/
	dd=fabs(d[m])+fabs(d[m+1]);
	if ((double)(fabs(e[m])+dd) == dd) 
	  break;
      }
      if (m != l) {
	if (iter++ == 30) 
	  FatalError("Too many iterations in tqli");
	g=(d[l+1]-d[l])/(2.0*e[l]);   /*Form shift.*/
	r=pythag(g,1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));  /*This is dm − ks.*/
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) {    /*A plane rotation as in the original QL, followed by Givens rotations to restore tridiagonal form.*/
	  f=s*e[i];
	  b=c*e[i];
	  e[i+1]=(r=pythag(f,g));
	  if (r == 0.0) {     /*Recover from underflow.*/
	    d[i+1] -= p;
	    e[m]=0.0;
	    break;
	  }
	  s=f/r;
	  c=g/r;
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  d[i+1]=g+(p=s*r);
	  g=c*r-b;
	}
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l] = g;
	e[m] = 0.0;
      }
    } while (m != l);
  }
  return( fmax( d[0], fmax(d[1],d[2]) ) ); 
}
