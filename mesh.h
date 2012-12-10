/*
 *  mesh.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#ifndef INC_MESH_H
#define INC_MESH_H

void LoadMeshData(void);
void FreeMeshData(void);
void LoadCartMeshData(void);
void LoadUnstructMeshData(void);
void LoadBoundaryElementFlags(void);
void LoadMeshSurfaceNormals(void);
int Get_Element_Global_Search(const double *X);
int Get_Element_Global_Search_Aux(const double *X);
int Get_Element_Local_Search(const double *X, int guess);
int Get_Element_Local_Search_Aux(const double *X, int guess);

#endif