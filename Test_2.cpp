#include<iostream>
#include <math.h>
#include<cmath>
#include <fstream>
using namespace std;
/*////////////////Begin_Header///////////////////////////////////////////////////
================================================================================
// This C++ code was written by Behrooz Jadidi_Student ID:501001145
// The first version was generated on 07_15_2020
// This code was written as a solution to AE/ME8112 - Problem_Set5_Question1
// This code solves vorticity-velocity in a laminar gaseous flow through a cylindrical pipe
================================================================================
/*///////////////End_Header//////////////////////////////////////////////////////

//Global variables
//Geometry
#define nr 3
#define nz 3
#define lr 0.01
#define lz 0.1
#define N nr*nz
long double Vor_1[nr][nz],Vr_1[nr][nz],Vz_1[nr][nz],Vor_2[nr][nz],Vr_2[nr][nz],Vz_2[nr][nz];
long double R[nr],Z[nz],A[N][5],b[N],Velocity[N];
long double v[N],r[N],tt[N],z[N],y[N];
long double deltar,deltaz,deltat,Vor_0,Vr_0,Vz_0,d,vis,C1_r,C1_z,C1,C2;
int i,j;
fstream output_file ("PS5.xls",ios::out);
