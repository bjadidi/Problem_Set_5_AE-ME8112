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
#define nr 4
#define nz 6
#define lr 0.01
#define lz 0.1
#define N nr*nz
long double Vor_1[nr][nz],Vr_1[nr][nz],Vz_1[nr][nz],Vor_2[nr][nz],Vr_2[nr][nz],Vz_2[nr][nz],R[nr],Z[nz];
long double deltar,deltaz,deltat,Vor_0,Vr_0,Vz_0,d,vis,C1_r,C1_z;
int i,j;
fstream output_file ("PS5.xls",ios::out);



//Decleration of Functions
void Mesh_Generation ();
void Initialize ();
void Output_Function ();
void Vor_Middle_Points ();
void Vor_Boundary ();



/*////////////////Begin_#///////////////////////////////////////////////////

/*///////////////End_#//////////////////////////////////////////////////////



////////////////Begin_Main code///////////////////////////////////////////////////

int main()
{
	
vis = 1.4e-4; //viscosity
d = 8.3e-4; //density
deltat = 0.01;

//Initialize

Vor_0 = 0.0;
Vr_0 = 0.0;
Vz_0 = 0.0;
Initialize ();

//Mesh
Mesh_Generation ();


//Vorticity
Vor_Middle_Points ();
//Vor_Boundary ();


//Radial Velocity



//Axial Velocity



//Output_Function
Output_Function ();

return 0;

}

///////////////End_Main code//////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////
//===========================FUNCTIONS===========================================
/////////////////////////////////////////////////////////////////////////////////


////////////////Begin_///////////////////////////////////////////////////

///////////////End_//////////////////////////////////////////////////////

////////////////Begin_Mesh_Generation///////////////////////////////////////////////////

void Mesh_Generation()
{
	R[0] = 0;
	deltar = lr/(nr-1) ;
	for (int i=1;i<nr;i++)
	{
		R[i] = R[i-1] + deltar;
	}
	
	Z[0] = 0;
	deltaz = lz/(nz-1) ;
	for (int j=1;j<nz;j++)
	{
		Z[j] = Z[j-1] + deltaz;
	}
}

///////////////End_Mesh_Generation//////////////////////////////////////////////////////


////////////////Begin_Initialize///////////////////////////////////////////////////
void Initialize ()
{
	for (int i=0;i<nr;i++)
	{
		for (int j=0;j<nz;j++)
		{
			Vor_1[i][j] = Vor_0;
			Vr_1[i][j] = Vr_0;
			Vz_1[i][j] = Vz_0;
		}
	}

}


///////////////End_Initialize//////////////////////////////////////////////////////


////////////////Begin_Vor_Middle_Points///////////////////////////////////////////////////

void Vor_Middle_Points ()
{
	C1_r = (deltat * vis / d)/(deltar * deltar);
	C1_z = (deltat * vis / d)/(deltaz * deltaz);
	for (int i=1;i<nr-1;i++)
	{
		for (int j=1;j<nz-1;j++)
		{
			Vor_2[i][j] = Vor_1[i][j] * (1 - 2*C1_r - 2*C1_z + (deltat * Vr_1[i][j]/R[i])) +\
						  Vor_1[i+1][j]	* (C1_r + (deltat * vis)/(d * 2 * deltar) - (deltat * Vr_1[i][j])/(2*deltar)) +\
						  Vor_1[i-1][j]	* (C1_r - (deltat * vis)/(d * 2 * deltar) + (deltat * Vr_1[i][j])/(2*deltar)) +\
						  Vor_1[i][j+1]	* (C1_z - (deltat * Vz_1[i][j])/(2* deltaz)) +\
						  Vor_1[i][j-1]	* (C1_z + (deltat * Vz_1[i][j])/(2* deltaz)) ;
		}
	}	
	
	
}

///////////////End_Vor_Middle_Points//////////////////////////////////////////////////////

////////////////Begin_Vor_Boundary///////////////////////////////////////////////////

//void Vor_Boundary ()
//{
//// axis of symmetry
//for (int j=0;j<nz;j++)
//{
//	Vor_2[0][j] = 0.0;
//}
//
//// inflow
//for (int i=0;i<nr;i++)
//{
//	Vor_2[i][0] = ((Vr_1[i][j+1] - Vr_1[i][j-1])/(2*deltaz)) - ((Vz_1[i+1][j] - Vz_1[i-1][j])/(2*deltar));
//}
//
////wall
//for (int j=0;j<nz;j++)
//{
//	Vor_2[nr][j] = - ((Vz_1[i+1][j] - Vz_1[i-1][j])/(2*deltar));	
//}	
//
////outflow
//for (int i=0;i<nr;i++)
//{
//	Vor_2[i][nz] = Vor_2[i][nz-1]; 
//}
//	
//}

///////////////End_Vor_Boundary//////////////////////////////////////////////////////


////////////////Begin_///////////////////////////////////////////////////
void Output_Function ()
{
	cout<<endl<<"R ="<<endl;
	for (int i=0;i<nr;i++)
	{
		cout<<R[i]<<"	";
	}
	
	cout<<endl<<"Z ="<<endl;
	for (int j=0;j<nz;j++)
	{
		cout<<Z[j]<<"	";
	}
	
	cout<<endl<<"Vorticity in n+1 ="<<endl;
	for (int i=0;i<nr;i++)
	{
		for (int j=0;j<nz;j++)
		{
			cout<<Vor_2[i][j]<<"	";
		}
		cout<<endl;
	}
}
///////////////End_//////////////////////////////////////////////////////





















