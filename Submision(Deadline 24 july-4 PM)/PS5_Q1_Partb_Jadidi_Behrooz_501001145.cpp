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
#define nr 21
#define nz 201
#define lr 1.0
#define lz 10.0
#define N nr*nz

long double Vor_1[nr][nz],Vr_1[nr][nz],Vz_1[nr][nz],Vorticity[nr][nz],m[nz];
long double R[nr],Z[nz],A[N][5],b[N],x[N];
int t;
long double deltar,deltaz,deltat,d,vis,C1_r,C1_z,C1,C2,clvz,max_m,avg_m,Tol,Tol_limit,sum_m;
fstream output_file_1 ("PS5_Vorticity.xls",ios::out);
fstream output_file_2 ("PS5_Vr.xls",ios::out);
fstream output_file_3 ("PS5_Vz.xls",ios::out);
fstream output_file_4 ("PS5_max_m.xls",ios::out);
fstream output_file_5 ("PS5_avg_m.xls",ios::out);






//Decleration of Functions
void Mesh_Generation ();
void Initialize (int IC);
void Output(string label);
void Vor_Points ();
void Zero_A_b ();
void Vr_Middle_Points ();
void Vr_Boundary ();
void P_Bi_CGSTAB ();
long double dot_product (long double a_1[], long double b_1[]);
void Vz_Points ();
void Change_x_Vr ();
void Update_Vor ();
long double Integration(int j);



////////////////Begin_Main code///////////////////////////////////////////////////

int main()
{
clvz = 	5.0;
vis = 1.4e-4; //viscosity
d = 8.3e-4; //density
deltat = 10e-6;
Tol_limit = 1e-9;

//Mesh
Mesh_Generation ();


//Initialize
Initialize (1);

Tol = 10.0;
t =1;
while (Tol>Tol_limit)
//while (t<10000)
{

//Vorticity
Vor_Points ();

//Radial velocity
Zero_A_b ();
Vr_Middle_Points ();
Vr_Boundary ();
P_Bi_CGSTAB ();
Change_x_Vr ();

//Axial Velocity
Vz_Points ();

Update_Vor ();



for (int j=0;j<nz;j++)
{
	m[j] = Integration(j);
}

// Max_m & Avg_m
max_m = abs(m[0]-m[0]);
sum_m = m[0];
for (int j=1;j<nz;j++)
{
	if (abs(m[j]-m[0]) > max_m)
	{
		max_m = abs(m[j]-m[0]);
	}
	sum_m = sum_m + m[j];
}
avg_m = (sum_m/nz)/m[0];
Tol = max_m / m[0];

if (t%100 == 0)
{
	cout<<"Iteration ="<<t<<endl;
	cout<<"Error = "<<Tol<<endl;
	Output("Avg");
	Output("Max");
}

t=t+1;

//cout <<"m1 = "<<m[0]<<endl;
//cout<<"m_max = "<<max_m<<endl;
}

Output("Vor");
Output("Vr");
Output("Vz");


return 0;

}

///////////////End_Main code//////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////
//===========================FUNCTIONS===========================================
/////////////////////////////////////////////////////////////////////////////////

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
void Initialize (int IC)
{
	if (IC == 1)
	{
	
	for (int i=0;i<nr;i++)
	{
		for (int j=0;j<nz;j++)
		{
			Vor_1[i][j] = 0.0;
			Vr_1[i][j] = 0.0;
			Vz_1[i][j] = 0.0 ;
		}
	}
	
	}
	if (IC == 2)
	{
	
	for (int i=0;i<nr;i++)
	{
		for (int j=0;j<nz;j++)
		{
			Vor_1[i][j] = 2*clvz*R[i]/(lr*lr);
			Vr_1[i][j] = 0.0;
			Vz_1[i][j] = clvz*(1 - (R[i]/lr)*(R[i]/lr)) ;
		}
	}
	
	}
	
}
///////////////End_Initialize//////////////////////////////////////////////////////

////////////////Begin_Output_Function///////////////////////////////////////////////////
void Output (string label)
{
	if (label == "R")
	{

		cout<<endl<<"R ="<<endl;
		for (int i=0;i<nr;i++)
		{
			cout<<R[i]<<"	";
		}
		
	}
//////////////////////////////////////////	
	if (label == "Z")
	{

		cout<<endl<<"Z ="<<endl;
		for (int j=0;j<nz;j++)
		{
			cout<<Z[j]<<"	";
		}
	}
//////////////////////////////////////////
	if (label == "Vor")
	{
	
		cout<<endl<<"Vorticity ="<<endl;
		for (int i=0;i<nr;i++)
		{
			for (int j=0;j<nz;j++)
			{
				cout<<Vorticity[i][j]<<"	";
				output_file_1<<Vorticity[i][j]<<"	";
			}
			cout<<endl;
			output_file_1<<endl;
		}
	}
	
//////////////////////////////////////////	
	if (label == "A")
	{
			
		cout<<endl<<"A ="<<endl;
		for (int i=0;i<N;i++)
		{
			for (int j=0;j<5;j++)
			{
				cout<<A[i][j]<<"	";
			}
			cout<<endl;
		}
	}
	
//////////////////////////////////////////	
	if (label == "b")
	{
	
		cout<<endl<<"b ="<<endl;
		for (int i=0;i<N;i++)
		{
			cout<<b[i];
			cout<<endl;
		}
	}
//////////////////////////////////////////	
	if (label == "Vr")
	{
	
	cout<<endl<<"Radial Velocity="<<endl;
	for (int i=0;i<nr;i++)
		{
			for (int j=0;j<nz;j++)
			{
				cout<<Vr_1[i][j]<<"	";
				output_file_2<<Vr_1[i][j]<<"	";
			}
			cout<<endl;
			output_file_2<<endl;
		}
	}
//////////////////////////////////////////	
	if (label == "Vz")
	{
	cout<<endl<<"Axial Velocity="<<endl;
	for (int i=0;i<nr;i++)
		{
			for (int j=0;j<nz;j++)
			{
				cout<<Vz_1[i][j]<<"	";
				output_file_3<<Vz_1[i][j]<<"	";
			}
			cout<<endl;
			output_file_3<<endl;
		}
	}
//////////////////////////////////////////	
	if (label == "Max")
	{
		output_file_4<<max_m<<"	";
		output_file_4<<t<<"	";
		output_file_4<<endl;
		
	}	
//////////////////////////////////////////	
	if (label == "Avg")
	{
		output_file_5<<avg_m<<"	";
		output_file_5<<t<<"	";
		output_file_5<<endl;
	}
}
///////////////End_Output_Function//////////////////////////////////////////////////////

////////////////Begin_Vor_Points///////////////////////////////////////////////////

void Vor_Points ()
{
//Midle_Points
	C1_r = (deltat * vis / d)/(deltar * deltar);
	C1_z = (deltat * vis / d)/(deltaz * deltaz);
	for (int i=1;i<nr-1;i++)
	{
		for (int j=1;j<nz-1;j++)
		{
			Vorticity[i][j] = Vor_1[i][j] * (1 - 2*C1_r - 2*C1_z - (deltat * vis / d)/(R[i]*R[i])+ (deltat * Vr_1[i][j]/R[i])) +\
						  Vor_1[i+1][j]	* (C1_r + (deltat * vis)/(d * 2 * deltar*R[i]) - (deltat * Vr_1[i][j])/(2*deltar)) +\
						  Vor_1[i-1][j]	* (C1_r - (deltat * vis)/(d * 2 * deltar*R[i]) + (deltat * Vr_1[i][j])/(2*deltar)) +\
						  Vor_1[i][j+1]	* (C1_z - (deltat * Vz_1[i][j])/(2* deltaz)) +\
						  Vor_1[i][j-1]	* (C1_z + (deltat * Vz_1[i][j])/(2* deltaz)) ;
		}
	}	
	
// axis of symmetry
	int i,j;
	i =1;
	for (j=0;j<nz;j++)
	{
		Vorticity[i][j] = 0.0;
	}

//outflow
	j = nz-1;
	for (int i=0;i<nr;i++)
	{
		Vorticity[i][j] = Vorticity[i][j-1]; 
	}
	
////wall
	i = nr-1;
	for (int j=0;j<nz;j++)
	{
		Vorticity[i][j] = -((Vz_1[i][j] - Vz_1[i-1][j])/(deltar));
	}	

//inflow 
	j = 0;
	for (int i=1;i<nr-1;i++)
	{
		Vorticity[i][j] = ((Vr_1[i][j+1] - Vr_1[i][j])/(deltaz)) - ((Vz_1[i+1][j] - Vz_1[i-1][j])/(2*deltar));
	}
	
}

///////////////End_Vor_Points//////////////////////////////////////////////////////


////////////////Begin_Zero_A_b///////////////////////////////////////////////////
void Zero_A_b ()
{
	for (int i=0;i<N;i++)
	{
		for (int j=0;j<5;j++)
		{
			A[i][j]=0.0;
		}
		b[i]=0.0;
	}
}
///////////////End_Zero_A_b//////////////////////////////////////////////////////

////////////////Begin_Vr_Middle_Points///////////////////////////////////////////////////
void Vr_Middle_Points ()
{
	C1 = 1.0/(deltar * deltar);
	C2 = 1.0/(deltaz * deltaz);

	for (int i=1;i<nr-1;i++)
	{
		for (int j=1;j<nz-1;j++)
		{
			A[i+j*nr][2]= (-2.0*C1 -2.0*C2 -1.0/(R[i]*R[i]));
			A[i+j*nr][3]= (C1 + 1.0/(2.0*R[i]*deltar)) ;
			A[i+j*nr][1]= (C1 - 1.0/(2.0*R[i]*deltar)) ;
			A[i+j*nr][4]= C2 ;
			A[i+j*nr][0]= C2 ;
			
			b[i+j*nr]= (Vorticity[i][j+1] - Vorticity[i][j-1])/(2*deltaz);

		}
		
	}	
	
}
 
///////////////End_Vr_Middle_Points//////////////////////////////////////////////////////


////////////////Begin_Vr_Boundary///////////////////////////////////////////////////

void Vr_Boundary ()
{

// inflow, wall, axis of symmetry
// b is already Zero
for (int j=0;j<nz;j++)
		{
			A[0+j*nr][2]= 1.0 ;
			A[nr-1+j*nr][2]= 1.0;	
		}
for (int i=1;i<nr;i++)
		{
			A[i][2]= 1.0 ;	
		}

// outflow	
for (int i=1;i<nr-1;i++)
	{
		A[(nz-1)*nr+i][2]=1.0 ;
		A[(nz-1)*nr+i][0]=-1.0;
			
	}
}

///////////////End_Vr_Boundary//////////////////////////////////////////////////////

////////////////Begin_Linear_Solver///////////////////////////////////////////////////
// This function solves a tri-diagonal linear system using preconditioned Bi-CGSTAB
void P_Bi_CGSTAB ()
{
	long double r[N],v[N],p[N],r_bar[N],vim1[N],pim1[N],rim1[N],xim1[N];
	long double y[N],K[N],z[N],s[N],t[N],kinvs[N],kinvt[N],esym[N],fsym[N],gsym[N],bsym[N];
	int i,j,k;
	long double errsum,rho,alpha,omega,CGTOL,rhoim1,omegaim1,beta;
	
	// Error limit
	CGTOL = 1e-5;

	// initial guess
	for (i=0;i<N;i++)
	{
		//xim1[i]=b[i]/A[i][2];
		xim1[i] = 0.00001;
	}
	
	//Preconditioner
	for (i=0;i<N;i++)
	{
		K[i] = A[i][2];
	}
	
	//Calculate r_0
	i=0;
	rim1[i]=b[i] - A[i][4]*xim1[i+nr] - A[i][2]*xim1[i] - A[i][3]*xim1[i+1];
	for (i=1;i<nr;i++)
	{
		rim1[i]=b[i] - A[i][4]*xim1[i+nr] - A[i][1]*xim1[i-1] - A[i][2]*xim1[i] - A[i][3]*xim1[i+1];
	}
	for (i=nr;i<N-nr;i++)
	{
		rim1[i]=b[i] - A[i][0]*xim1[i-nr] - A[i][4]*xim1[i+nr] - A[i][1]*xim1[i-1] - A[i][2]*xim1[i] - A[i][3]*xim1[i+1];
	}
	for (i=N-nr;i<N-1;i++)
	{
		rim1[i]=b[i] - A[i][0]*xim1[i-nr] - A[i][1]*xim1[i-1] - A[i][2]*xim1[i] - A[i][3]*xim1[i+1];
	}
	i = N-1;
	rim1[i]=b[i] - A[i][0]*xim1[i-nr] - A[i][1]*xim1[i-1] - A[i][2]*xim1[i];
	
	//Set r_0_bar = r_0
	for (i=0;i<N;i++)
	{
		r_bar[i]=rim1[i];
	}
	
	

// Initialize constants
 
k=0;
rhoim1 = 1.0;
alpha = 1.0;
omegaim1 = 1.0;
for (i=0;i<N;i++)
{
	vim1[i] = 0.0;
    pim1[i] = 0.0;
}


// Begin main loop

errsum=1.0;

while (errsum > CGTOL)
{
	k = k+1;
		
	rho = dot_product(rim1,r_bar);
	beta=rho*alpha/(rhoim1*omegaim1);
	
	for (i=0;i<N;i++)
	{
		p[i]=rim1[i]+beta*(pim1[i]-omegaim1*vim1[i]);		
	}
	
	
	// Solve for y from Ky=p
	
	for (i=0;i<N;i++)
	{
		y[i]=p[i]/K[i];	
	}
	
	// v = Ay (matrix vector multiplication
	
	i=0;
	v[i] = A[i][4]*y[i+nr] + A[i][2]*y[i] + A[i][3]*y[i+1];
	for (i=1;i<nr;i++)
	{
		v[i] = A[i][4]*y[i+nr] + A[i][1]*y[i-1] + A[i][2]*y[i] + A[i][3]*y[i+1];
	}
	for (i=nr;i<N-nr;i++)
	{
		v[i] = A[i][0]*y[i-nr] + A[i][4]*y[i+nr] + A[i][1]*y[i-1] + A[i][2]*y[i] + A[i][3]*y[i+1];
	}
	for (i=N-nr;i<N-1;i++)
	{
		v[i]=A[i][0]*y[i-nr] + A[i][1]*y[i-1] + A[i][2]*y[i] + A[i][3]*y[i+1];
	}
	i = N-1;
	v[i] = A[i][0]*y[i-nr] + A[i][1]*y[i-1] + A[i][2]*y[i];
	
	//update alpha
	alpha=rho/dot_product(r_bar,v);
	
	for (i=0;i<N;i++)
	{
		s[i]=rim1[i]-alpha*v[i];
	}
	
	
	// Solve for z from Kz=s
	
	for (i=0;i<N;i++)
	{
		z[i]=s[i]/K[i];	
	}
	
	// t = Az (matrix vector multiplication
	
	i=0;
	t[i] = A[i][4]*z[i+nr] + A[i][2]*z[i] + A[i][3]*z[i+1];
	for (i=1;i<nr;i++)
	{
		t[i] = A[i][4]*z[i+nr] + A[i][1]*z[i-1] + A[i][2]*z[i] + A[i][3]*z[i+1];
	}
	for (i=nr;i<N-nr;i++)
	{
		t[i] = A[i][0]*z[i-nr] + A[i][4]*z[i+nr] + A[i][1]*z[i-1] + A[i][2]*z[i] + A[i][3]*z[i+1];
	}
	for (i=N-nr;i<N-1;i++)
	{
		t[i]=A[i][0]*z[i-nr] + A[i][1]*z[i-1] + A[i][2]*z[i] + A[i][3]*z[i+1];
	}
	i = N-1;
	t[i] = A[i][0]*z[i-nr] + A[i][1]*z[i-1] + A[i][2]*z[i];
	
	
	for (i=0;i<N;i++)
	{
		kinvs[i]=s[i]/K[i];	
		kinvt[i]=s[i]/K[i];	
		
	}
	
// Update omega
	omega = (dot_product(kinvt,kinvs))/(dot_product(kinvt,kinvt));

// Final answer	
	for (i=0;i<N;i++)
	{
		x[i]=xim1[i]+alpha*y[i]+omega*z[i];
	}
	
	for (i=0;i<N;i++)
	{
		r[i]=s[i]-omega*t[i];
	}
	
	
	// Calculation of Error
	
	errsum = 0.0;
	
	i=0;
	errsum = errsum + pow((b[i] - A[i][2]*x[i] - A[i][3]*x[i+1] - A[i][4]*x[i+nr]),2.0);
	for (i=1;i<nr;i++)
	{
		errsum = errsum + pow((b[i] - A[i][1]*x[i-1] - A[i][2]*x[i] - A[i][3]*x[i+1] - A[i][4]*x[i+nr]),2.0);
	}
	for (i=nr;i<N-nr;i++)
	{
		errsum = errsum + pow((b[i] -  A[i][0]*x[i-nr] - A[i][1]*x[i-1] - A[i][2]*x[i] - A[i][3]*x[i+1] - A[i][4]*x[i+nr]),2.0);
	}
	for (i=N-nr;i<N-1;i++)
	{
		errsum = errsum + pow((b[i] -  A[i][0]*x[i-nr] - A[i][1]*x[i-1] - A[i][2]*x[i] - A[i][3]*x[i+1]),2.0);
	}
	i = N-1;
	errsum = errsum + pow((b[i] -  A[i][0]*x[i-nr] - A[i][1]*x[i-1] - A[i][2]*x[i]),2.0);
	
	errsum = sqrt(errsum)/N;
	
	
// Update all the variables 
	rhoim1 = rho;
	omegaim1 = omega;
	
	for (i=0;i<N;i++)
	{
		rim1[i]=r[i];
		vim1[i]=v[i];
		pim1[i]=p[i];
		xim1[i]=x[i];
	}	
}

}

////////////////End_Linear_Solver///////////////////////////////////////////////////

////////////////Begin_dot_product///////////////////////////////////////////////////
long double dot_product (long double a_1[], long double b_1[] )
{
	long double Dot = 0.0;
	
	for (int i=0;i<N;i++)
	{
		Dot = Dot + a_1[i]*b_1[i];
	}
	
	return Dot;
}

////////////////End_dot_product///////////////////////////////////////////////////

////////////////Begin_Vz_Points///////////////////////////////////////////////////
void Vz_Points()
{
int i,j;	

// Inflow
j = 0;
for (int i=0;i<nr;i++)
{
	Vz_1[i][j] = clvz*(1 - (R[i]/lr)*(R[i]/lr));
}

// Wall
i = nr-1;
for (int j=0;j<nz;j++)
{
	Vz_1[i][j] = 0;	
}	

// Axis of symmetry
i = 0;
for (int j=1;j<nz;j++)
{
	Vz_1[i][j] = Vz_1[i][j-1];
}

// Middle points
	for (int i=1;i<nr-1;i++)
	{
		for (int j=1;j<nz-1;j++)
		{
			Vz_1[i][j] = (-(Vr_1[i][j]/R[i]) - ((Vr_1[i+1][j]-Vr_1[i-1][j])/(2*deltar)))*deltaz + Vz_1[i][j-1];
		}
	}
	
//outflow
j = nz-1;
for (int i=0;i<nr;i++)
{
	Vz_1[i][j] = Vz_1[i][j-1]; 
}

}



///////////////End_Vz_Points//////////////////////////////////////////////////////

////////////////Begin_Change_x_Vr///////////////////////////////////////////////////
void Change_x_Vr ()
{
	int k = 0;
	for (int i=0;i<nr;i++)
	{
		for (int j=0;j<nz;j++)
		{
			Vr_1[i][j] = x[k];
			k = k+1;
		}
	}
}
///////////////End_Change_x_Vr//////////////////////////////////////////////////////


////////////////Begin_Update_Vor///////////////////////////////////////////////////
void Update_Vor ()
{
	for (int i=0;i<nr;i++)
	{
		for (int j=0;j<nz;j++)
		{
			Vor_1[i][j] = Vorticity[i][j];
		}
	}
}
///////////////End_Update_Vor//////////////////////////////////////////////////////

////////////////Begin_Update_Vor///////////////////////////////////////////////////
long double Integration (int j)
{
	long double integration =0.0;
	
	
	integration = (Vz_1[0][j]*R[0]) + (Vz_1[nr-1][j]*R[0]);
	
	for (int i=1;i<nr-1;i++)
	{
		
		if(i%2==0)
  		{
   			 integration = integration + 2 * (Vz_1[i][j]*R[i]);
  		}
  		else
  		{
    		integration = integration + 4 * (Vz_1[i][j]*R[i]);
 		}

	}
	
	integration = 2*3.14*d*integration * deltar/3;
	
	return integration;

}
///////////////End_Update_Vor//////////////////////////////////////////////////////


