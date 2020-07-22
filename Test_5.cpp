#include<iostream>
#include <math.h>
#include<cmath>
#include <fstream>
using namespace std;

/*////////////////Begin_Header///////////////////////////////////////////////////
================================================================================
// This C++ code was written by Behrooz Jadidi_Student ID:501001145
// The first version was generated on 07_22_2020
// This code was written as a test fo Bi-CGSTAB
// This code use preconditioned Bi-CGSTAB to solve Ax = b 

================================================================================
/*///////////////End_Header//////////////////////////////////////////////////////

#define nr 3
#define nz 3
#define lr 1.0
#define lz 10.0
#define N nr*nz

long double b [N]={17,33,28,49,85,63,39,75,45};
long double A[N][5]={0,0,1,2,3,\
			         0,2,1,3,4,\
			         0,2,1,0,3,\
			         2,0,1,3,4,\
			         2,3,1,4,5,\
			         2,3,1,0,4,\
			         2,0,1,3,0,\
			         2,3,1,4,0,\
			         2,3,1,0,0};
long double x [N]={0,0,0,0,0,0,0,0,0};


//Decleration of Functions
void Cout_Ax_b (long double A[][5],long double b[]);
void P_Bi_CGSTAB ();
long double dot_product (long double a_1[], long double b_1[] );

//

////////////////Begin_Main code///////////////////////////////////////////////////
int main()
{	
//int N,nr,nz;
//nr = 3;
//nz =3;
//N = nr*nz;





Cout_Ax_b (A,b);
P_Bi_CGSTAB ();

cout<<"x = "<<endl;
for (int i=0;i<N;i++)
{
	cout<<x[i]<<endl;
}


return 0;
}
///////////////End_Main code//////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////
//===========================FUNCTIONS===========================================
/////////////////////////////////////////////////////////////////////////////////

/*////////////////Begin_#///////////////////////////////////////////////////

/*///////////////End_#//////////////////////////////////////////////////////




////////////////Begin_Cout Ax=b///////////////////////////////////////////////////
void Cout_Ax_b (long double A[][5],long double b[])
{
	cout<<endl<<"Ax=b"<<endl;
	for (int i=0;i<N;i++)
	{
		for (int j=0;j<5;j++)
		{
			cout<<A[i][j]<<" ";
		}
		cout<<" "<<"x"<<"="<<" "<<b[i];
		cout<<endl;
	}
	
}
///////////////End_Cout Ax=b//////////////////////////////////////////////////////


////////////////Begin_Linear_Solver///////////////////////////////////////////////////
// This function solves a tri-diagonal linear system using preconditioned Bi-CGSTAB
void P_Bi_CGSTAB ()
{
	long double r[N],v[N],p[N],r_bar[N],vim1[N],pim1[N],rim1[N],xim1[N];
	long double y[N],K[N],z[N],s[N],t[N],kinvs[N],kinvt[N],esym[N],fsym[N],gsym[N],bsym[N];
	int i,j,k;
	long double errsum,rho,alpha,omega,CGTOL,rhoim1,omegaim1,beta;
	
	
	CGTOL = 1e-9;

	// initial guess
	for (i=0;i<N;i++)
	{
		xim1[i]=b[i]/A[i][2];
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
	
	//Track operating of solver
	if (k%1000 == 0)
	{
		cout<<"BiCGSTAB solver iteration ="<<k<<endl<<"Residual ="<<errsum<<endl;
	}
	
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



