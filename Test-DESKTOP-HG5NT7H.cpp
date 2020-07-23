#include<iostream>
#include <math.h>
#include<cmath>
#include <fstream>
using namespace std;

/*////////////////Begin_Header///////////////////////////////////////////////////
================================================================================
// This C++ code was written by Behrooz Jadidi_Student ID:501001145
// The first version was generated on MM_DD_YYYY
// This code was written as a solution to AE/ME8112 - Problem_Set3_Question2
// This code solves dT/dt=alpha*(d2T/dx2 + d2T/dy2)+S
// This code Uses a timestep of 0.01 s, a fully implicit scheme, and 80 x 80 equispaced control volumes 
// Also it can solve with unequal dx and dy
================================================================================
/*///////////////End_Header//////////////////////////////////////////////////////

//Global variables
//Geometry

#define N 6
float T[N],TF[N];
float b [N]={1,-1,2,0,3,-4};
float A[N][N]={3,1,0,0,2,0,1,3,1,0,0,2,0,1,3,1,0,0,0,0,1,3,1,0,2,0,0,1,3,1,0,2,0,0,1,3};



//Decleration of Functions
void Zero_A_b ();
void Cout_Ax_b ();
void Cout_T ();
void Linear_Solver ();
void matmul(float rr[N], float AA[N][N], float TT[N],float bb[N]);
float dot_product(float rr_hat[N],float rr[N]);
void matmul_2(float vv[N], float AA[N][N], float yy[N]);

//

////////////////Begin_Main code///////////////////////////////////////////////////
int main()
{	


//Zero_A_b ();
Cout_Ax_b ();

Linear_Solver ();

Cout_T ();

return 0;
}
///////////////End_Main code//////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////
//===========================FUNCTIONS===========================================
/////////////////////////////////////////////////////////////////////////////////

/*////////////////Begin_#///////////////////////////////////////////////////

/*///////////////End_#//////////////////////////////////////////////////////

////////////////Begin_Zero N*N///////////////////////////////////////////////////
void Zero_A_b ()
{
	for (int i=0;i<N;i++)
	{
		for (int j=0;j<N;j++)
		{
			A[i][j]=0.0;
		}
		b[i]=0.0;
		T[i]=0.0;
	}
	
}
///////////////End_Zero N*N//////////////////////////////////////////////////////



////////////////Begin_Cout Ax=b///////////////////////////////////////////////////
void Cout_Ax_b ()
{
	cout<<endl<<"Ax=b"<<endl;
	for (int i=0;i<N;i++)
	{
		for (int j=0;j<N;j++)
		{
			cout<<A[i][j]<<" ";
		}
		cout<<" "<<"x"<<"="<<" "<<b[i];
		cout<<endl;
	}
	
}
///////////////End_Cout Ax=b//////////////////////////////////////////////////////

////////////////Begin_Cout_T///////////////////////////////////////////////////
void Cout_T ()
{
	cout<<endl<<"T="<<endl;
	for (int i=0;i<N;i++)
	{
		cout<<T[i];
		cout<<endl;
	}
	
}
///////////////End_Cout_T//////////////////////////////////////////////////////


////////////////Begin_Linear_Solver///////////////////////////////////////////////////

void Linear_Solver ()
{
	//declerimg parameters
	float r[N],r_hat[N],y[N],tt[N],s[N],k[N],p[N],v[N],z[N];
	float rho,rho_old,alpha,omega,beta,resid;
	int i,j;
	bool converged ;
	
	// initial guess according to lecture content

	for (int i=0;i<N;i++)
	{
		T[i]=b[i]/A[i][i];
	}
	
	// variable initialize
	matmul (r,A,T,b);
	
	for (int i=0;i<N;i++)
	{
		r_hat[i]=r[i];
	}
	rho=1.0;
	alpha=1.0;
	omega=1.0;
	for (int i=0;i<N;i++)
	{
		v[i]=0.0;
		p[i]=0.0;
	}
	
	// precondition vector K=diag(A)
	for (int i=0;i<N;i++)
	{
		k[i]=A[i][i];
	} 
	
	// Preconditioned BICGSTAB algorithm main body
	
	converged = false;
	
	while (converged == false)// check if the norm is satisfied 
	{
	
	rho_old=rho;
	
	rho=dot_product(r_hat,r);
	beta= (rho/rho_old) * (alpha/omega);
	
	for (int i=0;i<N;i++)
	{
		p[i]=r[i] + beta * (p[i] - omega*v[i]);
		y[i]=p[i]/k[i];
	}  
  
  	matmul_2(v,A,y);
  	alpha= rho/dot_product(r_hat,v);
  	
  	for (int i=0;i<N;i++)
	{
		s[i]=r[i]-alpha*v[i];
		z[i]=s[i]/k[i];
	}  
	
  	matmul_2(tt,A,z);
  	omega= dot_product(tt,s)/dot_product(tt,tt);
  	
  	for (int i=0;i<N;i++)
	{
		T[i]=T[i] + alpha *y[i] + omega * z[i];
		r[i]=s[i] - omega * tt[i];
	} 
 	
 	resid = 0.0;
  
 	for (int i=0;i<N;i++)
 	{
 		resid=resid + r[i]*r[i];
	}
  	resid = sqrt(resid)/N;
	if( resid < 1e-9)
     converged= true;

  
	}


}
////////////////End_Linear_Solver///////////////////////////////////////////////////

////////////////Begin_matmul///////////////////////////////////////////////////

void matmul(float rr[N], float AA[N][N], float TT[N],float bb[N])
{
	float sum[N];
	for (int i=0;i<N;i++)
	{
		sum[i]=0.0;
		for (int j=0;j<N;j++)
		{
		
		sum[i]=sum[i]+AA[i][j]*TT[j];
		
		}
		rr[i]=b[i]-sum[i];
	}
}

///////////////End_matmul//////////////////////////////////////////////////////

////////////////Begin_matmul_2///////////////////////////////////////////////////

void matmul_2(float vv[N], float AA[N][N], float yy[N])
{
	float sum[N];
	for (int i=0;i<N;i++)
	{
		sum[i]=0.0;
		for (int j=0;j<N;j++)
		{
		
		sum[i]=sum[i]+AA[i][j]*yy[j];
		
		}
		vv[i]=sum[i];
	}
}

///////////////End_matmul_2//////////////////////////////////////////////////////

////////////////Begin_dot_product///////////////////////////////////////////////////

float dot_product(float rr_hat[N],float rr[N])
{
	float sum=0.0;
	for (int i=0;i<N;i++)
		{
		
		sum=sum+rr_hat[i]*rr[i];
		
		}
	
	return sum;
}

///////////////End_dot_product//////////////////////////////////////////////////////
