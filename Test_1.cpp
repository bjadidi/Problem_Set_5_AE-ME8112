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
#define nr 40
#define nz 40
#define lr 0.01
#define lz 0.1
#define N nr*nz
long double Vor_1[nr][nz],Vr_1[nr][nz],Vz_1[nr][nz],Vor_2[nr][nz],Vr_2[nr][nz],Vz_2[nr][nz],R[nr],Z[nz],Velocity[N],v[N],r[N],tt[N],z[N],y[N];
float T[N],TF[N],A[N][N],b[N];
long double deltar,deltaz,deltat,Vor_0,Vr_0,Vz_0,d,vis,C1_r,C1_z,C1,C2;
int i,j,t;
fstream output_file ("PS5.xls",ios::out);



//Decleration of Functions
void Mesh_Generation ();
void Initialize ();
void Output_Function ();
void Vor_Middle_Points ();
void Vor_Boundary ();
void Vr_Middle_Points ();
void Zero_A_b ();
void Vr_Boundary ();
void Linear_Solver ();
void matmul(float rr[N], float AA[N][N], float TT[N],float bb[N]);
float dot_product(float rr_hat[N],float rr[N]);
void matmul_2(float vv[N], float AA[N][N], float yy[N]);







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
Vor_Boundary ();

//Radial velocity
Zero_A_b ();
Vr_Middle_Points ();
Vr_Boundary ();

Linear_Solver ();
//Change_Vr ();


//Radial Velocity



//Axial Velocity



//Output_Function
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

////////////////Begin_Zero_A_b///////////////////////////////////////////////////
void Zero_A_b ()
{
	for (int i=0;i<N;i++)
	{
		for (int j=0;j<N;j++)
		{
			A[i][j]=0.0;
		}
		b[i]=0.0;
	}
}
///////////////End_Zero_A_b//////////////////////////////////////////////////////

////////////////Begin_Vor_Middle_Points///////////////////////////////////////////////////

void Vor_Middle_Points ()
{
	C1_r = (deltat * vis / d)/(deltar * deltar);
	C1_z = (deltat * vis / d)/(deltaz * deltaz);
	for (int i=1;i<nr-1;i++)
	{
		for (int j=1;j<nz-1;j++)
		{
			Vor_2[i][j] = Vor_1[i][j] * (1 - 2*C1_r - 2*C1_z - (deltat * vis / d)/(R[i]*R[i])+ (deltat * Vr_1[i][j]/R[i])) +\
						  Vor_1[i+1][j]	* (C1_r + (deltat * vis)/(d * 2 * deltar*R[i]) - (deltat * Vr_1[i][j])/(2*deltar)) +\
						  Vor_1[i-1][j]	* (C1_r - (deltat * vis)/(d * 2 * deltar*R[i]) + (deltat * Vr_1[i][j])/(2*deltar)) +\
						  Vor_1[i][j+1]	* (C1_z - (deltat * Vz_1[i][j])/(2* deltaz)) +\
						  Vor_1[i][j-1]	* (C1_z + (deltat * Vz_1[i][j])/(2* deltaz)) ;
		}
	}	
	
	
}

///////////////End_Vor_Middle_Points//////////////////////////////////////////////////////

////////////////Begin_Vor_Boundary///////////////////////////////////////////////////

void Vor_Boundary ()
{
// axis of symmetry
for (int j=0;j<nz;j++)
{
	Vor_2[0][j] = 0.0;
}

// inflow
for (int i=0;i<nr;i++)
{
	Vor_2[i][0] = ((Vr_1[i][j+1] - Vr_1[i][j-1])/(2*deltaz)) - ((Vz_1[i+1][j] - Vz_1[i-1][j])/(2*deltar));
}

//wall
for (int j=0;j<nz;j++)
{
	Vor_2[nr][j] = - ((Vz_1[i+1][j] - Vz_1[i-1][j])/(2*deltar));	
}	

//outflow
for (int i=0;i<nr;i++)
{
	Vor_2[i][nz] = Vor_2[i][nz-1]; 
}
	
}

///////////////End_Vor_Boundary//////////////////////////////////////////////////////


////////////////Begin_Output_Function///////////////////////////////////////////////////
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
	
	cout<<endl<<"Vorticity in t ="<<t<<endl;
	for (int i=0;i<nr;i++)
	{
		for (int j=0;j<nz;j++)
		{
			cout<<Vor_2[i][j]<<"	";
		}
		cout<<endl;
	}
	
	
	cout<<endl<<"A ="<<t<<endl;
	for (int i=0;i<N;i++)
	{
		for (int j=0;j<N;j++)
		{
			cout<<A[i][j]<<"	";
		}
		cout<<endl;
	}
	
	cout<<endl<<"T="<<endl;
	for (int i=0;i<N;i++)
	{
		cout<<T[i];
		cout<<endl;
	}
}
///////////////End_Output_Function//////////////////////////////////////////////////////

////////////////Begin_Vr_Middle_Points///////////////////////////////////////////////////
void Vr_Middle_Points ()
{
	C1 = 1/(deltar * deltar);
	C2 = 1/(deltaz * deltaz);
	for (int i=1;i<nr-1;i++)
	{

		for (int j=1;j<nz-1;j++)
		{
			A[i+j*nr][i+j*nr]= (-2*C1 -2*C2 -1/R[i]) ;
			A[i+j*nr][i+j*nr+1]= (C1 + 1/(2*R[i]*deltar)) ;
			A[i+j*nr][i+j*nr-1]= (C1 - 1/(2*R[i]*deltar)) ;
			A[i+j*nr][i+j*nr+nr]= C2 ;
			A[i+j*nr][i+j*nr-nr]= C2 ;
			
			b[i+j*nr]= (Vor_2[i][j+1] - Vor_2[i][j-1])/(2*deltaz);
		}
	}	
}
 
///////////////End_Vr_Middle_Points//////////////////////////////////////////////////////


////////////////Begin_Vr_Boundary///////////////////////////////////////////////////

void Vr_Boundary ()
{

// inflow, wall, axis of symmetry
// it is already Zero


// outflow	
for (int i=1;i<nr-1;i++)
	{
		A[(nz-1)*nr+i][(nz-1)*nr+i]=1.0 ;
		A[(nz-1)*nr+i][(nz-1)*nr+i-nr]=-1.0;
			
	}
}

///////////////End_Vr_Boundary//////////////////////////////////////////////////////

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
		//T[i]=b[i]/A[i][i];
		T[i]=0.00001;
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
  	alpha= rho/dot_product( r_hat,v);
  	
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
	if( resid < 1e-3)
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
		
		sum[i]=sum[i]+A[i][j]*T[j];
		
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






