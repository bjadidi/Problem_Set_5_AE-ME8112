#include<iostream>
#include <math.h>
#include<cmath>
#include <fstream>
using namespace std;
/*////////////////Begin_Header///////////////////////////////////////////////////
================================================================================
// This C++ code was written by Behrooz Jadidi_Student ID:501001145
// The first version was generated on 06_20_2020
// This code was written as a solution to AE/ME8112 - Problem_Set3_Question1_Part a and b
// This code solves 1D convection/diffusion problem governed by: d/dx(puf)=d/dx(Ldf/dx)
// In the domain 0<=x<=l , f(0)=1 and f(l)=0
================================================================================
/*///////////////End_Header//////////////////////////////////////////////////////

//Global variables
//Geometry
#define N 22
#define l 1.0
long double A[N][N],b[N],T[N],T_exact[N],x[N];
long double deltax,alpha,Pe,L_2;
int time_step;
fstream output_file ("T_Question1.xls",ios::out);



//Decleration of Functions
void Cout_Zero_A_T_b ();
void Interior_Points ();
void West_Wall_Points ();
void East_Wall_Points ();
void Cout_Ax_b ();
void Linear_Solver ();
void Cout_T ();

void T_exact_Calculation ();
void Error_Calculation ();



//

////////////////Begin_Main code///////////////////////////////////////////////////
int main()
{	
//Variables_Decleration
deltax=l/N;
cout<<"deltax="<<deltax;
Pe=50.0;
//alpha=pow(Pe,2)/(pow(Pe,2)+5);
alpha=0.0;
//alpha=1.0;
x[0] = deltax/2;
for (int i=0;i<N-1;i++)
{
	x[i+1]=x[i]+deltax;
}

Cout_Zero_A_T_b ();
Interior_Points ();
West_Wall_Points ();
East_Wall_Points ();
Cout_Ax_b ();
Linear_Solver ();
//Cout_T ();
T_exact_Calculation ();
Error_Calculation ();
Cout_T ();

//cout<<endl<<"L_2="<<L_2;

return 0;
}
///////////////End_Main code//////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////
//===========================FUNCTIONS===========================================
/////////////////////////////////////////////////////////////////////////////////

/*////////////////Begin_#///////////////////////////////////////////////////

/*///////////////End_#//////////////////////////////////////////////////////

////////////////Begin_Zero_A_T_b///////////////////////////////////////////////////
void Cout_Zero_A_T_b ()
{
	for (int i=0;i<N;i++)
	{
		for (int j=0;j<N;j++)
		{
			A[i][j] = 0.0;
		}
		b[i] = 0.0;
		T[i] = 0.0;
	}
	
}
///////////////End_Zero_A_T_b//////////////////////////////////////////////////////

////////////////Begin_Interior_Points///////////////////////////////////////////////////
void Interior_Points ()
{
	for (int i=1;i<N-1;i++)
		{
			A[i][i]= 2.0/deltax - (Pe/2.0)*(1.0-alpha) + (Pe/2.0)*(1.0+alpha);
			A[i][i-1]= -1.0/deltax - (Pe/2.0)*(1.0+alpha);
			A[i][i+1]= -1.0/deltax + (Pe/2.0)*(1.0-alpha);
			
			b[i]=0.0;
		}
}
///////////////End_Interior_Points//////////////////////////////////////////////////////

////////////////Begin_West_Wall_Points///////////////////////////////////////////////////
void West_Wall_Points ()
{
	int i=0;
			A[i][i]= 3.0/deltax - (Pe/2.0)*(1.0-alpha) + (Pe)*(1.0+alpha);
			
			A[i][i+1]= -1.0/deltax + (Pe/2.0)*(1.0-alpha);
			
			b[i]=2.0/deltax + Pe*(1.0+alpha);
		
}
///////////////End_West_Wall_Points//////////////////////////////////////////////////////

////////////////Begin_East_Wall_Points///////////////////////////////////////////////////
void East_Wall_Points ()
{
	int i=N-1;
			A[i][i]= 3.0/deltax - (Pe)*(1.0-alpha) + (Pe/2.0)*(1.0+alpha);
			A[i][i-1]= -1.0/deltax - (Pe/2.0)*(1.0+alpha);
			
			b[i]=0.0;

}
///////////////End_East_Wall_Points//////////////////////////////////////////////////////

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
		cout<<" "<<"@"<<"="<<" "<<b[i];
		cout<<endl;
	}
	
}
///////////////End_Cout Ax=b//////////////////////////////////////////////////////

////////////////Begin_Linear_Solver///////////////////////////////////////////////////
void Linear_Solver ()
{
//Decomposition 
	for (int i=1;i<N;i++)
	{
		A[i][i-1] = A[i][i-1]/A[i-1][i-1];
		A[i][i] = A[i][i] - A[i][i-1]*A[i-1][i];
	}

//Forward Substitution//
	for (int i=1;i<N;i++)
	{
		b[i] = b[i] - A[i][i-1]*b[i-1];
	}

//Backward Substitution//
	T[N-1]=b[N-1]/A[N-1][N-1];
	for (int i=N-2;i>-1;i--)
	{
		T[i] = (b[i] - A[i][i+1]*T[i+1])/A[i][i];
	}

}
///////////////End_Linear_Solver//////////////////////////////////////////////////////

////////////////Begin_Cout_T///////////////////////////////////////////////////
void Cout_T ()
{
//	cout<<endl<<"T="<<endl;
//	output_file<<"T"<<"	"<<"x"<<endl;
//	for (int i=0;i<N;i++)
//	{
//		cout<<"T"<<i+1<<"="<<T[i]<<endl;
//		output_file<<T[i]<<"	"<<x[i]<<endl;
//
//	}
//	
	cout<<endl<<"For N="<<N<<endl;	
cout<<endl<<"L_2="<<L_2<<endl<<endl;

//cout<<"x[i]"<<"		"<<"T[i]"<<"		"<<"T_exact[i]"<<"		"<<"Diff"<<endl;
cout<<"Diff(T_exact - T)"<<endl;
//cout<<"-----"<<"		"<<"-----"<<"		"<<"-----"<<"			"<<"-----"<<endl;
cout<<"-----------"<<endl;;
for (int i=0;i<N;i++)
{
	//cout<<x[i]<<"		"<<T[i]<<"		"<<T_exact[i]<<"		"<<(-T[i]+T_exact[i])<<endl;
	cout<<(-T[i]+T_exact[i])<<endl;
}
	
}
///////////////End_Cout_T//////////////////////////////////////////////////////

////////////////Begin_T_exact_Calculation///////////////////////////////////////////////////
void T_exact_Calculation ()
{
	for (int i=0;i<N;i++)
	{
		T_exact[i] = 1.0 - ((exp(Pe*x[i]) - 1.0) / (exp(Pe) - 1.0));
	}
	
}
///////////////End_T_exact_Calculation//////////////////////////////////////////////////////

////////////////Begin_Error_Calculation///////////////////////////////////////////////////
void Error_Calculation ()
{
	L_2 = 0.0;
	for (int i=0;i<N;i++)
	{
		L_2 = L_2 + (T[i] - T_exact[i])*(T[i] - T_exact[i]);
	}
	L_2 = sqrt(L_2)/N;
}
///////////////End_Error_Calculation//////////////////////////////////////////////////////
