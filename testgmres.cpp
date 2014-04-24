#include"gmres.h"

//=========================================================
void findAX(double *AX,double *X,int n)
{
	for(int i=0;i<n;i++)
	{
		if(i == 0)
		{
			AX[i] = -2*X[i]+X[i+1];
		}
		else
		if( i == n-1)
		{
			AX[i]=-2*X[i]+X[i-1];
		}
		else
		{
			AX[i]=X[i+1]+X[i-1]-2*X[i];
		}
	}
}
//=========================================================
void noprecond(double *MinvX,double *X,int n)
{
	for(int i=0;i<n;i++)
	{
		MinvX[i]=X[i];
	}
}
//=========================================================
void gauss_seidel(double *MinvX,double *X,int n)
{
	double cf1,cf2,cf3;
	double left,right;
	int numit;

	cf1=1.0; cf2=-2.0; cf3=1.0;
	numit=500;

	for(int it=0;it<numit;it++)
	{
		for(int i=0;i<n;i++)
		{
			left   = (i==0)?0.0:MinvX[i-1];
			right  = (i==n-1)?0.0:MinvX[i+1];

			MinvX[i] = (X[i] -cf1*left -cf3*right)/cf2;
		}
	}
}
//=========================================================
void findbvec(double *b,int n)
{
	b[0]   = 0.0;
	b[n-1] = -1.0;
}
//=========================================================
int main(int argc,char *argv[])
{
	double *x0,*x,*b;
	int m,n,it;
	double *testvec;
	std::string pctype;
	std::ofstream outfile("soln.dat");

	if(argc == 2)
	{
		pctype=argv[1];
	}
	else
	{
		pctype="none";
	}

	m=40; n=50; it=20;

	x0 = new double[n]();
	b  = new double[n]();
	x  = new double[n]();

	solvergmres obj;

	findbvec(b,n);

	obj.setkspvalues(m,n,it);

	if(pctype == "gs")
	{
		obj.performgmres(b,x0,x,&findAX,&gauss_seidel);
	}
	else
	{
	 obj.performgmres(b,x0,x,&findAX,&noprecond);
	}

	//gauss_seidel(x,b,n);

	for(int i=0;i<n;i++)
	{
		outfile<<i<<"\t"<<x[i]<<"\n";
	}

	return(0);
}
//=========================================================
