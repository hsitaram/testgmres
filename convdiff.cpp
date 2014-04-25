#include"convdiff.h"

double convdiff::m_c; 
double convdiff::m_k; 
double convdiff::m_h; 
double convdiff::m_phi0; 
double convdiff::m_phiL; 

//=======================================================================
void convdiff::findAXcentral(double *AX,double *X,int n)
{
	double c,k,h;
	double cf1,cf2,cf3;
	double left,right,myself;

	c=m_c; k=m_k; h=m_h;

	cf1 = k+0.5*c*h;
	cf2 = -2.0*k;
	cf3 = k-0.5*c*h;

	for(int i=0;i<n;i++)
	{

		left   = (i==0)?0.0:X[i-1];
		myself = X[i];
		right  = (i==n-1)?0.0:X[i+1];

		AX[i] = cf1*left + cf2*myself + cf3*right;
	}
}
//=======================================================================
void convdiff::findAXupwind(double *AX,double *X,int n)
{
	double c,k,h;
	double cf1,cf2,cf3;
	double left,right,myself;

	c=m_c; k=m_k; h=m_h;

	cf1 = k + c*h;
	cf2 = -2.0*k - c*h;
	cf3 = k;

	for(int i=0;i<n;i++)
	{
		left   = (i==0)?0.0:X[i-1];
		myself = X[i];
		right  = (i==n-1)?0.0:X[i+1];
		
		AX[i] = cf1*left + cf2*myself + cf3*right;
	}
}
//=======================================================================
void convdiff::findBcentral(double *b,int n)
{
	b[0]   =   -(m_k + 0.5*m_c*m_h)*m_phi0;
	b[n-1] =   -(m_k - 0.5*m_c*m_h)*m_phiL;
}
//=======================================================================
void convdiff::findBupwind(double *b,int n)
{
	b[0] = -(m_k + m_c*m_h)*m_phi0;
	b[n-1] = -m_k*m_phiL;
}
//=======================================================================
void convdiff::noprecond(double *MinvX,double *X,int n)
{
	for(int i=0;i<n;i++)
	{
		MinvX[i]=X[i];
	}
}
//=======================================================================
void convdiff::printexactsoln(double xmin,double xmax,double minval,double maxval)
{
	std::ofstream outfile("exactsoln.dat");
	int npoints=100;
	double dx,Pe_per_x,PeL,Pe,x;
	double val;

	dx=(xmax-xmin)/float(npoints-1);

	Pe_per_x = m_c/m_k;
	PeL = (xmax-xmin)*Pe_per_x;

	
	for(int i=0;i<npoints;i++)
	{
		x=xmin+i*dx;
		Pe = x*Pe_per_x;
		
		val=minval+(maxval-minval)*(1.0-exp(Pe))/(1.0-exp(PeL));

		outfile<<x<<"\t"<<val<<"\n";
	}

	outfile.close();
}
