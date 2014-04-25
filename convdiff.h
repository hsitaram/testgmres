#include<iostream>
#include<cmath>
#include<fstream>

class convdiff
{
	private:

	static double m_c,m_k;
	static double m_h;
	static double m_phi0,m_phiL;

	public:

	static void setparams(double vel,double diff,double dx,double phi0,double phiL)
	{m_c=vel; m_k=diff; m_h=dx; m_phi0=phi0; m_phiL=phiL;}
	static void findAXcentral(double *,double *,int n);
	static void findAXupwind(double *,double *,int n);
	static void findBcentral(double *b,int n);
	static void findBupwind(double *b,int n);

	static void noprecond(double *MinvX,double *X,int n);
	static void printexactsoln(double xmin,double xmax,double minval,double maxval);

};
