#include<iostream>
#include<cmath>
#include<string>
#include<fstream>

class solvergmres
{
	private:

		int m_kspdim;
		int m_numit;
		int m_matsize;

		double *m_kspvectors;
		double *m_Hessbergmat;

		double m_innerproduct(double *v1,double *v2,int n);
		double m_findnorm(double *v1,int n);
		
		void   m_getkspvector(double *v1,int vecnum);
		void   m_setkspvector(double *vec,int vecnum);
		void   m_addvectors(double *v1,double *v2,double *v12,int n,double a,double b);
		void   m_copyvector(double *v1,double *v2,int n); //(dest,source,size)
		bool   m_arnoldialgorithm(double *v1,void (*findAX)(double *,double *,int),
				void (*precond)(double *,double *,int));
		void   m_leastsqminimize(double *y,double beta);

		void   m_printvec(std::string str,double *v,int n);
		void   m_printmat(std::string str,double *mat,int m,int n);

	public:
		void setkspvalues(int m,int n,int it);
		void performgmres(double *b,double *x0,double *x,
		void (* findAX)(double *,double *,int ),void(*precond)(double *,double *,int));  

};
