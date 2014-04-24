#include"gmres.h"

//Private functions
void solvergmres::m_printvec(std::string str,double *v,int n)
{
	std::cout<<"\n"<<str<<"\t";
	for(int i=0;i<n;i++)
	{
		std::cout<<v[i]<<"\t";
	}
	std::cout<<"\n";

}
//========================================================================
void solvergmres::m_printmat(std::string str,double *mat,int m,int n)
{
	std::cout<<"\n"<<str<<"\n";

	for(int i=0;i<m;i++)
	{
		for(int j=0;j<n;j++)
		{
			std::cout<<mat[i*n+j]<<"\t";
		}
		std::cout<<"\n";
	}
	std::cout<<"\n";

}
//========================================================================
double solvergmres::m_findnorm(double *v1,int n)
{
	double norm=0;

	for(int i=0;i<n;i++)
	{
		norm=norm+v1[i]*v1[i];
	}

	return(sqrt(norm));
}
//========================================================================
double solvergmres::m_innerproduct(double *v1,double *v2,int n)
{
	double innerprod=0;

	for(int i=0;i<n;i++)
	{
		innerprod=innerprod+v1[i]*v2[i];
	}

	return(innerprod);
}
//========================================================================
void solvergmres::m_getkspvector(double *v1,int vecnum)
{
	int m,n;
	n=m_matsize;

	for(int i=0;i<n;i++)
	{
		v1[i]=m_kspvectors[vecnum*n+i];
	}
}
//========================================================================
void solvergmres::m_setkspvector(double *vec,int vecnum)
{
	int n;
	n=m_matsize;

	for(int i=0;i<n;i++)
	{
		m_kspvectors[vecnum*n+i]=vec[i];
	}
}
//========================================================================
void solvergmres::m_addvectors(double *v1,double *v2,double *v12,int n,double a,double b)
{
	for(int i=0;i<n;i++)
	{
		v12[i]=a*v1[i]+b*v2[i];
	}
}
//========================================================================
void solvergmres::m_copyvector(double *v1,double *v2,int n) //(dest,source,size)
{
	for(int i=0;i<n;i++)
	{
		v1[i]=v2[i];
	}
}
//========================================================================
bool solvergmres::m_arnoldialgorithm(double *v1, void (* findAX)(double *,double *,int ),void (*precond)(double *,double *,int))
{
	int m,i,j,index,n;
	double *Avj,*vj,*vi;
	double *wj,*tempvec;
	double *MinvAvj;

	bool lucky; //when norm becomes 0, KSP 
	//is no longer linearly independent.
	//we would have got the best solution.

	m = m_kspdim;
        n = m_matsize;

	Avj     = new double[n]();
	MinvAvj = new double[n]();
	vj      = new double[n]();
	vi      = new double[n]();
	wj      = new double[n]();
	tempvec = new double[n]();

	m_copyvector(vj,v1,n); 
	m_setkspvector(vj,0);

	lucky=false;
	for(j=0;j<m;j++)
	{
		m_getkspvector(vj,j);
		findAX(Avj,vj,n);
		precond(MinvAvj,Avj,n);

		m_copyvector(Avj,MinvAvj,n);
		//Avj is now M^-1 A vj
		//remember we are solving M^-1 A X = M^-1 b
		
		for(i=0;i<=j;i++)
		{
			index=i*m+j;
		        m_getkspvector(vi,i);
			m_Hessbergmat[index]=m_innerproduct(Avj,vi,n);
		}

		m_copyvector(wj,Avj,n);

		for(i=0;i<=j;i++)
		{
			index=i*m+j;
			m_getkspvector(vi,i);
			m_addvectors(wj,vi,tempvec,n,1.0,-m_Hessbergmat[index]);
			m_copyvector(wj,tempvec,n);
		}
		

		m_Hessbergmat[(j+1)*m+j] = m_findnorm(wj,n);

		if(m_Hessbergmat[(j+1)*m+j] > 0.0)
		{
			for(int i=0;i<n;i++)
			{
				wj[i]=wj[i]/m_Hessbergmat[(j+1)*m+j];
			}
		}
		else
		{
			lucky=true;
			break;
		}		

		m_setkspvector(wj,j+1);

	}	

	return(lucky);
}
//========================================================================
void solvergmres::m_leastsqminimize(double *y,double beta)
{
	int m,index;
	double *beta_e1;
	double c,s,h_up,h_down,dtr;
	double val1,val2;

	m=m_kspdim;
	beta_e1 = new double[m+1]();

	beta_e1[0] = beta;

	//convert H into QR
	for(int i=0;i<m;i++)
	{
		h_up   = m_Hessbergmat[i*m + i    ];
		h_down = m_Hessbergmat[(i+1)*m + i];

		dtr = sqrt(h_up*h_up + h_down*h_down);
		
		c=h_up/dtr; s=h_down/dtr;

		for(int j=0;j<m;j++)
		{
			h_up   = m_Hessbergmat[i*m+j];
			h_down = m_Hessbergmat[(i+1)*m+j];

			//perform rotations
			//ith row
			m_Hessbergmat[i*m + j    ] =  c*h_up+s*h_down;
			//(i+1)th row
			m_Hessbergmat[(i+1)*m + j] = -s*h_up+c*h_down;

		}

		val1 =  c*beta_e1[i] + s*beta_e1[i+1];
		val2 = -s*beta_e1[i] + c*beta_e1[i+1];

		beta_e1[i]=val1; beta_e1[i+1]=val2;
	
	}


	// ||Hm y - beta e1|| = || QR y - Q Q^T beta e1||
	// || Q ( Ry - Q^T beta e1) || = || Ry - Q^T beta e1||

	//solve least squares problem
	y[m-1] = beta_e1[m-1]/m_Hessbergmat[(m-1)*m+(m-1)];

	for(int i=m-2;i>=0;i--)
	{
		y[i]=beta_e1[i];

		for(int j=i+1;j<m;j++)
		{
			y[i]=y[i]-m_Hessbergmat[i*m+j]*y[i+1];
		}

		y[i]=y[i]/m_Hessbergmat[i*m+i];
	}

}
//========================================================================
//Public functions
//========================================================================
void solvergmres::setkspvalues(int m,int n,int it)
{
	m_kspdim = m;       //dimension of Krylov subspace
	m_numit	   = it;     //number of restart iterations
	m_matsize  = n;

	//each ksp vector has dimension n and there are m+1 of them
	//the vectors themselves are rows here, not columns

	// v1.....
	// v2.....
	// .
	// .
	// vm+1....
	

	m_kspvectors =  new double[n*(m+1)]();

	//m+1 rows with m columns
	m_Hessbergmat = new double[(m+1)*m]();

}
//========================================================================
void solvergmres::performgmres(double *b,double *x0,double *x,
		void (*findAX)(double *,double *,int),void (*precond)(double *,double *,int ))
{
	int n,m;
	double beta;
	double *v1;
	double *tempvec;
	double *v,*y;
	bool arnoldistopped;
	double *r,*r0,*Ax0,*Ax;
	double *Minvr;
	
	n = m_matsize;
	m = m_kspdim;
	
	v1      = new double[n]();
	tempvec = new double[n]();
	v       = new double[n]();
	r	= new double[n]();
	r0	= new double[n]();
	Ax0	= new double[n]();
	Ax	= new double[n]();
	Minvr   = new double[n]();
	
	y       = new double[m]();

	//finding r0
	findAX(Ax0,x0,n);
	m_addvectors(b,Ax0,r0,n,1.0,-1.0);
	precond(Minvr,r0,n);
	m_copyvector(r0,Minvr,n);

	//initial residual is r0=M^-1(b-Ax0)
	//we are solving M^-1 A X = b


	m_copyvector(r,r0,n);
	m_copyvector(x,x0,n);

	for(int it=0;it<m_numit;it++)
	{
		std::cout<<"restart iteration:"<<it<<"\t";
		beta = m_findnorm(r,n);

		for(int i=0;i<n;i++)
		{
			v1[i]=r[i]/beta;
		}

		arnoldistopped = m_arnoldialgorithm(v1,findAX,precond);

		if(arnoldistopped)
		{
			std::cout<<"lucky condition\n";
			break;
		}

		m_leastsqminimize(y,beta);

		for(int i=0;i<m;i++)
		{
			m_getkspvector(v,i);
			m_addvectors(x,v,tempvec,n,1.0,y[i]);
			m_copyvector(x,tempvec,n);
		}	

		//finding new residual
		findAX(Ax,x,n);
		m_addvectors(b,Ax,r,n,1.0,-1.0);
		precond(Minvr,r,n);
		m_copyvector(r,Minvr,n);

		//m_printvec("x",x,n);
		std::cout<<"norm of residual:"<<m_findnorm(r,n)<<"\n";
		//std::cout<<"******************\n";
	}

}
//========================================================================
