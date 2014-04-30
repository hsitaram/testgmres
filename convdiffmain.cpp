#include"convdiff.h"
#include"gmres.h"

int main()
{
	double *x0,*x,*b,dx;
	int np=101;
	int n,m,it;
	solvergmres kspsolver;
	std::ofstream outfile("soln.dat");

	dx=1.0/float(np-1);

	convdiff::setparams(1,0.1,dx,0.0,1.0);

	n=np-2;
	m=n/2;
	it=n;

	x0 = new double[n]();
	x  = new double[n]();
	b  = new double[n]();

	convdiff::findBupwind(b,n);

	kspsolver.setkspvalues(m,n,it);
	kspsolver.performgmres(b,x0,x,&convdiff::findAXupwind,&convdiff::noprecond);

	convdiff::printexactsoln(0.0,1.0,0.0,1.0);
		
	outfile<<"0.0\t0.0\n";
	for(int i=0;i<n;i++)
	{
		outfile<<(i+1)*dx<<"\t"<<x[i]<<"\n";
	}
	outfile<<"1.0\t1.0";

	return(0);
}
