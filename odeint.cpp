void	rkck(double *y, double *dydx, int n, const double x,
	      const double h, double *yout, double *yerr,
	      void derivs(const double, double *, double *))
{
	static const double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
	b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42 = -0.9,
	b43=1.2, b51 = -11.0/54.0, b52=2.5, b53 = -70.0/27.0,
	b54=35.0/27.0, b61=1631.0/55296.0, b62=175.0/512.0,
	b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0,
	c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0,
	dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0,
	dc4=c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6=c6-0.25;
	int i;
	
	double ak2[n],ak3[n],ak4[n],ak5[n],ak6[n],ytemp[n];
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+b21*h*dydx[i];
	derivs(x+a2*h,ytemp,ak2);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	derivs(x+a3*h,ytemp,ak3);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	derivs(x+a4*h,ytemp,ak4);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	derivs(x+a5*h,ytemp,ak5);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	derivs(x+a6*h,ytemp,ak6);
	for (i=0;i<n;i++)
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=0;i<n;i++)
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
}

bool	rkqs(double *y, double *dydx, int n, double &x, const double htry,
	      const double eps, double *yscal, double &hdid, double &hnext,
	      void derivs(const double, double *, double *))
{
	const double SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;
	int i;
	double errmax,h,htemp,xnew;
	
	h=htry;
	double yerr[n],ytemp[n];
	for (;;) {
		rkck(y,dydx,n,x,h,ytemp,yerr,derivs);
		errmax=0.0;
		for (i=0;i<n;i++) errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=(h >= 0.0 ? MAX(htemp,0.1*h) : MIN(htemp,0.1*h));
		xnew=x+h;
		if (xnew == x) return true;//nrerror("stepsize underflow in rkqs");
	}
	if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);
	else hnext=5.0*h;
	x += (hdid=h);
	for (i=0;i<n;i++) y[i]=ytemp[i];
	return false;
}

double dxsav;
int kmax,kount;
double *xp_p;
double **yp_p;

bool	odeint(double *ystart, int nvar, const double x1, const double x2, const double eps,
	       const double h1, const double hmin, int &nok, int &nbad,
	       void derivs(const double, double *, double *),
	       bool rkqs(double *, double *, int, double &, const double, const double,
			 double *, double &, double &, void (*)(const double, double *, double *)))
{
	const int MAXSTP=1000000;
	const double TINY=1.0e-30;
	int i,nstp;
	double xsav,x,hnext,hdid,h;
	
	double yscal[nvar],y[nvar],dydx[nvar];
	double *xp=xp_p;
	double **yp=yp_p;
	x=x1;
	h=SIGN(h1,x2-x1);
	nok = nbad = kount = 0;
	for (i=0;i<nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=0;nstp<MAXSTP;nstp++) {
		derivs(x,y,dydx);
		for (i=0;i<nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
			for (i=0;i<nvar;i++) yp[i][kount]=y[i];
			xp[kount++]=x;
			xsav=x;
		}
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		if(rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs))
			return true;
		if (hdid == h) ++nok; else ++nbad;
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=0;i<nvar;i++) ystart[i]=y[i];
			if (kmax != 0) {
				for (i=0;i<nvar;i++) yp[i][kount]=y[i];
				xp[kount++]=x;
			}
			return false;
		}
		if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
		h=hnext;
	}
	nrerror("Too many steps in routine odeint");
	return true;
}

bool	odeint(void f(double x, double *y, double *f),
		    double a, double b, double *y, int n, double err=1e-6) {
	int	dum;
	kmax = 0;
	return odeint(y,n,a,b, err, (b-a)*1.e-9, 0,dum,dum,f,rkqs);
}
