#include "nrecipe.cpp"
#include "odeint.cpp"

struct BAO {
	int Th_Lmax, Nu_Lmax, Th_idx, Nu_idx, SIZE;
	double omega_ph, omega_nu, omega_dm, omega_br, omega_r, omega_m, lambda, omega_k;
	double k, a0, x1, c1, c2, Psi0, x, *y;
	enum { Phi_idx, delta_dm_idx, delta_br_idx, v_dm_idx, v_br_idx };
	BAO(double om_ph, double om_nu, double om_dm, double om_br, double lm);
	~BAO();
	double differential_optical_depth(double a);
	double horizon_radius(double a);
	void boltzmann_eq(double x, double *y, double *dydx, bool tight);
	void init(double k);
	void evolv(double a);
	double& delta_br() { return y[delta_br_idx]; }
	double& delta_dm() { return y[delta_dm_idx]; }
	double& v_br() { return y[v_br_idx]; }
	double& v_dm() { return y[v_dm_idx]; }
	double& Phi()  { return y[Phi_idx]; }
	double& Nu(int i) { return y[Nu_idx + i]; }
	double& Th(int i) { return y[Th_idx + i]; }
} *bao_p;

BAO::BAO(double om_ph, double om_nu, double om_dm, double om_br, double lm=1e30) {
	bao_p = this;
	double d_tau_thresh(0.75e2*om_br/om_ph);
	Th_Lmax = 20;
	Nu_Lmax = 20;
	a0 = 1.e-9;
	Psi0 = 1;
	omega_ph = om_ph;
	omega_nu = om_nu;
	omega_dm = om_dm;
	omega_br = om_br;
	omega_r = omega_ph + omega_nu;
	omega_m = omega_dm + omega_br;
	lambda  = (lm > 0.99e30 ? 1 - omega_r - omega_m : lm);
	omega_k = (lm > 0.99e30 ? 0 : 1 - omega_r - omega_m - lm);
	c1 = 15 + 1.8*log(omega_br);
	c2 = pow(omega_br, 0.43)*(c1+1)*1.e-3;
	x1 = -log(d_tau_thresh/c2*pow(1e3,c1))/(c1+1);
	Nu_idx = v_br_idx + 1;
	Th_idx = Nu_idx + Nu_Lmax;
	SIZE = Th_idx + Th_Lmax;
	y = new double[SIZE];
}

BAO::~BAO() { delete[] y; }

double BAO::differential_optical_depth(double a) {
	double z((1/a-1)*1.e-3);
	return (z<0 ? 0 : c2*pow(z,c1)/a);
}

double BAO::horizon_radius(double a) {
	return a/sqrt(omega_r + a*(omega_m + a*(omega_k + lambda*a*a)));
}

void blzm_p(double x, double *y, double *f) { bao_p->boltzmann_eq(x,y,f,false); }
void tcpl_p(double x, double *y, double *f) { bao_p->boltzmann_eq(x,y,f,true); }

void BAO::boltzmann_eq(double x, double *y, double *f, bool tight) {
	double a = exp(x);
	double& Phi = y[Phi_idx];
	double& delta_dm = y[delta_dm_idx];
	double& delta_br = y[delta_br_idx];
	double& v_dm = y[v_dm_idx];
	double& v_br = y[v_br_idx];
	double *Th = &y[Th_idx];
	double *Nu = &y[Nu_idx];
	double& d_Phi = f[Phi_idx];
	double& d_delta_dm = f[delta_dm_idx];
	double& d_delta_br = f[delta_br_idx];
	double& d_v_dm = f[v_dm_idx];
	double& d_v_br = f[v_br_idx];
	double *d_Nu = &f[Nu_idx];
	double *d_Th = &f[Th_idx];
	double R = 0.75*omega_br*a/omega_ph;
	double T0i = 1.5*(omega_dm*v_dm + omega_br*v_br
			+ 4*(omega_ph*Th[1] + omega_nu*Nu[1])/a)/a/k/k;
	double Tij = 12*(omega_ph*Th[2] + omega_nu*Nu[2])/a/a/k/k;
	double Psi = -Tij - Phi;
	double kR = k*horizon_radius(a);
	d_Phi = -T0i*kR + Psi;
	d_delta_dm = -kR*v_dm - 3*d_Phi;
	d_delta_br = -kR*v_br - 3*d_Phi;
	d_v_dm = kR*Psi - v_dm;
	d_Nu[0] = -kR*Nu[1] - d_Phi;
	d_Nu[1] = kR*(Nu[0] - 2*Nu[2] + Psi)/3;
	d_Th[0] = -kR*Th[1] - d_Phi;
	for(int l=2; l<Nu_Lmax; l++) {
		if(l+1 < Nu_Lmax) d_Nu[l] = Nu[l+1];
		else d_Nu[l] = (2*l+1)*Nu[l]/kR - Nu[l-1];
		d_Nu[l] = kR*(l*Nu[l-1] - (l+1)*d_Nu[l])/(2*l+1);
	}
	if(tight) {
		d_v_br = kR*Psi + (kR*Th[0] - R*v_br)/(1+R);
		d_Th[1] = d_v_br/3;
		for(int l=2; l<Th_Lmax; l++) d_Th[l] = 0;
	}
	else {	double d_tau = differential_optical_depth(a);
		d_v_br = kR*Psi - v_br + d_tau*(3*Th[1] - v_br)/R;
		d_Th[1] = kR*(Th[0] - 2*Th[2] + Psi)/3 + d_tau*(v_br/3 - Th[1]);
		d_Th[2] = kR*(2*Th[1] - 3*Th[3])/5 - 0.9*d_tau*Th[2];
		for(int l=3; l<Th_Lmax; l++) {
			if(l+1 < Th_Lmax) d_Th[l] = Th[l+1];
			else d_Th[l] = (2*l+1)*Th[l]/kR - Th[l-1];
			d_Th[l] = kR*(l*Th[l-1] - (l+1)*d_Th[l])/(2*l+1) - d_tau*Th[l];
		}
	}
}

void BAO::init(double k1) {
	k = k1;
	double kR = k*horizon_radius(a0);
	x = log(a0);
	Phi() = -(1 + 0.4*omega_nu/omega_r)*Psi0;
	delta_dm() = delta_br() = -1.5*Psi0;
	Th(0) = Nu(0) = -0.5*Psi0;
	v_dm() = v_br() = 0.5*kR*Psi0;
	Th(1) = Nu(1) = v_br()/3;
	Nu(2) = kR*kR*Psi0/30.;
	for(int l=2; l<Th_Lmax; l++) Th(l) = 0;
	for(int l=3; l<Nu_Lmax; l++) Nu(l) = 0;
}

void BAO::evolv(double a) {
	double x0(x);
	x = log(a);
	if(x0 <= x1 && x <= x1)
		odeint(tcpl_p, x0, x, y, SIZE);
	else if(x0 >= x1 && x >= x1)
		odeint(blzm_p, x0, x, y, SIZE);
	else if(x0 <= x1 && x >= x1) {
		odeint(tcpl_p, x0, x1, y, SIZE);
		odeint(blzm_p, x1, x, y, SIZE);
	}
	else {	odeint(blzm_p, x0, x1, y, SIZE);
		odeint(tcpl_p, x1, x, y, SIZE);
	}	
}

#include<cstdio>

main() {
	FILE *fp = fopen("bao.txt", "w");
	double h(0.7), om_dm(0.233), om_br(0.0463), ns(0.972);
	double om_ph = 2.471e-5/h/h;
	double om_nu = 3*7/8.*pow(4/11., 4/3.)*om_ph;
	BAO bao(om_ph, om_nu, om_dm, om_br);
	int n(100);
	double x1(1.e-3), x2(1e0), dx(pow(x2/x1, 1./n));
	for(int i=0; i<=n; i++) {
		double x = x1*pow(dx,i); // wave number / h/Mpc
		double k = x*2997.9;	  // wave number / H0/c
		bao.init(k);
		bao.evolv(1);
		double	y = pow(bao.delta_dm(), 2);//*pow(k,ns-1);
		fprintf(fp, "%e\t%e\n", x, y);
		printf("%e\t%e\n", x, y);
	}
}
