#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double a[5] = {.771058495001320e-04, -.133733772997339e-02, .323076579225834e-01, .479137145607681e-01, .128379167095513e+00};
double b[3] = {.301048631703895e-02, .538971687740286e-01, .375795757275549e+00};
double p[8] = {-1.36864857382717e-07, 5.64195517478974e-01, 7.21175825088309e+00, 4.31622272220567e+01, 1.52989285046940e+02, 3.39320816734344e+02, 4.51918953711873e+02, 3.00459261020162e+02};
double q[8] = {1.00000000000000e+00, 1.27827273196294e+01, 7.70001529352295e+01, 2.77585444743988e+02, 6.38980264465631e+02, 9.31354094850610e+02, 7.90950925327898e+02, 3.00459260956983e+02};
double r[5] = {2.10144126479064e+00, 2.62370141675169e+01, 2.13688200555087e+01, 4.65807828718470e+00, 2.82094791773523e-01};
double s[4] = {9.41537750555460e+01, 1.87114811799590e+02, 9.90191814623914e+01, 1.80124575948747e+01};

double cephes_erf(double x){
	/* This function is ripped from Cephes (https://www.netlib.org/cephes/) (MIT licensed)

	MIT License
	Permission is hereby granted, free of charge, to any person obtaining a copy of
	this software and associated documentation files (the "Software"), to deal in
	the Software without restriction, including without limitation the rights to
	use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
	of the Software, and to permit persons to whom the Software is furnished to
	do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
	PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
	OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
	SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
	*/

	double ax,bot,t,top,x2,out;
	double c = .564189583547756e0;

	ax = fabs(x);
	if (ax  > 0.5){ // line 45
		if (ax > 4.0){  // line 52
			if (ax >= 5.8) { // line 61
				// line 71-72
				out = copysign(1.0,x);
				return out;
			} else { // line 62
				// 62-69
				x2 = x*x;
				t = 1.0/x2;
				top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t + r[4];
				bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t + 1.;
				out = (c- top/(x2*bot))/ax;
				out = .5 + (.5-exp(-x2)*out);
				if (x < 0.){out = -out;}
				return out;
			}
		} else { //line 53
			// 53-59
			top = ((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax+p[6])*ax+p[7];
			bot = ((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax+q[6])*ax + q[7];
			out = .5 + (.5-exp(-x*x)*top/bot);
			if (x < 0.){out = -out;}
			return out;
		}
	} else{ // line 46-50
		t = x*x;
		top = ((((a[0]*t+a[1])*t + a[2])*t + a[3])*t + a[4]) + 1.0;
		bot = ((b[0]*t + b[1])*t + b[2])*t + 1.0;
		out = x*(top/bot);
		return out;
	}
}

void cephes_erf_nd(int N, double *x, double *y){
	for (int i=0; i<N; i++){
		y[i] = cephes_erf(x[i]);
	}
}


/*
 * Cephes Math Library, Release 2.3:  March, 1995
 * Copyright 1984, 1995 by Stephen L. Moshier
 */

/*
rewriting incbet from Cephes into numba jit-able function
CKT Feb 02, 2022
*/

// ### UNK(nown) mode from const.c
const double MAXGAM = 171.624376956302725;
const double MACHEP =  1.38777878078144567553E-17; //2**-56
const double MAXLOG =  8.8029691931113054295988E1;  //log(2**127)
const double MINLOG = -8.872283911167299960540E1; //log(2**-128)
const double MAXNUM =  1.701411834604692317316873e38; //2**127
const double big = 4.503599627370496e15;
const double biginv = 2.22044604925031308085e-16;


/* Power series for incomplete beta integral.
 * Use when b*x is small and x not too close to 1.  */
double pseries(double a, double b, double x){
	double s, t, u, v, n, t1, z, ai;

	ai = 1.0 / a;
	u = (1.0 - b) * x;
	v = u / (a + 1.0);
	t1 = v;
	t = u;
	n = 2.0;
	s = 0.0;
	z = MACHEP * ai;
	while (fabs(v) > z) {
		u = (n - b) * x / n;
		t *= u;
		v = t / (a + n);
		s += v;
		n += 1.0;
	}
	s += t1;
	s += ai;

	u = a * log(x);
	if ((a + b) < MAXGAM && fabs(u) < MAXLOG) {
		t = tgamma(a+b)/tgamma(a)/tgamma(b);
		s = s * t * pow(x, a);
	} else {
		t = lgamma(a+b)-lgamma(a)-lgamma(b) + u + log(s);
		if (t < MINLOG) {
			s = 0.0;
		} else {
			s = exp(t);
		}
	}
	return s;
}


/* Continued fraction expansion #1
 * for incomplete beta integral
 */

double incbcf(double a, double b, double x) {
	double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
	double k1, k2, k3, k4, k5, k6, k7, k8;
	double r, t, ans, thresh;
	int n;

	k1 = a;
	k2 = a + b;
	k3 = a;
	k4 = a + 1.0;
	k5 = 1.0;
	k6 = b - 1.0;
	k7 = k4;
	k8 = a + 2.0;

	pkm2 = 0.0;
	qkm2 = 1.0;
	pkm1 = 1.0;
	qkm1 = 1.0;
	ans = 1.0;
	r = 1.0;
	n = 0;
	thresh = 3.0 * MACHEP;
	while (n < 300){
		n += 1;

		xk = -(x * k1 * k2) / (k3 * k4);
		pk = pkm1 + pkm2 * xk;
		qk = qkm1 + qkm2 * xk;
		pkm2 = pkm1;
		pkm1 = pk;
		qkm2 = qkm1;
		qkm1 = qk;

		xk = (x * k5 * k6) / (k7 * k8);
		pk = pkm1 + pkm2 * xk;
		qk = qkm1 + qkm2 * xk;
		pkm2 = pkm1;
		pkm1 = pk;
		qkm2 = qkm1;
		qkm1 = qk;

		if (qk != 0){
			r = pk / qk;
		}

		if (r != 0){
			t = fabs((ans - r) / r);
			ans = r;
		} else {
			t = 1.0;
		}

		if (t < thresh){
			n = 301;
		}

		k1 += 1.0;
		k2 += 1.0;
		k3 += 2.0;
		k4 += 2.0;
		k5 += 1.0;
		k6 -= 1.0;
		k7 += 2.0;
		k8 += 2.0;

		if ((fabs(qk) + fabs(pk)) > big) {
			pkm2 *= biginv;
			pkm1 *= biginv;
			qkm2 *= biginv;
			qkm1 *= biginv;
		}

		if ((fabs(qk) < biginv) || (fabs(pk) < biginv)) {
			pkm2 *= big;
			pkm1 *= big;
			qkm2 *= big;
			qkm1 *= big;
		}
	}
	return ans;
}


/* Continued fraction expansion #2
 * for incomplete beta integral
 */

double incbd(double a, double b,double  x){
	double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
	double k1, k2, k3, k4, k5, k6, k7, k8;
	double r, t, ans, z, thresh;
	int n;

	k1 = a;
	k2 = b - 1.0;
	k3 = a;
	k4 = a + 1.0;
	k5 = 1.0;
	k6 = a + b;
	k7 = a + 1.0;
	k8 = a + 2.0;

	pkm2 = 0.0;
	qkm2 = 1.0;
	pkm1 = 1.0;
	qkm1 = 1.0;
	z = x / (1.0 - x);
	ans = 1.0;
	r = 1.0;
	n = 0;
	thresh = 3.0 * MACHEP;
	while (n < 300){
		n += 1;

		xk = -(z * k1 * k2) / (k3 * k4);
		pk = pkm1 + pkm2 * xk;
		qk = qkm1 + qkm2 * xk;
		pkm2 = pkm1;
		pkm1 = pk;
		qkm2 = qkm1;
		qkm1 = qk;

		xk = (z * k5 * k6) / (k7 * k8);
		pk = pkm1 + pkm2 * xk;
		qk = qkm1 + qkm2 * xk;
		pkm2 = pkm1;
		pkm1 = pk;
		qkm2 = qkm1;
		qkm1 = qk;

		if (qk != 0){
			r = pk / qk;
		}

		if (r != 0){
			t = fabs((ans - r) / r);
			ans = r;
		} else {
			t = 1.0;
		}

		if (t < thresh){
			n = 301;
		}

		k1 += 1.0;
		k2 -= 1.0;
		k3 += 2.0;
		k4 += 2.0;
		k5 += 1.0;
		k6 += 1.0;
		k7 += 2.0;
		k8 += 2.0;

		if ((fabs(qk) + fabs(pk)) > big){
			pkm2 *= biginv;
			pkm1 *= biginv;
			qkm2 *= biginv;
			qkm1 *= biginv;
		}
		if ((fabs(qk) < biginv) || (fabs(pk) < biginv)){
			pkm2 *= big;
			pkm1 *= big;
			qkm2 *= big;
			qkm1 *= big;
		}
	}
	return ans;
}

double betainc(double aa, double bb, double xx){
	double a, b, t, x, xc, w, y;
	int flag;

	if (aa <= 0.0 || bb <= 0.0){
		return NAN;
	}

	if ((xx <= 0.0) || (xx >= 1.0)){
		if (xx == 0.0){
			return 0.0;
		}
		if (xx == 1.0){
			return 1.0;
		}
		return NAN;
	}

	flag = 0;
	if ((bb * xx) <= 1.0 && xx <= 0.95){
		t = pseries(aa, bb, xx);
		return t;
	}

	w = 1.0 - xx;

	/* Reverse a and b if x is greater than the mean. */
	if (xx > (aa / (aa + bb))){
		flag = 1;
		a = bb;
		b = aa;
		xc = xx;
		x = w;
	} else {
		a = aa;
		b = bb;
		xc = w;
		x = xx;
	}

	if (flag == 1 && (b * x) <= 1.0 && x <= 0.95){
		t = pseries(a, b, x);
		if (t <= MACHEP){
			t = 1.0 - MACHEP;
		} else {
			t = 1.0 - t;
		}
		return t;
	}

	/* Choose expansion for better convergence. */
	y = x * (a + b - 2.0) - (a - 1.0);
	if (y < 0.0){
		w = incbcf(a, b, x);
	} else {
		w = incbd(a, b, x) / xc;
	}

	/* Multiply w by the factor
	 * a	  b   _			 _	 _
	 * x  (1-x)   | (a+b) / ( a | (a) | (b) ) .   */

	y = a * log(x);
	t = b * log(xc);
	if ((a + b) < MAXGAM && fabs(y) < MAXLOG && fabs(t) < MAXLOG){
		t = pow(xc, b);
		t *= pow(x, a);
		t /= a;
		t *= w;
		t *= tgamma(a+b)/tgamma(a)/tgamma(b);
		if (flag == 1){
			if (t <= MACHEP){
				t = 1.0 - MACHEP;
			} else {
				t = 1.0 - t;
			}
		}
		return t;
	}

	/* Resort to logarithms.  */
	y += t - lgamma(a) - lgamma(b) + lgamma(a+b);
	y += log(w / a);
	if (y < MINLOG){
		t = 0.0;
	} else {
		t = exp(y);
	}

	// done
	if (flag == 1){
		if (t <= MACHEP){
			t = 1.0 - MACHEP;
		} else {
			t = 1.0 - t;
		}
	}
	return t;
}

void zero_out(double *density, long *nxyz, int *dboundsijk) {
	for (int ii = dboundsijk[0]; ii < dboundsijk[3]+1; ii++){
		for (int ji = dboundsijk[1]; ji < dboundsijk[4]+1; ji++){
			for (int ki = dboundsijk[2]; ki < dboundsijk[5]+1; ki++){
				density[ii*(int)nxyz[1]*(int)nxyz[2]+ji*(int)nxyz[2]+ki] = 0.;
			}
		}
	}
}

void render_atoms_gaussian(int natoms, double *dmodel, int *dboundsijk, double *weights, double *origin, double *dxyz, long *nxyz, double *xyz, double *sigmas, int nsigma_cutoff, double offset){
	// dmodel (*nxyz) - should be made in wrapper
	// dboundsijk (6) where :,0:3 is min and :,3: is max of rendered box - should be made in wrapper
	// weights (natoms)
	// origin (3)
	// dxyz (3)
	// nxyz (3)
	// xyz (natoms,3)
	// sigmas (natoms)
	// nsigma_cutoff (1))
	// offset (1))

	double r2ss = 0.;

	// face or edge centered...
	double offsethigh = 1. - offset;
	double offsetlow  = 0. - offset;
	double oxh = offsethigh*dxyz[0];
	double oxl = offsetlow*dxyz[0];
	double oyh = offsethigh*dxyz[1];
	double oyl = offsetlow*dxyz[1];
	double ozh = offsethigh*dxyz[2];
	double ozl = offsetlow*dxyz[2];

	double xyzi[3] = {0.,0.,0.};
	double xyzimin[3] = {0.,0.,0.};
	double xyzimax[3] = {0.,0.,0.};
	int ijkimin[3] = {0,0,0};
	int ijkimax[3] = {0,0,0};
	double di,dj,dk,ci,cj,ck,c;

	// reset residue boundaries
	for (int i=0; i<3; i++){
		dboundsijk[i] = nxyz[i]-1;
		dboundsijk[3+i] = 0; // yes it should start out backwards...
	}

	for (int atomi = 0; atomi < natoms; atomi++){

		// figure out atom boundaries
		for (int i = 0; i<3; i++){
			xyzi[i] = xyz[atomi*3 + i];
			xyzimin[i] = xyzi[i] - sigmas[atomi] * (double)nsigma_cutoff;
			xyzimax[i] = xyzi[i] + sigmas[atomi] * (double)nsigma_cutoff;

			ijkimin[i] = (int)((xyzimin[i]-origin[i])/dxyz[i]);
			ijkimax[i] = (int)((xyzimax[i]-origin[i])/dxyz[i]);

			if (ijkimin[i] < 0){
				ijkimin[i] = 0;
			}
			if (ijkimax[i] > nxyz[i]-1){
				ijkimax[i] = nxyz[i]-1;
			}

			if (ijkimin[i] < dboundsijk[i]){
				dboundsijk[i] = ijkimin[i];
			}
			if (ijkimax[i] > dboundsijk[3+i]){
				dboundsijk[3+i] = ijkimax[i];
			}
		}

		r2ss = sqrt(2.*sigmas[atomi]*sigmas[atomi]);
		for (int ii = ijkimin[0]; ii < ijkimax[0]+1; ii++){
			di = ((double)ii*dxyz[0] + origin[0]) - xyzi[0]; // distances from mu
			ci = (cephes_erf((di+oxh)/r2ss) - cephes_erf((di+oxl)/r2ss));
			for (int ji = ijkimin[1]; ji < ijkimax[1]+1; ji++){
				dj = ((double)ji*dxyz[1] + origin[1]) - xyzi[1];
				cj = (cephes_erf((dj+oyh)/r2ss) - cephes_erf((dj+oyl)/r2ss));
				for (int ki = ijkimin[2]; ki < ijkimax[2]+1; ki++){
					dk = ((double)ki*dxyz[2] + origin[2]) - xyzi[2];
					ck = (cephes_erf((dk+ozh)/r2ss) - cephes_erf((dk+ozl)/r2ss));
					c = ci*cj*ck;
					if (isfinite(c)){
						dmodel[ii*(int)nxyz[1]*(int)nxyz[2]+ji*(int)nxyz[2]+ki] +=  .125*c*weights[atomi];
					}
				}
			}
		}
	}
	//## Returns:
	//## dmodel is nxyz sized
	//## dboundsijk is (6) where values are min/mx bins containing density. ranges need max+1
}



double ln_evidence(double *x, double *y, long *nxyz, double lnprior_loc, double lnprior_scale) {
	// nxyz is for the full array
	// x and y are full data
	// calculation is over local ROI
	// evidence calculation 2.2.2

	// Initialize
	double N = (double)(nxyz[0]*nxyz[1]*nxyz[2]);
	double M = N/2.-1.;
	double Ex = 0.;
	double Ey = 0.;
	double Exx = 0.;
	double Eyy = 0.;
	double Exy = 0.;
	double xi,yi,normalization;

	for (int ii = 0; ii < (int)nxyz[0]; ii++){
		for (int jj = 0; jj < (int)nxyz[1]; jj++){
			for (int kk = 0; kk < (int)nxyz[2]; kk++){
				xi = x[ii*(int)nxyz[1]*(int)nxyz[2]+jj*(int)nxyz[2]+kk];
				yi = y[ii*(int)nxyz[1]*(int)nxyz[2]+jj*(int)nxyz[2]+kk];
				Ex += xi;
				Ey += yi;
				Exx += xi*xi;
				Exy += xi*yi;
				Eyy += yi*yi;
			}
		}
	}

	// the above are just sums right now
	normalization = 1.;//sqrt(Exx/N - Ex*Ex/N/N);
	Ex /= N*normalization;
	Ey /= N;
	Exx /= N*normalization*normalization;
	Exy /= N*normalization;
	Eyy /= N;

	double vx = Exx-Ex*Ex;
	double vy = Eyy-Ey*Ey;
	double vxy = Exy-Ex*Ey;

	if (vx <=0. || vy <= 0. || vy == 0. || isnan(vx)){
		return -INFINITY;
	}

	double r = vxy/sqrt(vx*vy);
	double r2 = r*r;

	double out;
	out = lgamma(M) -.5*N*log(N) -.5*log(vx) -log(2.) -2.*lnprior_loc -lnprior_scale;
	out += -M*log(M_PI) - M*log(vy) - M*log(1.-r2) + log(1+r/fabs(r)*betainc(.5,M,r2));

	return out;
}

void harp(int natoms, int nadfs, int nblobs, double *ln_out, double *density, double *com, double *origin, double *dxyz, long *nxyz, double *xyz, double *weights, double *blobs, double *adfs, int nsigma_cutoff, double offset) {
	// ln_out (nblobs+1)
	// density (*nxyz)
	// origin (3)
	// dxyz (3)
	// nxyz (3)
	// blob_sigmas (nblobs)
	// atom_sigma (1)
	// nsigma_cutoff (1)

	double lnprior_loc = 5.*log(10.)+log(2.);
	double lnprior_scale = log(2*3.*log(10));

	// initialize
	double *dmodel, *sigmas_atomic;
	dmodel = (double *)malloc(sizeof(double)*((int)nxyz[0]+1)*((int)nxyz[1]+1)*((int)nxyz[2]+1));
	sigmas_atomic = (double *)malloc(sizeof(double)*natoms);
	double sigmas_spheric[1] = {0.};
	int dboundsijk[6] = {0,0,0,(int)nxyz[0],(int)nxyz[1],(int)nxyz[2]};

	// set some values
	zero_out(dmodel,nxyz,dboundsijk); // maybe this isn't necessary... probably is.
	double SA_weight[1] = {0.}; // the weight for a superatom is the sum of all the weights for the individual atoms
	for (int atomi=0; atomi<natoms; atomi++){
		SA_weight[0] += weights[atomi];
	}


	// do the atomic models
	for (int adfi = 0; adfi < nadfs; adfi++){
		// prep
		for (int atomi = 0; atomi < natoms; atomi++){
			sigmas_atomic[atomi] = adfs[adfi];
		}
		for (int i=0; i<3; i++){
			dboundsijk[i] = 0;
			dboundsijk[i+3] = (int)nxyz[i];
		}
		zero_out(dmodel,nxyz,dboundsijk);
		// calc
		render_atoms_gaussian(natoms,dmodel,dboundsijk,weights,origin,dxyz,nxyz,xyz,sigmas_atomic,nsigma_cutoff,offset);
		ln_out[adfi] = ln_evidence(dmodel,density,nxyz,lnprior_loc,lnprior_scale);
	}

	// do the spheric models
	for (int blobi = 0; blobi < nblobs; blobi++){
		// prep
		sigmas_spheric[0] = blobs[blobi];
		for (int i=0; i<3; i++){
			dboundsijk[i] = 0;
			dboundsijk[i+3] = (int)nxyz[i];
		}
		zero_out(dmodel,nxyz,dboundsijk);
		// calc
		render_atoms_gaussian(1,dmodel,dboundsijk,SA_weight,origin,dxyz,nxyz,com,sigmas_spheric,nsigma_cutoff,offset); // This would grow, but the density is clipped in the wrapper, so it should be constant.
		ln_out[nadfs+blobi] = ln_evidence(dmodel,density,nxyz,lnprior_loc,lnprior_scale);
	}

	 // -------  DEBUGGING -------- //
	 //for (int i=0;i<nblobs+1;i++){
		 //printf("%.8f,",ln_out[i]);
	 //} printf("\n");
	 //exit(1);

	// free up some stuff
	free(dmodel);
	free(sigmas_atomic);
}
