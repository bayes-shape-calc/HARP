import numpy as np
np.seterr(all='ignore')
from .jit_math import betainc, gammaln
import numba as nb

lnprior_factor_location = 5.*np.log(10) +np.log(2.)## +10^5 to -10^5: 1/(\delta \mu)
lnprior_factor_scale = np.log(2 * 3*np.log(10.)) ## +/-3 decades: 1/(\delta ln(\tau))

@nb.njit(cache=True)
def ln_evidence(x,y):
	### 2.2.2

	N = float(x.size)
	M = N/2.-1.
	Ex = np.mean(x)
	Ey = np.mean(y)
	Exx = np.mean(x*x)
	Eyy = np.mean(y*y)
	Exy = np.mean(x*y)
	vx = Exx-Ex*Ex
	vy = Eyy-Ey*Ey
	vxy = Exy-Ex*Ey

	## underflow protection
	if vx <= 0 or vy <= 0 or vxy == 0: ## vxy can be negative
		return -np.inf

	r = vxy/np.sqrt(vx*vy)
	r2 = r*r
	out = gammaln(M) -N/2.*np.log(N) -.5*np.log(vx) -np.log(2) -2.*lnprior_factor_location -lnprior_factor_scale -M*np.log(np.pi) -M*np.log(vy) -M*np.log(1.-r2) + np.log(1.+r/np.abs(r)*betainc(.5,M,r2))
	return out


def ln_evidence_scipy(x,y):
	from scipy.special import betainc as ssbetainc
	### 2.2.2

	N = float(x.size)
	M = N/2.-1.
	Ex = np.mean(x)
	Ey = np.mean(y)
	Exx = np.mean(x*x)
	Eyy = np.mean(y*y)
	Exy = np.mean(x*y)
	vx = Exx-Ex*Ex
	vy = Eyy-Ey*Ey
	vxy = Exy-Ex*Ey

	## underflow protection
	if vx <= 0 or vy <= 0 or vxy == 0: ## vxy can be negative
		return -np.inf

	r = vxy/np.sqrt(vx*vy)
	r2 = r*r
	out = gammaln(M) -N/2.*np.log(N) -.5*np.log(vx) -np.log(2) -2.*lnprior_factor_location -lnprior_factor_scale -M*np.log(np.pi) -M*np.log(vy) -M*np.log(1.-r2) + np.log(1.+r/np.abs(r)*ssbetainc(.5,M,r2))
	return out


@nb.njit(cache=True)
def null_ln_evidence(y):
	## 2.2.4

	N = float(y.size)
	M = N/2.-1.
	Ey = np.mean(y)
	Eyy = np.mean(y*y)
	vy = Eyy-Ey*Ey

	## underflow protection
	if vy <= 0:
		return -np.inf

	out = gammaln(M+.5) -N/2.*np.log(N) -lnprior_factor_location -lnprior_factor_scale -(M+.5)*np.log(np.pi) -(M+.5)*np.log(vy)
	return out
