import numpy as np
np.seterr(all='ignore')
import numba as nb
from math import erf
import ctypes
import os
lib_path = os.path.dirname(os.path.abspath(__file__))

_render_model = None
version = ''

## Load the Python version

@nb.njit(cache=True) ## don't parallelize this b/c otherwise race conditions
def _render_model_python(origin,dxyz,nxyz,xyz,weights,sigma,nsigma_cutoff,offset):
	'''
	sum up for several spots: integrate gaussian over voxels close to spot
	gives a decent speed up
	'''

	dmodel = np.zeros((nb.int64(nxyz[0]),nb.int64(nxyz[1]),nb.int64(nxyz[2])))


	## face or edge centered...
	offsethigh = 1. - offset
	offsetlow  = 0. - offset
	oxh = offsethigh*dxyz[0]
	oxl = offsetlow*dxyz[0]
	oyh = offsethigh*dxyz[1]
	oyl = offsetlow*dxyz[1]
	ozh = offsethigh*dxyz[2]
	ozl = offsetlow*dxyz[2]

	for atomi in range(xyz.shape[0]):
		r2ss = np.sqrt(2.*sigma[atomi]*sigma[atomi])
		xyzi = xyz[atomi]
		xyzimin = xyzi - nsigma_cutoff*sigma[atomi]
		xyzimax = xyzi + nsigma_cutoff*sigma[atomi]

		ijkimin = ((xyzimin - origin)//dxyz)//1
		ijkimax = ((xyzimax - origin)//dxyz)//1

		for ii in range(max(0,ijkimin[0]),min(ijkimax[0]+1,dmodel.shape[0])):
				di = (ii*dxyz[0] + origin[0]) - xyzi[0] ## distances from mu
				for ji in range(max(0,ijkimin[1]),min(ijkimax[1]+1,dmodel.shape[1])):
						dj = (ji*dxyz[1] + origin[1]) - xyzi[1]
						for ki in range(max(0,ijkimin[2]),min(ijkimax[2]+1,dmodel.shape[2])):
								dk = (ki*dxyz[2] + origin[2]) - xyzi[2]
								c = (erf((di+oxh)/r2ss) - erf((di+oxl)/r2ss))
								c *= (erf((dj+oyh)/r2ss) - erf((dj+oyl)/r2ss))
								c *= (erf((dk+ozh)/r2ss) - erf((dk+ozl)/r2ss))
								dmodel[ii,ji,ki] += .125*c*weights[atomi]
	return dmodel

def use_python():
	global _render_model,version
	version = 'python'
	_render_model = _render_model_python
use_python()



### Try to load the C version, which might not be compiled
try:
	sopath = os.path.join(lib_path,'render_model.so')
	lib_c = np.ctypeslib.load_library(sopath, '.') ## future-self: the library has to end in .so ....

	lib_c.render_atoms_gaussian.argtypes = [
		ctypes.c_int,
		np.ctypeslib.ndpointer(dtype = np.float64),
		np.ctypeslib.ndpointer(dtype = np.int32),
		np.ctypeslib.ndpointer(dtype = np.float64),
		np.ctypeslib.ndpointer(dtype = np.float64),
		np.ctypeslib.ndpointer(dtype = np.float64),
		np.ctypeslib.ndpointer(dtype = np.int64),
		np.ctypeslib.ndpointer(dtype = np.float64),
		np.ctypeslib.ndpointer(dtype = np.float64),
		ctypes.c_int,
		ctypes.c_double
		]
	lib_c.render_atoms_gaussian.restype  = ctypes.c_void_p

	lib_c.harp.argtypes = [
		ctypes.c_int,
		ctypes.c_int,
		ctypes.c_int,
		np.ctypeslib.ndpointer(dtype = np.float64),
		np.ctypeslib.ndpointer(dtype = np.float64),
		np.ctypeslib.ndpointer(dtype = np.float64),
		np.ctypeslib.ndpointer(dtype = np.float64),
		np.ctypeslib.ndpointer(dtype = np.float64),
		np.ctypeslib.ndpointer(dtype = np.int64),
		np.ctypeslib.ndpointer(dtype = np.float64),
		np.ctypeslib.ndpointer(dtype = np.float64),
		np.ctypeslib.ndpointer(dtype = np.float64),
		np.ctypeslib.ndpointer(dtype = np.float64),
		ctypes.c_int,
		ctypes.c_double
		]
	lib_c.harp.restype  = ctypes.c_void_p

	lib_c.ln_evidence.argtypes = [
		np.ctypeslib.ndpointer(dtype = np.float64),
		np.ctypeslib.ndpointer(dtype = np.float64),
		np.ctypeslib.ndpointer(dtype = np.int64),
		ctypes.c_double,
		ctypes.c_double
		]
	lib_c.ln_evidence.restype  = ctypes.c_double

	def _ln_evidence_c(x,y):
		from .evidence import lnprior_factor_location, lnprior_factor_scale
		nxyz = np.array(x.shape)
		return lib_c.ln_evidence(x,y,nxyz,lnprior_factor_location,lnprior_factor_scale)

	def _render_model_c(origin,dxyz,nxyz,xyz,weights,sigmas,nsigma_cutoff,offset):
		natoms,_ = xyz.shape
		out = np.zeros(nxyz,dtype='double',order='C')
		if not xyz.flags.c_contiguous: ## FML.... this was a terrible error to track down
			xyz = np.ascontiguousarray(xyz)
		dbounds = np.zeros(6,dtype='int32')
		lib_c.render_atoms_gaussian(natoms,out,dbounds,weights,origin,dxyz,nxyz,xyz,sigmas,int(nsigma_cutoff),float(offset))
		return out

	def _harp(density,com,origin,dxyz,nxyz,xyz,weights,adfs,blobs,nsigma_cutoff,offset):
		nblobs = blobs.size
		nadfs = adfs.size
		natoms,_ = xyz.shape
		ln_ev = np.zeros(nblobs+nadfs)

		if not xyz.flags.c_contiguous: ## FML.... this was a terrible error to track down
			xyz = np.ascontiguousarray(xyz)

		lib_c.harp(natoms,nadfs,nblobs,ln_ev,density,com,origin,dxyz,nxyz,xyz,weights,blobs,adfs,int(nsigma_cutoff),float(offset))
		return ln_ev

	def use_c():
		global _render_model,version
		version = 'C'
		_render_model = _render_model_c
	use_c()

except Exception as e:
	print(e)
	print('failed to load c')


def density_point(grid,point,weight=1.,sigma=.7,nsigma=5,offset=0.5):
	'''
	single point calculation
	note: It seems like face centered is used -- 4f35 gives higher <b> of .47 to .29 for face vs edge.
	'''
	weights = np.array((weight,))
	sigmas = np.array((sigma,))
	return _render_model(grid.origin,grid.dxyz,grid.nxyz,point[None,:],weights,sigmas,nsigma,offset)

def density_atoms(grid,xyzs,weights=None,sigma=.7,nsigma=5,offset=0.5):
	'''
	multiple point calculation
	'''
	natoms,_ = xyzs.shape
	if weights is None:
		weights = np.ones(natoms)
	if np.size(sigma)==1:
		sigmas = np.zeros((natoms))+sigma
	else:
		sigmas = sigma
	return _render_model(grid.origin,grid.dxyz,grid.nxyz,xyzs,weights,sigmas,nsigma,offset)
