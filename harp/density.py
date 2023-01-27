import numpy as np

class gridclass(object):
	'''
	Convienience class to keep track of grided density data
	'''
	def __init__(self,origin,dxyz,nxyz,*args,**kwargs):
		self.origin = origin.astype('double')
		self.dxyz = dxyz.astype('double')
		self.nxyz = nxyz.astype('int64')

	def xyzmin(self):
		''' np.ndarray of min corner '''
		return self.origin
	def xyzmax(self):
		''' np.ndarray of max corner '''
		return self.origin + self.nxyz*self.dxyz

	def val(self):
		'''Returns all information (o,d,n)'''
		return self.origin,self.dxyz,self.nxyz

	def com(self):
		'''center of grid -- not weighted by density....'''
		return (self.xyzmin()+self.xyzmax())*.5

	def grid_ijk(self):
		'''return indices for grid points in ijk space'''
		gi,gj,gk = np.ogrid[:self.nxyz[0],:self.nxyz[1],:self.nxyz[2]]
		return gi,gj,gk
	def grid_xyz(self):
		'''return indices for grid points in xyz space'''
		gx,gy,gz = self.grid_ijk()
		gx = gx.astype('float64')*self.dxyz[0] + self.origin[0]
		gy = gy.astype('float64')*self.dxyz[1] + self.origin[1]
		gz = gz.astype('float64')*self.dxyz[2] + self.origin[2]
		return gx,gy,gz

def load(fname,options={}):
	'''
	Load densities (map, mrc, SF(cif), SF(mtz))
	input: string fname
	returns: grid(gridclass), density(np.ndarray)
	'''
	from . import io
	if fname.endswith('.map') or fname.endswith('.MAP') or fname.endswith('.mrc') or fname.endswith('.MRC') or fname.endswith('.map.gz') or fname.endswith('.MAP.gz'):
		grid,density = io.load_mrc(fname,options=options)
	elif fname.endswith('.cif') or fname.endswith('.CIF') or fname.endswith('.mtz') or fname.endswith('.MTZ') or fname.endswith('.ENT') or fname.endswith('.ent') or fname.endswith('.ENT.gz') or fname.endswith('.ent.gz'):
		grid,density = io.load_sf(fname,options=options)
	else:
		raise Exception('Do not know how to load %s'%(fname))
	density = np.ascontiguousarray(density).astype('double')
	return grid,density


def x2i(grid,point):
	''' transfer xyz to ijk on grid '''
	return ((point-grid.origin)//grid.dxyz).astype('int64')
def i2x(grid,point):
	''' transfer ijk on grid to xyz'''
	return (point*grid.dxyz + grid.origin).astype('float64')

def xyz_in_grid(grid,point):
	ijk = x2i(grid,point)
	if np.all(ijk < grid.nxyz) and np.all(ijk >= 0):
		return True
	return False

def subgrid_center(grid,centerxyz,halfwidthxyz):
	'''
	create a cubic grid centered around centerxyz padded by halfwidthxyz
	'''
	if not xyz_in_grid(grid,centerxyz):
		print('Centering: %s not in grid'%(str(centerxyz)))
	ijkmin = x2i(grid,centerxyz-halfwidthxyz)
	ijkmax = x2i(grid,centerxyz+halfwidthxyz)+1
	ijkmin[ijkmin < 0] = 0
	upper = ijkmax >= grid.nxyz
	ijkmax[upper] = grid.nxyz[upper]
	newnxyz = ijkmax - ijkmin
	neworigin = i2x(grid,ijkmin)
	return gridclass(neworigin,grid.dxyz,newnxyz)

def subgrid_extract(grid,data,subgrid):
	'''
	extract density for subgrid from (grid,data) pair
	'''
	ijkmin = x2i(grid,subgrid.origin)
	ijkmax = ijkmin+ subgrid.nxyz
	ijkmin[ijkmin<0] = 0
	ijkmax[ijkmax<0] = 0
	return data[ijkmin[0]:ijkmax[0],ijkmin[1]:ijkmax[1],ijkmin[2]:ijkmax[2]].copy()

def trim_density_to_mol(grid,data,mol,pad=8.):
	'''
	finds box that emcompasses mol's xyz positions, pads it isotropically with pad
	exports (grid,data) pair (using periodic boundary conditions) for new padded box from old data
	'''
	import numba as nb
	mmin = mol.xyz.min(0)-pad
	mmax = mol.xyz.max(0)+pad

	ijkmin = x2i(grid,mmin)
	ijkmax = x2i(grid,mmax)

	nxyz = ijkmax-ijkmin
	origin = i2x(grid,ijkmin)


	@nb.njit
	def fast_pbc_copy(data,ijkmin,nxyz,oldnxyz):
		d2 = np.zeros((int(nxyz[0]),int(nxyz[1]),int(nxyz[2])))

		for ii in range(nxyz[0]):
			i = ijkmin[0] + ii
			if i >= oldnxyz[0]:
				i -= oldnxyz[0]
			elif i < 0:
				i += oldnxyz[0]
			for jj in range(nxyz[1]):
				j = ijkmin[1] + jj
				if j >= oldnxyz[1]:
					j -= oldnxyz[1]
				elif j < 0:
					j += oldnxyz[1]
				for kk in range(nxyz[2]):
					k = ijkmin[2] + kk
					if k >= oldnxyz[2]:
						k -= oldnxyz[2]
					elif k < 0:
						k += oldnxyz[2]
					out = data[i,j,k]
					d2[ii,jj,kk] = data[i,j,k]
		return d2

	data = fast_pbc_copy(data,ijkmin,nxyz,grid.nxyz)
	grid.origin = origin
	grid.nxyz = nxyz
	return grid,data
