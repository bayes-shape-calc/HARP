import numpy as np
from ..density import gridclass

def load_mrc(fname_mrc,options={}):
	import mrcfile
	## Load file
	mrcf = mrcfile.open(fname_mrc)
	## Get header info
	h = mrcf.header
	## I'm not handling non-squre grids yet
	if not np.all([h.cellb[angle] == 90 for angle in ['alpha','beta','gamma']]):
		raise Exception('MRC: non-square grid')

	## get information from header
	ncrs = np.array((h['nx'],h['ny'],h['nz'])).astype('int64') ## whole volume
	ncrs_start = np.array((h['nxstart'],h['nystart'],h['nzstart'])).astype('int64') ## whole volume
	unitcell_lengths = np.array([h['cella'][i] for i in ['x','y','z']]) ## unitcell
	unitcell_nxyz = np.array((h['mx'],h['my'],h['mz'])) ## unitcell
	## Origin is often messed up according to MRC paper Cheng et al JStructBiol 2015 - so ignore it
	# origin = np.array([h['origin'][i] for i in ['x','y','z']])

	## reordering from CRS to XYZ
	crs = np.array((h['mapc'],h['mapr'],h['maps'])) ## order
	ijk = crs - 1 ## use zero index
	asort = ijk.argsort() ## get reordering vector

	## calculate xyz grid info
	dxyz = unitcell_lengths/unitcell_nxyz.astype('double')
	nxyz_start = ncrs_start[asort]
	origin = (nxyz_start*dxyz).astype('double')
	nxyz = ncrs[asort]

	## reorder data from CRS to XYZ! Note: mrcfile library provides .data as [z,y,x]..... so rev to get x,y,z
	# data = mrcf.data.copy()
	# del mrcf
	swap = np.array((0,1,2))[::-1] ## reverse because mrcfile gives zyx and we want xyz
	data = np.moveaxis(mrcf.data,swap,ijk) ##  swap to account for the map order in mrc standard
	del mrcf

	if not np.all(nxyz==np.array(data.shape)):
		raise Exception("Density loading failed check")

	return gridclass(origin,dxyz,nxyz),data

def save_mrc(fname_mrc,grid,data,origin_at_zero=False,reverse=True):
	import mrcfile
	with mrcfile.new(fname_mrc,overwrite='True') as mrcf:
		nstart = np.array((0,0,0)).astype('int64') if origin_at_zero else (grid.origin//grid.dxyz).astype('int64')
		origin = np.array((0.,0.,0.)) if origin_at_zero else grid.origin
		dxyz = grid.dxyz
		nxyz = grid.nxyz
		# if extend_to_zero: ##extend_to_zero
		# 	if np.any(nstart < 0):
		# 		raise Exception('Save MRC: can pad to zero b/c origin < 0: %s'%(str(nstart)))
		# 	data2 = np.zeros(nstart+grid.nxyz,dtype=data.dtype)
		# 	data2[nstart[0]:,nstart[1]:,nstart[2]:] = data
		# else:
		# 	data2 = data

		order = np.arange(3)
		if reverse:
			order = order[::-1]

		if reverse:
			mrcf.set_data(np.moveaxis(data,np.arange(3),order)) ## reverse for mrcfile wonkyness
		else:
			mrcf.set_data(data)

		for li,nsi,ni in zip(('x','y','z'),('nxstart','nystart','nzstart'),(0,1,2)):
			mrcf.header.origin[li] = origin[ni]
			mrcf.header.cella[li] = dxyz[ni] * nxyz[ni]
			mrcf.header.__setattr__(nsi,nstart[ni])

		# 	## required for (old) chimera
		# 	mrcf.header.nxstart = nstart[0]
		# 	mrcf.header.nystart = nstart[1]
		# 	mrcf.header.nzstart = nstart[2]

		# else:
		# 	mrcf.set_data(data)
		# 	for li,ni in zip(('x','y','z'),(0,1,2)):
		# 		mrcf.header.origin[li] = grid.origin[ni]
		# 		mrcf.header.cella[li] = grid.dxyz[ni] * grid.nxyz[ni] ## technically should be mxyz....
		# 	## required for (old) chimera
		# 	mrcf.header.nxstart = nstart[0]
		# 	mrcf.header.nystart = nstart[1]
		# 	mrcf.header.nzstart = nstart[2]
		# 	mrcf.set_data(np.moveaxis(data,np.arange(3),np.arange(3)[::-1])) ## reverse for mrcfile wonkyness
		mrcf.flush()

def print_mrc_header(mrcf):
	h = mrcf.header
	for k,v in zip(h.dtype.names,h.item()):
		print(k,':\t',v)
