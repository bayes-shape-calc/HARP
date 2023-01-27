import numpy as np

def save_ccp4(fname,grid,data):
	import gemmi
	ccp4 = gemmi.Ccp4Map()
	g = gemmi.FloatGrid(*grid.nxyz)
	uc = (grid.nxyz*grid.dxyz).tolist() + [90.,90.,90.]
	g.set_unit_cell(gemmi.UnitCell(*uc))
	v = np.array(g,copy=False)
	v[:,:,:] = np.asfortranarray(data) ## there might be a C/Fortran issue here
	ccp4.grid = g
	ccp4.update_ccp4_header(2, True)
	ccp4.write_ccp4_map(fname)
