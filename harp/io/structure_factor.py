import numpy as np
from ..density import gridclass

def load_sf(fname,options={}):
	import gemmi
	if not 'sample_rate' in options:
		options['sample_rate'] = 3.
	if not 'factor label' in options:
		options['factor label'] = 'F_meas_au'
	if not 'phase label' in options:
		options['phase label'] = 'phase_calc'

	# half_l = True
	# order = gemmi.HklOrient.HKL ## gemmi.HklOrient.LKH

	if fname.endswith('.ent') or fname.endswith('.ENT') or fname.endswith('.cif') or fname.endswith('.CIF') or fname.endswith('.ent.gz') or fname.endswith('.ENT.gz'):
		doc = gemmi.cif.read(fname)
		rblock = gemmi.as_refln_blocks(doc)[0]
		size = rblock.get_size_for_hkl(sample_rate=options['sample_rate'])
		if options['factor label'] in rblock.column_labels() and options['phase label'] in rblock.column_labels():
			_emap = rblock.transform_f_phi_to_map(options['factor label'],options['phase label'],size)
		else:

			raise Exception('Do not know how to load these structure factors\nAvailable options in file are: %s'%(rblock.column_labels()))

	elif fname.endswith('.mtz') or fname.endswith('.MTZ'):
		mtz = gemmi.read_mtz_file(fname)
		column_labels = np.array(mtz.column_labels())
		size = mtz.get_size_for_hkl(sample_rate=options['sample_rate'])
		if options['factor label'] in column_labels and options['phase label'] in column_labels:
			_emap = mtz.transform_f_phi_to_map(options['factor label'],options['phase label'],size)
		else:
			raise Exception('Do not know how to load these structure factors\nAvailable options in file are: %s'%(column_labels))
	else:
		raise Exception('Do not know how to load these structure factors')

	uc = _emap.unit_cell
	emap = np.array(_emap,copy=False)

	if np.any(np.array([uc.alpha,uc.beta,uc.gamma])!=90.0):

		## Do it the hard way....
		if options['minmax'] is None:
			raise Exception('Not a rectilinear grid.... Do more symmetry coding colin...')

			# # origin = np.zeros(3)
			# # xyzmin = uc.orthogonalize(gemmi.Fractional(0.,0.,0.))
			# # xyzmax = uc.orthogonalize(gemmi.Fractional(1.,1.,1.))
			# xyzmin = np.array((xyzmin.x,xyzmin.y,xyzmin.z))
			# xyzmax = np.array((xyzmax.x,xyzmax.y,xyzmax.z))
			# deltaxyz = (xyzmax-xyzmin)
			# nxyz = np.array(emap.shape)
			# # nxyz += nxyz//2
			# dxyz = deltaxyz/nxyz
			#
			# emin = emap.min()
			# emax = emap.max()
			# edelta = emin-emax
			# emin -= edelta*.05
			# emax += edelta*.05
			#
			# emap = np.zeros(nxyz)
			# for xi in range(nxyz[0]):
			# 	for yi in range(nxyz[1]):
			# 		for zi in range(nxyz[2]):
			# 			eval = _emap.interpolate_value(gemmi.Position(dxyz[0]*xi+xyzmin[0],dxyz[1]*yi+xyzmin[1],dxyz[2]*zi+xyzmin[2]))
			# 			if eval >= emin and eval <= emax:
			# 				emap[xi,yi,zi] = eval
			# grid = gridclass(origin,dxyz,nxyz)
		minmax = options['minmax']
		origin = minmax[0]
		deltaxyz = (minmax[1]-minmax[0])
		nxyz = np.array(emap.shape)
		# nxyz += nxyz//2
		dxyz = deltaxyz/nxyz

		emin = emap.min()
		emax = emap.max()
		edelta = emin-emax
		emin -= edelta*.05
		emax += edelta*.05

		emap = np.zeros(nxyz)
		for xi in range(nxyz[0]):
			for yi in range(nxyz[1]):
				for zi in range(nxyz[2]):
					pos = gemmi.Position(dxyz[0]*xi+minmax[0,0], dxyz[1]*yi+minmax[0,1], dxyz[2]*zi+minmax[0,2])
					eval = _emap.interpolate_value(pos)
					if eval >= emin and eval <= emax:
						emap[xi,yi,zi] = eval
		grid = gridclass(origin,dxyz,nxyz)

	else:
		origin = np.zeros(3)
		nxyz = np.array(emap.shape)
		cellxyz = np.array((uc.a,uc.b,uc.c))
		dxyz = cellxyz / nxyz
		grid = gridclass(origin,dxyz,nxyz)

	return grid,emap.astype('double')
