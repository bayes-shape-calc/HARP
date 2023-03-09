#!/usr/bin/env python

import os
import argparse
import numpy as np
import sys
from .. import models
from .. import io
from .. import molecule
from .. import density
from .. import bayes_model_select
try:
	import blobview
	flag_blobview = True
except:
	flag_blobview = False

def dir_path(string):
	if os.path.isdir(string):
		return string
	else:
		err = "\n\nError: output directory \"%s\" is not a directory that exists yet. Make it?\n"%(string)
		raise NotADirectoryError(err)

def file_path(string):
	if os.path.exists(string):
		return string
	else:
		err = "\n\nError: file %s does not exist\n\n"%(string)
		raise Exception(err)

def harpcalc(pdbid, basedir, label_sf='pdbx_FWT', label_phase='pdbx_PHWT', adfs=None, blobs=None, offset=.5, chains=None, verbose=True, quiet=False, emit=print, overwrite=False, end_view=False, view_threshold=.5, input_files=None, skip_calc=False, skip_load=False, reduced=False, only_polymers=True):
	if verbose: emit('Using %s library'%(models.version))

	## Try to download files from the wwPDB
	success, path_mol, path_density, flag_density = io.rcsb.get_pdb(pdbid,basedir,overwrite,verbose,emit)

	## Load molecule
	if not skip_load:
		if not quiet: emit('Loading %s'%(pdbid))
		mol = molecule.load(path_mol,only_polymers)

		## Load X-ray density from structure factors
		if flag_density == 'xray':
			try:
				pad = 8.
				options = {
					'minmax':np.array((mol.xyz.min(0)-pad,mol.xyz.max(0)+pad)),
					'factor label':label_sf,
					'phase label':label_phase,
				}
				## gemmi library handles .gz too
				grid,dens = density.load(path_density,options)
				if not np.all([density.xyz_in_grid(grid,mol.xyz[i]) for i in range(mol.natoms)]):
					if verbose: emit('Not all atoms fit in the density. Often an X ray problem. Wrapping density using PBC to fit them all')
					grid,dens = density.trim_density_to_mol(grid,dens,mol)
			except Exception as e:
				emit('%s\nCould not load structure factors -- check labels'%(e))
				return

		## Load EM density
		elif flag_density == 'em':
			## mrcfile library  handles .gz too
			grid,dens = density.load(path_density)

		## Finished loading
		if verbose: emit('--> Loaded %s - %s'%(pdbid,flag_density))
		if verbose: emit('--> Grid: %s, %s, %s'%(str(grid.origin),str(grid.nxyz),str(grid.dxyz)))

	## Perform the HARP model selection on all residues in all chains of molecule.
	probs = np.zeros(mol.xyz.shape[0])
	if not skip_calc:
		probs, ln_evidence = bayes_model_select.bms_molecule(grid, dens, mol, chains = chains, adfs=adfs , blobs=blobs, offset=offset, reduced=reduced)
		keep = np.bitwise_or(mol.atomname=='\"C4\'\"',mol.atomname=='CA')
		last_resid = -1
		last_chain = ''
		for i in range(keep.size):
			if keep[i]:
				if last_resid == mol.resid[i] and last_chain == mol.chain[i]:
					keep[i] = False
				else:
					last_resid = mol.resid[i]
					last_chain = mol.chain[i]


		## Return a little information
		if not quiet:
			emit("<P_res>: %.4f"%(np.nanmean(probs[keep])))
			emit(" Finite: %d / %d"%(np.isfinite(probs[keep]).sum(),keep.sum()))

		## Write info to a CSV file
		path_out = os.path.join(basedir,'results_%s.csv'%(pdbid))
		if not os.path.exists(path_out) or overwrite:
			with open(path_out,'w') as f:
				# f.write('Chain,Residue ID,Auth ID,Residue Name,P_res,sigma_ADF MAP,sigma_blob MAP\n')
				f.write('Chain,Residue ID,Auth ID,Residue Name,P_res\n')
				for i in range(keep.sum()):
					# out = '%s,%s,%s,%s,%.5f,%.5f,%.5f\n'%(mol.chain[keep][i],mol.resid[keep][i],mol.authresid[keep][i],mol.resname[keep][i],probs[keep][i], adfs[np.nanargmax(ln_evidence[keep][i][:adfs.size])], blobs[np.nanargmax(ln_evidence[keep][i][adfs.size:])])
					out = '%s,%s,%s,%s,%.5f\n'%(mol.chain[keep][i],mol.resid[keep][i],mol.authresid[keep][i],mol.resname[keep][i],probs[keep][i])
					f.write(out)
			if verbose: emit('Results for %s saved in %s'%(pdbid,path_out))

	## View the molecule and density -- only if blobview is installed
	if end_view and flag_blobview:
		if verbose: emit('Launching Molecule Viewer')
		options = blobview.default_options
		dens -= dens.mean()
		dens /= dens.max()
		options['iso_thresh'] = view_threshold
		# options['iso_thresh'] = density.mean()+density.std()
		blobview.view(dens,grid,mol,probs,options=options,verbose=verbose)

def main():
	parser = argparse.ArgumentParser(description="Use HARP to calculate trustiworthness of model")

	parser.add_argument('pdb', type=str, help='The ID of the structure in the wwPDB to process')
	parser.add_argument('-o','--output',type=dir_path,default='./',help='Directory to store results and files')
	parser.add_argument('--overwrite',action='store_true',default=False,help='Overwrite old data files or not')

	parser.add_argument('--mmcif_location',type=file_path,default=None,help='Location of mmcif file. Must include density_location too')
	parser.add_argument('--density_location',type=file_path,default=None,help='Location of density file. Must include mmcif location too')
	parser.add_argument('--label_sf',type=str,default='pdbx_FWT',help='Structure factor label in X-ray file')
	parser.add_argument('--label_phase',type=str,default='pdbx_PHWT',help='Phase label in X-ray file')

	########### Working on implementing
	# parser.add_argument('--atom_width',type=float,default=.6,help='Width of atoms in density models (standard deviation of 3D normal distribution)')
	parser.add_argument('--atoms_min',type=float,default=.25,help='Minimum width of atoms to check (standard deviation of 3D normal distribution)')
	parser.add_argument('--atoms_max',type=float,default=1.0,help='Maximum width of atoms to check (standard deviation of 3D normal distribution)')
	parser.add_argument('--atoms_num',type=int,default=10,help='Number of atoms models to use. Spaced evenly between atom_min and atom_max')

	parser.add_argument('--voxel_offset',type=float,default=.5,help='XYZ coordinate offset for each density voxel to relate to atomic model coordinates. A value of 0.0 is edge centered, 0.5 is face centered.')

	parser.add_argument('--blobs_min',type=float,default=.25,help='Minimum width of a blob to check (standard deviation of 3D normal distribution)')
	parser.add_argument('--blobs_max',type=float,default=3.0,help='Maximum width of a blob to check (standard deviation of 3D normal distribution)')
	parser.add_argument('--blobs_num',type=int,default=20,help='Number of blobs to use. Spaced evenly between blob_min and blob_max')

	parser.add_argument('--skip_calc',action='store_true',default=False,help='Skip the trustworthiness calculation, still download and visualize')
	parser.add_argument('--skip_load',action='store_true',default=False,help='Skip loading molecule and density')
	parser.add_argument('--only_polymers',action='store_true',default=True,help='Only calculate chains that come from entities that are polymers')

	if flag_blobview:
		## you have to have blobview
		parser.add_argument('--view_threshold',type=float,default=.2,help='View model - isosurface threshold value. Mean is 0, Max is 1.')
		parser.add_argument('--view',action='store_true',default=False,help='View model after calculation is completed')

	group_loud = parser.add_mutually_exclusive_group()
	group_loud.add_argument("--verbose", action="store_true")
	group_loud.add_argument("--quiet", action="store_true")
	group_loud.add_argument("--normal", action="store_true",default=True)

	group_fullreduced = parser.add_mutually_exclusive_group()
	group_fullreduced.add_argument("--reduced", action="store_true")
	group_fullreduced.add_argument("--full", action="store_true",default=True)

	group_lib = parser.add_mutually_exclusive_group()
	group_lib.add_argument('--use_c',action='store_true',default=False,help='Use the C library to render models')
	group_lib.add_argument('--use_python',action='store_true',default=True,help='Use Python to render models')

	parser.add_argument('-c','--chains', nargs='+', help='Chains in mmcif to run the calculation for. Separate by spaces. If not provided, all chains are used.')

	args = parser.parse_args()

	if args.use_c:
		models.use_c()
	elif args.use_python:
		models.use_python()

	if flag_blobview:
		aview = args.view
		aviewt = args.view_threshold
	else:
		aview = None
		aviewt = None


	adfs,blobs = bayes_model_select.gen_adfsblobs(args.atoms_min,args.atoms_max,args.atoms_num,args.blobs_min,args.blobs_max,args.blobs_num)


	harpcalc(
		pdbid = args.pdb,
		basedir = args.output,
		verbose = args.verbose,
		quiet = args.quiet,
		overwrite = args.overwrite,
		end_view = aview,
		view_threshold = aviewt,
		input_files = [args.mmcif_location, args.density_location],
		label_sf = args.label_sf,
		label_phase = args.label_phase,
		adfs = adfs,
		blobs = blobs,
		offset = args.voxel_offset,
		chains = args.chains,
		skip_calc=args.skip_calc,
		skip_load=args.skip_load,
		reduced=args.reduced,
		only_polymers=args.only_polymers,
	)

if __name__ == '__main__':
	main()
