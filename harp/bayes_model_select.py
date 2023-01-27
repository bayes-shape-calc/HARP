'''
Notes:
======
- You call bms_molecule, which will loop over residues and run bms_residue on reach residue, which is kind of a wrapper to provide all the information into bms_xyz, which is the base case.
'''

import numpy as np
from . import density as dc
from . import models
from . import evidence
from .jit_math import erf
import time


def bms_residue(grid, data, subresidue, adfs, blobs, subgrid_size=8., sigmacutoff=5, offset=0.5, atom_types=None, atom_weights=None):
	'''
	Run the BMS calculation on a (sub)residue
		- adfs are \sigma (std dev) of atoms in M0 -- atomic DISTORTION factor... not exclusively DISPLACEMENT. e.g. includes bad reconstruction noise.
		- blobs are \sigma (std dev) of superatom in M1
		- removes hydrogen atoms before doing any calculation
		- results are extended to all atoms in the residue (including the removed hydrogens!)
		- moment matching a weight mixture of gaussians to a single gaussian gives that the super-atom \mu is located at the average \mu of the atoms. The sigma (std dev) of super atom is \sqrt(\Sigma) = \sqrt(<\Sigma> + var(\mu))
	'''

	nadf = adfs.size
	nblob = blobs.size
	ncalc = nadf+nblob

	## this includes the hydrogens not used in the calculation. make sure to keep the +,- down below!
	record_prob_good = np.zeros((subresidue.natoms))
	record_ln_ev = np.zeros((subresidue.natoms,ncalc))

	## if there are no atoms left in the residue, don't do the calculation
	if subresidue.natoms == 0:
		return record_prob_good,record_ln_ev

	com = subresidue.com()

	## only run the calculations using the local neighborhood (a la BITS definition)
	subgrid = dc.subgrid_center(grid,com,subgrid_size)
	subdata = dc.subgrid_extract(grid,data,subgrid).astype('double')

	## if the position is outside of the density grid, we can't calculate anything
	if subdata.size == 0:
		record_prob_good += 0. # model priors set to zero if outside.
		record_ln_ev -= 2*evidence.lnprior_factor_location + evidence.lnprior_factor_scale ## keep the + and - signs!
		return record_prob_good,record_ln_ev

	#### Initialize atom weights
	if atom_weights is None or atom_types is None:
		## usually normalized to carbon. Could be anything b/c a shape analysis, but needs to be something
		raise Exception('no weights provided!')

	weights = np.zeros(subresidue.natoms)
	for ati in range(atom_types.size):
		weights[subresidue.element == atom_types[ati]] = atom_weights[ati] ## assign weights
	if subresidue.occupancy.sum() > 0:
		## if multiple conformations, usually the occupancy is split, so do this. doesn't really matter unless trying to compare absolute values between residues, but that's non-sense anyway because the local neighborhoods might be different sizes.
		weights *= subresidue.occupancy

	#### Calculations
	if  models.version == 'C':
		ln_evidence_out = models._harp(subdata,com,subgrid.origin,subgrid.dxyz,subgrid.nxyz,subresidue.xyz,weights,adfs,blobs,sigmacutoff,offset)

	else:
		ln_evidence_out = np.zeros((ncalc))

		## Run the atomic model - adfs
		for i in range(nadf):
			atomic_model = models.density_atoms(subgrid, subresidue.xyz, weights, adfs[i], sigmacutoff, offset)
			ln_evidence_out[i] = evidence.ln_evidence(atomic_model, subdata)

		## Run the spheric  model - blobs
		superatom_weight = weights.sum()
		for i in range(nblob):
			spheric_model = models.density_point(subgrid, com, superatom_weight, blobs[i], sigmacutoff, offset)
			ln_evidence_out[nadf+i] = evidence.ln_evidence(spheric_model, subdata)

	## Model selection
	if np.all(~np.isfinite(ln_evidence_out)): ## something went wrong.
		record_prob_good += np.nan
	else:
		# ###### 1/K priors
		# ## both M0 and M1 are equal a priori at .5, therefore, e.g., individual blobs are .5/nblobs
		# lnpriors = np.zeros(ncalc) + 0.5
		# lnpriors[:nadf] /= float(nadf)
		# lnpriors[nadf:] /= float(nblob)

		###### max-ent mag priors
		lnpriors = np.zeros(ncalc) + .5
		if nadf > 1:
			dlnx = np.log(adfs[-1])-np.log(adfs[0])
			adfbounds = (adfs[1:]+adfs[:-1])/2.
			lnpriors[0] *= (np.log(adfbounds[0]) - np.log(adfs[0]))/dlnx
			lnpriors[nadf-1] *= (np.log(adfs[-1]) - np.log(adfbounds[-1]))/dlnx
			if nadf > 2:
				lnpriors[1:nadf-1] *= (np.log(adfbounds[1:]) - np.log(adfbounds[:-1]))/dlnx
		if nblob > 1:
			dlnx = np.log(blobs[-1])-np.log(blobs[0])
			blobbounds = (blobs[1:]+blobs[:-1])/2.
			lnpriors[nadf+0] *= (np.log(blobbounds[0]) - np.log(blobs[0]))/dlnx
			lnpriors[nadf+nblob-1] *= (np.log(blobs[-1]) - np.log(blobbounds[-1]))/dlnx
			if nblob > 2:
				lnpriors[nadf+1:nadf+nblob-1] *= (np.log(blobbounds[1:]) - np.log(blobbounds[:-1]))/dlnx

		lnpriors = np.log(lnpriors)
		R = ln_evidence_out+lnpriors
		R = np.exp(R-R.max()) ## remove the maximum to save the numerical precision -- it cancels in the prob calc.
		record_prob_good += np.nansum(R[:nadf])/np.nansum(R)
	record_ln_ev += ln_evidence_out[None,:] ## keep the plus sign!

	return record_prob_good, record_ln_ev

def bms_residue_reduced(grid, data, subresidue, adf =.25, subgrid_size=8.,sigmacutoff=5,offset=0.5,atom_types=None,atom_weights=None):
	'''
	Run the reduced version of the BMS calculation on a (sub)residue
		- only one M0 and one M1 are calculated!!!
		- adf is the \sigma (std dev) of atoms in M0 -- atomic DISTORTION factor... not exclusively DISPLACEMENT. e.g. includes bad reconstruction noise.
		- blob is calculated as sqrt(adf^2+var(\mu))
		- removes hydrogen atoms before doing any calculation
		- results are extended to all atoms in the residue (including the removed hydrogens!)
	'''

	## if there are no atoms left in the residue, don't do the calculation
	if subresidue.natoms == 0:
		record_prob_good = np.zeros((subresidue.natoms))
		record_ln_ev = np.zeros((subresidue.natoms,2))
		return record_prob_good,record_ln_ev

	# r = np.sqrt(np.sum(subresidue.xyz**2.,axis=1))
	r = np.sqrt(np.sum((subresidue.xyz-subresidue.com()[None,:])**2.,axis=1))
	E_r = np.mean(r)
	E_rr = np.mean(r*r)
	blob = np.sqrt(adf**2.+E_rr-E_r**2.)
	return bms_residue(grid, data, subresidue, np.array((adf,)), np.array((blob,)), subgrid_size, sigmacutoff, offset,atom_types,atom_weights)


def optimal_global_adf(grid,data,mol,sigmacutoff,offset,sigmas = np.linspace(0.3,1.2,20)):
	mol_nohet = mol.remove_hetatoms()
	weights = np.ones(mol_nohet.natoms)

	ln_ev = np.zeros_like(sigmas)
	for i in range(sigmas.size):
		model = models.density_atoms(grid, mol_nohet.xyz, weights, sigmas[i], sigmacutoff, offset)
		ln_ev[i] = evidence.ln_evidence(model, data)
	out = sigmas[np.nanargmax(ln_ev)]
	return out


def gen_adfsblobs(adf_low=None,adf_high=None,adf_n=None,blob_low=None,blob_high=None,blob_n=None):
	if adf_low is None or adf_high is None or adf_n is None:
		adf_low = .25
		adf_high = 2.5
		adf_n = 10
	if blob_low is None or blob_high is None or blob_n is None:
		blob_low = .25
		blob_high = 8.
		blob_n = 20

	adfs = np.logspace(np.log10(adf_low),np.log10(adf_high),adf_n)
	blobs = np.logspace(np.log10(blob_low),np.log10(blob_high),blob_n)
	return adfs,blobs


def bms_molecule(grid, data, mol, adfs = None, blobs = None, subgrid_size=8., sigmacutoff=5, offset=.5, emit=print, chains=None, reduced=False, atom_types=None, atom_weights=None):
	'''
	Run the BMS calculation on an entire molecule by looping over all of the residues in all of the chains
		- results are per residue but extended for all atoms
	'''

	if adfs is None:
		adfs = gen_adfsblobs()[0]
	if blobs is None:
		blobs = gen_adfsblobs()[1]

	# emit(adfs)
	# emit(blobs)
	emit('N_adfs = %d, N_blobs = %d'%(adfs.size,blobs.size))

	if atom_types is None or atom_weights is None:
		atom_weights = np.array([.05,1.,1.,1.,2.,2.])
		atom_types = np.array(['H','C','N','O','P','S'])
		emit('Using default weights: %s'%(str(atom_weights)))

	## initialize results
	nadf = adfs.size
	nblob = blobs.size
	if reduced:
		ncalc = 2
		adf = optimal_global_adf(grid,data,mol,sigmacutoff,offset)
		emit('Found Global ADF: %.2f'%(adf))
	else:
		ncalc = nadf+nblob

	record_prob_good = np.zeros((mol.natoms)) ## Natoms
	record_ln_ev = np.zeros((mol.natoms,ncalc)) ## Natoms x Nmodels (atomic, then blobs in sigmas_blob)

	## Loop over chains`
	for chain in mol.unique_chains:
		if not chains is None: ## Skip this chain?
			if not chain in chains:
				continue

		## Extract current chain
		subchain = mol.get_chain(chain)
		if np.all(subchain.hetatom):
			continue

		## Loop over residues in current chain
		t0 = time.time()
		avgprob = []
		for i in range(subchain.unique_residues.size):
			## Extract current residue
			resi = subchain.unique_residues[i]
			subresidue = subchain.get_residue(resi)

			## Run Calculation
			if reduced:
				out = bms_residue_reduced(grid, data, subresidue, adf, subgrid_size ,sigmacutoff, offset,atom_types,atom_weights)
			else:
				out = bms_residue(grid, data, subresidue, adfs, blobs, subgrid_size, sigmacutoff, offset,atom_types,atom_weights)

			## Extend results to all atoms in residue
			avgprob.append(out[0][0])
			which_atoms = np.isin(mol.atomid,subresidue.atomid)
			record_prob_good[which_atoms] = out[0]
			record_ln_ev[which_atoms] = out[1]

		## Reporting on per chain timing performance
		t1 = time.time()
		nres = subchain.unique_residues.size
		emit('Chain: %s, %d residues, <P_res> = %.4f, t = %.2f sec, <t/res> = %.2f msec'%(chain,nres,np.nanmean(avgprob),t1-t0,(t1-t0)/float(nres)*1000.))

	return record_prob_good,record_ln_ev
