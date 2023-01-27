'''
RCSB PDB
============================================================
	path_: help you get the filenames in the RCSB directory structure
	low-level: ftp,rcsb
	medium-level: mmcif,structure_factors
	high-level: xray,em
'''

import os

################################################################################
## These are how to get filenames within the RCSB structure

def path_cif_local(pdbid,fdir):
	# twocode = pdbid[1:3].lower()
	# local_path = os.path.join(fdir,twocode)
	# if not os.path.exists(local_path):
	# 	os.mkdir(local_path)
	# local_path = os.path.join(local_path,'%s.cif.gz'%(pdbid.lower()))
	local_path = os.path.join(fdir,'%s.cif.gz'%(pdbid.lower()))
	return local_path

def path_cif_ftp(pdbid,fdir):
	twocode = pdbid[1:3].lower()
	rcsb_ftp_mmcif = 'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF'
	ftp_path = os.path.join(rcsb_ftp_mmcif,twocode)
	ftp_path = os.path.join(ftp_path,'%s.cif.gz'%(pdbid.lower()))
	return ftp_path

def path_sf_local(pdbid,fdir):
	# twocode = pdbid[1:3].lower()
	fname = 'r%ssf.ent.gz'%(pdbid.lower())
	# local_path = os.path.join(fdir,twocode)
	# if not os.path.exists(local_path):
		# os.mkdir(local_path)
	local_sf = os.path.join(fdir,fname)
	return local_sf

def path_sf_ftp(pdbid,fdir):
	twocode = pdbid[1:3].lower()
	fname = 'r%ssf.ent.gz'%(pdbid.lower())
	rcsb_ftp_mmcif = 'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/structure_factors'
	ftp_path = os.path.join(rcsb_ftp_mmcif,twocode)
	ftp_path = os.path.join(ftp_path,fname)
	return ftp_path

def path_map_ftp(emdbid,fdir):
	fname = 'emd_%s.map.gz'%(emdbid)
	ftp_map = 'ftp://ftp.wwpdb.org/pub/emdb/structures/EMD-%s/map/%s'%(emdbid,fname)
	return ftp_map
def path_map_local(emdbid,fdir):
	fname = 'emd_%s.map.gz'%(emdbid)
	local_map = os.path.join(fdir,fname)
	return local_map

################################################################################

def checkpath(fdir):
	if not os.path.exists(fdir):
		os.mkdir(fdir)

def ftp(ftp_path,local_path,dry=False,verbose=True,emit=print):
	import shutil
	import urllib.request as request
	from contextlib import closing

	if not dry:
		with closing(request.urlopen(ftp_path,timeout=1)) as r:
			with open(local_path, 'wb') as f:
				if verbose: emit('FTP: %s --> %s'%(ftp_path,local_path))
				shutil.copyfileobj(r, f)
	else:
		if verbose: emit('Dry Run: %s --> %s'%(ftp_path,local_path))

def download_wrapper(fname_ftp,fname_local,overwrite=False,verbose=True,emit=print):
	if not os.path.exists(fname_local) or overwrite:
		try:
			ftp(fname_ftp,fname_local,verbose=verbose,emit=emit)
			return
		except:
			if verbose: emit('Failed ftp: %s, %s'%(fname_ftp,fname_local))
	else:
		if verbose: emit('Already exists: %s'%(fname_local))

def download_mmcif(pdbid,fdir,overwrite=False,verbose=True,emit=print):
	checkpath(fdir)
	ftp_cif = path_cif_ftp(pdbid,fdir)
	local_cif = path_cif_local(pdbid,fdir)
	download_wrapper(ftp_cif,local_cif,overwrite,verbose,emit)
	return local_cif


def download_structure_factors(pdbid,fdir,overwrite=False,verbose=True,emit=print):
	checkpath(fdir)
	ftp_sf = path_sf_ftp(pdbid,fdir)
	local_sf = path_sf_local(pdbid,fdir)
	download_wrapper(ftp_sf,local_sf,overwrite,verbose,emit)
	return local_sf

def download_em_map(pdbid,fdir,overwrite=False,verbose=True,emit=print):
	checkpath(fdir)
	from .mmcif import find_emdb
	local_cif = path_cif_local(pdbid,fdir)
	emdbid = find_emdb(local_cif)
	if (not emdbid is None) and (not emdbid == ''):
		ftp_map = path_map_ftp(emdbid,fdir)
		local_map = path_map_local(emdbid,fdir)
		download_wrapper(ftp_map,local_map,overwrite,verbose,emit)
		return local_map
	else:
		if verbose: emit("Could not identify corresponding EMDB entry to %s"%(pdbid))
	return None

def check_map_type(local_cif):
	from .mmcif import load_mmcif_dict
	d = load_mmcif_dict(local_cif,'_exptl.method') ## only bother with _exptl.method
	if d['_exptl.method'] == "\'X-RAY DIFFRACTION\'":
		return 'xray'
	return 'em'


def get_pdb(pdbid,fdir,overwrite=False,verbose=True,emit=print):
	checkpath(fdir)
	local_cif = download_mmcif(pdbid,fdir,overwrite,verbose,emit)

	flag = check_map_type(local_cif)
	if flag == 'xray':
		local_density = download_structure_factors(pdbid,fdir,overwrite,verbose,emit)
	else:
		local_density = download_em_map(pdbid,fdir,overwrite,verbose,emit)

	if os.path.exists(local_cif) and os.path.exists(local_density):
		return True,local_cif,local_density,flag
	if verbose: emit('Error: could not load files in \"%s\" for PDB ID: %s'%(fdir,pdbid))
	return False,None,None,None
