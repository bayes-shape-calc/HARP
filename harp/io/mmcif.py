import numpy as np
import gzip
from ..molecule import atomcollection

def pos_read(f):
	pos = f.tell()
	line = f.readline().decode('utf-8')
	line = line.lstrip(' ')
	line = line.lstrip('\t')
	return pos,line

def _get_full_loop(f,pos,line):
	tag_pos,loopline = pos_read(f)
	loopline = loopline
	tags = []
	prefix = loopline.split('.')[0]
	while loopline[0] == '_':
		tags.append(loopline.strip())
		tag_pos,loopline = pos_read(f)

	entries = []
	entry_pos = tag_pos
	while loopline:
		if (loopline == '\n') or loopline.startswith('#') or loopline.startswith('_') or loopline.startswith('loop_') or loopline.startswith('block_') or loopline.startswith('data_') or loopline.startswith('save_'):
			f.seek(entry_pos)
			break

		s = _special_split(f,loopline,tags)
		entries.append(s)

		## keep going to the next line... maybe
		entry_pos,loopline = pos_read(f)

	return prefix,tags,entries

def _special_split(f,loopline,tags):
	delimiters = ['\'','\"']
	skips = [' ','\t','\n']
	s = []
	recording = False
	delimiter = None
	current = ''
	while len(s) < len(tags):
		if loopline[0] == ';': ## go until the next line that startswith a semicolon...
			current += loopline[1:]
			while 1:
				pos,loopline = pos_read(f)
				current += loopline
				if loopline[0] == ';':
					s.append(current)
					current = ''
					break
				elif not loopline:
					raise Exception('Reached EOF trying to find ; in:\n%s'%(value))
		else:
			for l in loopline:
				if recording:
					if delimiter is None:
						if l in skips:
							s.append(current)
							current = ''
							recording = False
						else:
							current = current + l

					else:
						current = current + l
						if l == delimiter:
							s.append(current)
							current = ''
							recording = False
							delimiter = None

				elif not recording:
					if not l in skips:
						recording = True
						current = l
						if l in delimiters:
							delimiter = l

		if len(s) >= len(tags):
			break
		else:
			pos,loopline = pos_read(f)
	return s

def _get_single_entry(f,pos,line):
	s = line.split()
	entry = s[0]

	if len(s) > 1: ## if there was anything after the entry
		value = ' '.join(s[1:]) ## recombine the line
		success,value = _complete_delimiters(f,value)

	else: ## if nothing else on this line... look to see if the next is a "string..."
		pos,line = pos_read(f)
		value = '' + line
		if line[0] == ';': ## go until the next ;...
			while 1:
				pos,line = pos_read(f)
				value += line
				if line[0] == ';':
					break
				elif not line:
					raise Exception('Reached EOF trying to find ; in:\n%s'%(value))
		else:
			f.seek(pos) ## rewind...
	value = value.replace('\n','')
	return entry,value

def _complete_delimiters(f,value):
	## True if completed delimiters, False if no delimiter, value is string
	delimiters = ['\'','\"']
	delimiter = value[0]
	if delimiter in delimiters:
		while value.count(delimiter)%2 == 1:
			_pos,line = pos_read(f)
			if not line:
				raise Exception('Reached EOF trying to find %s in:\n%s'%(delimiter,value))
			value += line
		return True,value
	return False,value


def _load_mmcif_dict(f,exclusive=None):
	## mostly fixed now
	outs = []
	out = None
	pos,line = pos_read(f)

	armed = True
	if not exclusive is None:
		armed = False

	while line: 
		line = line[:-1]

		if line.startswith('data_'): ## it's a block
			if not out is None:
				outs.append(out)
			# print(repr(line))
			out = {}

		elif line.startswith('loop_') and not out is None: ## pull out a loop
			if not exclusive is None: ## peek at the next line
				pos,peek = pos_read(f)
				f.seek(pos)
				if peek.startswith(exclusive):
					armed = True
			if armed:
				prefix,tags,entries = _get_full_loop(f,pos,line)
				out[prefix] = [tags,entries]
				if not exclusive is None:
					armed = False
					break

		elif line.startswith('_') and not out is None: ## pull out a single entry
			if not exclusive is None:
				if line.startswith(exclusive):
					armed = True
			if armed:
				entry,value = _get_single_entry(f,pos,line)
				out[entry] = value
				if not exclusive is None:
					armed = False
					break

		pos,line = pos_read(f)

	outs.append(out)
	return outs

def _peek_headers_loop(f,lname):
	pos,line = pos_read(f)
	armed = False
	headers = []
	while line:
		if not armed:
			if line.startswith(lname):
				headers.append(line.strip())
				armed = True
		else:
			if line.startswith(lname):
				headers.append(line.strip())
			else:
				break
		pos,line = pos_read(f)
	return headers

def peek_headers_loop(fname,lname):
	if fname.endswith('.gz'):
		with gzip.open(fname, 'rb') as f:
			headers = _peek_headers_loop(f,lname)
	else:
		with open(fname,'rb') as f:
			headers = _peek_headers_loop(f,lname)
	return headers

def load_mmcif_dict(fname,exclusive=None):
	if fname.endswith('.gz'):
		with gzip.open(fname, 'rb') as f:
			dictionaries = _load_mmcif_dict(f,exclusive)
	else:
		with open(fname,'rb') as f:
			dictionaries = _load_mmcif_dict(f,exclusive)
	# if len(dictionaries) == 1:
		# return dictionaries[0]
	return dictionaries

def _find_emdb(f):
	stopper = '_atom_site.'
	for line in f:
		line = line.decode('utf-8')
		if line.startswith(stopper): ## ya went too far
			return None
		s = line.split()
		for ss in s:
			if ss.lower().startswith('emd-'):
				sss = ss[4:]
				return ''.join([sss for sss in ss[4:] if sss.isnumeric()])
				# return ss[4:]

def find_emdb(fname):
	if fname.endswith('.gz'):
		with gzip.open(fname, 'rb') as f:
			emdb = _find_emdb(f)
	else:
		with open(fname,'rb') as f:
			emdb = _find_emdb(f)
	return emdb

def gemmi_open(fname):
	import gemmi
	doc = gemmi.cif.read(fname)
	block = doc.sole_block()
	return block
def gemmi_kw_quick(block,kw):
	''' '_cell.length_a' '''
	return block.find_pair(kw)
def gemmi_loop_quick(block,kw):
	''' '_atom_type.symbol' '''
	return np.array(block.find_loop(kw))

def check_loop_header(fname,loopname,look_list):
	opener = open
	if fname.endswith('.gz'):
		opener = gzip.open
	with opener(fname,'r') as f:
		pos,line = pos_read(f)
		armed = False
		headers = [False,]*len(look_list)
		while line:
			if line.startswith(loopname):
				if not armed:
					armed = True
				header = line[len(loopname):]
				for i in range(len(look_list)):
					if header.startswith(look_list[i]):
						headers[i] = True
			else:
				if armed:
					break
			pos,line = pos_read(f)
	return all(headers)
	
################################################################################
################################################################################
################################################################################

def load_mmcif(fname, only_polymers=False, first_model=True,auth=False):
	'''
	into an atomcollection class
	'''

	tags,entries = load_mmcif_dict(fname,exclusive='_atom_site.')[0]['_atom_site']
	tags = np.array(tags)
	entries = np.array(entries).T

	atomid = entries[np.nonzero(tags=='_atom_site.id')].astype('int')[0]
	element = entries[np.nonzero(tags=='_atom_site.type_symbol')][0]
	atomname = entries[np.nonzero(tags=='_atom_site.label_atom_id')][0]
	resname = entries[np.nonzero(tags=='_atom_site.label_comp_id')][0]
	if auth:
		chain = entries[np.nonzero(tags=='_atom_site.auth_asym_id')][0]
	else:
		chain = entries[np.nonzero(tags=='_atom_site.label_asym_id')][0]
	conf = entries[np.nonzero(tags=='_atom_site.label_alt_id')][0]
	occupancy = entries[np.nonzero(tags=='_atom_site.occupancy')].astype('double')[0]
	bfactor = entries[np.nonzero(tags=='_atom_site.B_iso_or_equiv')].astype('double')[0]
	hetatom = entries[np.nonzero(tags=='_atom_site.group_PDB')][0] != 'ATOM'
	modelnum = entries[np.nonzero(tags=='_atom_site.pdbx_PDB_model_num')][0]
	xyz = np.array([entries[np.nonzero(tags=='_atom_site.Cartn_%s'%(xyz))][0] for xyz in ['x','y','z']]).astype('double').T

	######################################## resid... hard to cast to int.
	authresid = entries[np.nonzero(tags=='_atom_site.auth_seq_id')][0]
	bad = np.bitwise_or(authresid=='.',authresid=='?')
	authresid[bad] = np.array((-1),dtype=authresid.dtype)
	authresid = authresid.astype('int')

	labelresid = entries[np.nonzero(tags=='_atom_site.label_seq_id')][0]
	bad = np.bitwise_or(labelresid=='.',labelresid=='?')
	labelresid[bad] = np.array((-1),dtype=labelresid.dtype)
	labelresid = labelresid.astype('int')
	
	if auth:
		resid = authresid
	else:
		resid = labelresid

	# mol = atomcollection(atomid,labelresid,resname,atomname,chain,element,conf,xyz,occupancy,bfactor,hetatom,modelnum,authresid)
	mol = atomcollection(atomid,resid,resname,atomname,chain,element,conf,xyz,occupancy,bfactor,hetatom,modelnum,labelresid)
	
	if only_polymers: ### often PTMd residues are Hetatoms so... you can't remove bad stuff that way.
		try:
			## Pull out the _entity entries corresponding to polymers. store as strings in keep_entities
			method = 0
			try:
				tags2, entries2 = load_mmcif_dict(fname,'_entity.')[0]['_entity']
				tags2 = np.array(tags2)
				entries2 = np.array(entries2).T
				method = 1
			except:
				try:
					entity_id = load_mmcif_dict(fname,'_entity.id')[0]['_entity.id']
					entity_type = load_mmcif_dict(fname,'_entity.type')[0]['_entity.type']
					tags2 = np.array(['_entity.id','_entity.type'])
					entries2 = np.array([[entity_id,entity_type],]).T
					method = 2
				except:
					pass
			if method == 0:
				raise Exception('Cannot get _entity entry in %s!'%(fname))
			
			_entities = entries2[np.nonzero(tags2=='_entity.id')][0]
			_polymers = entries2[np.nonzero(tags2=='_entity.type')][0]
			keep_entities = _entities[np.nonzero(_polymers == 'polymer')]
			entities = entries[np.nonzero(tags=='_atom_site.label_entity_id')][0]
			
			mol = mol.get_set(np.isin(entities,keep_entities))

		except:
			## if everything fails, just do it the old fashioned way. This removes things like PTMs though.
			mol = mol.remove_hetatoms()
	
	if first_model: ## sometimes several are built into the same density; just grab the first model: _atom_site.pdbx_PDB_model_num = 1
		unique_models = np.unique(mol.modelnum)
		if unique_models.size > 1:
			mol = mol.get_set(mol.modelnum == unique_models[0])

	return mol

def write_bfactor(file_cif,mol,prob_good):
	out = ''
	with open(file_cif,'rb') as f:
		pos,line = pos_read(f)

		while line:
			if  line.startswith('loop_'):
				pos,peek = pos_read(f)
				f.seek(pos)
				if peek.startswith('_atom_site.'):
					prefix,tags,entries = _get_full_loop(f,pos,line)
								
					out +='loop_\n'
					for tag in tags:
						out += tag+'\n'
					b_ind = tags.index('_atom_site.B_iso_or_equiv')
					resid_ind = tags.index('_atom_site.auth_seq_id')
					chain_ind = tags.index('_atom_site.auth_asym_id')
				
					chains0 = np.array([e[chain_ind] for e in entries])
					resid0 = np.array([e[resid_ind] for e in entries])
					bad = np.bitwise_or(resid0=='.',resid0=='?')
					resid0[bad] = np.array((-1),dtype=resid0.dtype)
					resid0 = resid0.astype('int')
					
				
					for i in range(len(entries)):
						entries[i][b_ind] = '0.0000'
					# lookup = {}
					# for chain in np.unique(chains0):
					# 	lookup[chain] = {}
					# 	for resid in np.unique(resid0):
					# 		lookup[chain][resid] = 0.
					# mol.bfactor = prob_good.copy()
					# for chain in mol.unique_chains:
					# 	subchain = mol.get_chain(chain)
					# 	for residue in subchain.unique_residues:
					# 		subresidue = mol.get_residue(residue)
					# 		lookup[chain][residue] = subresidue.bfactor[0]
					# for i in range(chains0.size):
					# 	if chains0[i] in lookup.keys():
					# 		if resid0[i] in lookup[chains0[i]]:
					# 			bf = lookup[chains0[i]][resid0[i]]
					# 			entries[i][b_ind] = '%.4f'%((bf)*100.)
							
					
					keep = np.bitwise_and(np.isin(chains0,mol.chain),np.isin(resid0,mol.resid))
					for i in range(chains0.size):
						if keep[i]:
							kind =np.argmax(np.bitwise_and(np.array([mci==chains0[i] for mci in mol.chain]),np.array([mri==resid0[i] for mri in mol.resid])))
							# kind=np.argmax(np.bitwise_and(np.isin(chains0[i],mol.chain),np.isin(resid0[i],mol.resid)))

							entries[i][b_ind] = '%.4f'%((1-prob_good[kind])*100.)
						else:
							entries[i][b_ind] = '0.0000'
				
					for entry in entries:
						out += ' '.join(entry) + '\n'

				else:
					out += line
			else:
				out += line
			pos,line = pos_read(f)

	import os
	fout = os.path.splitext(file_cif)[0]+'_HARP.cif'
	print('wrote bfactor yes')
	print(fout)
	with open(fout,'w') as f:
		f.write(out)