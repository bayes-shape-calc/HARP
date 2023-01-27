import numpy as np
import gzip
from ..molecule import atomcollection

def pos_read(f):
	pos = f.tell()
	line = f.readline().decode('utf-8')
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
	if len(dictionaries) == 1:
		return dictionaries[0]
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

def load_mmcif(fname):
	'''
	into an atomcollection class
	'''

	tags,entries = load_mmcif_dict(fname,exclusive='_atom_site.')['_atom_site']
	tags = np.array(tags)
	entries = np.array(entries).T

	atomid = entries[np.nonzero(tags=='_atom_site.id')].astype('int')[0]
	element = entries[np.nonzero(tags=='_atom_site.type_symbol')][0]
	atomname = entries[np.nonzero(tags=='_atom_site.label_atom_id')][0]
	resname = entries[np.nonzero(tags=='_atom_site.label_comp_id')][0]
	chain = entries[np.nonzero(tags=='_atom_site.label_asym_id')][0]
	conf = entries[np.nonzero(tags=='_atom_site.label_alt_id')][0]
	occupancy = entries[np.nonzero(tags=='_atom_site.occupancy')].astype('double')[0]
	bfactor = entries[np.nonzero(tags=='_atom_site.B_iso_or_equiv')].astype('double')[0]
	hetatom = entries[np.nonzero(tags=='_atom_site.group_PDB')][0] != 'ATOM'

	xyz = np.array([entries[np.nonzero(tags=='_atom_site.Cartn_%s'%(xyz))][0] for xyz in ['x','y','z']]).astype('double').T

	######################################## resid... hard to cast to int.
	authresid = entries[np.nonzero(tags=='_atom_site.auth_seq_id')][0]
	bad = np.bitwise_or(authresid=='.',authresid=='?')
	authresid[bad] = np.array((-1),dtype=authresid.dtype)

	labelresid = entries[np.nonzero(tags=='_atom_site.label_seq_id')][0]
	bad = np.bitwise_or(labelresid=='.',labelresid=='?')
	labelresid[bad] = np.array((-1),dtype=labelresid.dtype)
	labelresid = labelresid.astype('int')

	mol = atomcollection(atomid,labelresid,resname,atomname,chain,element,conf,xyz,occupancy,bfactor,hetatom,authresid)
	return mol

def change_mmcif_xyzbfactor(fname,oname,xyzshift,bfactors,badval=1.):
	armed = False

	if fname.endswith('.gz'):
		f = gzip.open(fname, 'rb')
	else:
		f = open(fname,'rb')

	out = ""
	for line in f:
		line = line.decode('utf-8')
		lout = line
		if line.startswith('_atom_site.'):
			armed = True
		else:
			if line.startswith('_') or line.startswith('loop_'):
				armed = False
			elif armed and (not line.startswith('#') and not line.startswith('\n')):
				lout = line
				s = lout.split()
				atomid = int(s[1])
				s[10] = '%07.3f'%(float(s[10])+xyzshift[0])
				s[11] = '%07.3f'%(float(s[11])+xyzshift[1])
				s[12] = '%07.3f'%(float(s[12])+xyzshift[2])
				try:
					s[14] = '%06.2f'%(bfactors[atomid-1])
				except:
					s[14] = '%06.2f'%(badval)
				lout = ' '.join(s) + '\n'

		if not line.startswith('_audit_conform'):
			out += lout
	f.close()
	f = open(oname,'w')
	f.write(out)
	f.close()
