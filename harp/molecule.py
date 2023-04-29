import numpy as np

class atomcollection(object):
	'''
	Convienience class for keeping track of a group of atoms (i.e. from a molecule/system)
	Also good for extracting subgroups of atoms...

	ind is range(N)
	atomid is a number id to find particular atoms
	resid is the residue ID that the atom is part of
	resiname is the name of the residue (gly, ala etc...)
	chain is the chain ...
	element is (CNOFH etc..)
	atomname is CA CB N, etc...
	xyz is cartesian coordinates... (N,3)
	hetatom is boolean for whether it is a heteroatom entry or not
	authresid is author resid, or if not provided will default to resid
	'''
	def __init__(self,atomid,resid,resname,atomname,chain,element,conf,xyz,occupancy,bfactor,hetatom,modelnum,authresid=None):
		self.atomid = atomid
		self.resid = resid
		self.resname = resname
		self.atomname = atomname
		self.chain = chain
		self.element = element
		self.conf = conf
		self.xyz = np.ascontiguousarray(xyz).astype('double')
		self.occupancy = occupancy
		self.bfactor = bfactor
		self.hetatom = hetatom
		self.modelnum = modelnum
		self.shift = np.zeros(3)
		if authresid is None:
			self.authresid = self.resid.copy()
		else:
			self.authresid = authresid
		self.update()

	def update(self):
		self.unique_residues = np.unique(self.resid)
		self.unique_chains = np.unique(self.chain)
		self.unique_confs = np.unique(self.conf)
		self.natoms = self.atomid.size
		self.ind = np.arange(self.natoms)
		if self.natoms == 1:
			self.xyz = self.xyz.reshape((1,3))

	def get_set(self,keep):
		return atomcollection(self.atomid[keep], self.resid[keep], self.resname[keep], self.atomname[keep],self.chain[keep],self.element[keep], self.conf[keep], self.xyz[keep],self.occupancy[keep],self.bfactor[keep],self.hetatom[keep],self.modelnum[keep],self.authresid[keep])

	def get_residue(self,number):
		return self.get_set(self.resid==number)

	def get_chain(self,chain):
		return self.get_set(self.chain==chain)

	def get_chains(self,chains):
		return self.get_set(np.isin(self.chain,chains))

	def get_atomname(self,atomname):
		return self.get_set(self.atomname==atomname)

	def dehydrogen(self):
		return self.get_set(self.element != "H")

	def get_atomids(self,atomids):
		return self.get_set(np.isin(self.atomid,atomids))

	def get_conformation(self,conf):
		return self.get_set(np.bitwise_or(self.conf=='.',self.conf==conf))

	def remove_hetatoms(self):
		return self.get_set(np.bitwise_not(self.hetatom))

	def bool_NA(self):
		phosphate_atoms = np.array(['P','OP1','OP2','OXT','\"O5\'\"','HOP2','HOP3'])
		an = np.array([ani.upper() for ani in self.atomname])
		potential_phosphate = np.isin(an,phosphate_atoms)
		potential_ribose = np.array([ani.count('\'') > 0 for ani in an]) ## if it's a special one?
		potential_ribose = np.bitwise_and(potential_ribose,an!='\"O5\'\"') ## it gets eaten up to easily account for modifications
		potential_base = np.bitwise_not(np.bitwise_or(potential_phosphate,potential_ribose))
		keep_phosphate = np.zeros_like(potential_phosphate)
		keep_ribose = np.zeros_like(potential_ribose)
		keep_base = np.zeros_like(potential_base)
		for chain in self.unique_chains:
			keepchain = self.chain == chain
			NAs = np.unique(self.resid[np.bitwise_and(keepchain, self.atomname == '\"C4\'\"')])
			keep_na = np.bitwise_and(keepchain, np.isin(self.resid,NAs))
			keep_phosphate += np.bitwise_and(keep_na,potential_phosphate)
			keep_ribose += np.bitwise_and(keep_na,potential_ribose)
			keep_base += np.bitwise_and(keep_na,potential_base)
		return keep_phosphate,keep_ribose,keep_base

	def bool_AA(self):
		backbone_atoms = np.array(['N','H','CA','HCA','HCA1','HCA2','C','O','OXT']) ## backbone and gly
		an = np.array([ani.upper() for ani in self.atomname])
		potential_backbone = np.isin(an,backbone_atoms)
		potential_sidechain = np.bitwise_not(potential_backbone)
		keep_backbone = np.zeros_like(potential_backbone)
		keep_sidechain = np.zeros_like(potential_sidechain)
		for chain in self.unique_chains:
			keepchain = self.chain == chain
			AAs = np.unique(self.resid[np.bitwise_and(keepchain, self.atomname == 'CA')])
			keep_aa = np.bitwise_and(keepchain, np.isin(self.resid,AAs))
			keep_phosphate += np.bitwise_and(keep_aa,potential_backbone)
			keep_sidechain += np.bitwise_and(keep_aa,potential_sidechain)
		return keep_backbone,keep_sidechain

	def split_residue(self):
		if self.unique_residues.size > 1:
			raise Exception('Trying to split more than one residue')

		# if str(self.resname[0]).upper() in ['A','C','G','U','DA','DC','DG','DT']: ## now it doesn't matter b/c there's only one
		if np.any(self.atomname == '\"C4\'\"'): ## the definition of a nucleic acid in Chimera
			phosphate_atoms = np.array(['P','OP1','OP2','OXT','\"O5\'\"','HOP2','HOP3'])
			# ribose_atoms = np.array(['\"C5\'\"','\"C4\'\"','\"C3\'\"','\"O4\'\"','\"C1\'\"','\"C2\'\"','\"O2\'\"','\"O3\'\"','\"H1\'\"','\"H2\'\"','\"H2\'\'\"','\"HO2\'\"','\"H3\'\"','\"HO3\'\"','\"H4\'\"','\"H5\'\"','\"H5\'\'\"'])

			an = np.array([ani.upper() for ani in self.atomname])
			keep_phosphate = np.isin(an,phosphate_atoms)
			# keep_ribose = np.isin(an,ribose_atoms)
			keep_ribose = np.array([ani.count('\'') > 0 for ani in an]) ## if it's a special one?
			keep_ribose = np.bitwise_and(keep_ribose,an!='\"O5\'\"') ## it gets eaten up to easily account for modifications
			keep_base = np.bitwise_not(np.bitwise_or(keep_phosphate,keep_ribose))

			return [self.get_set(keep_phosphate),self.get_set(keep_ribose),self.get_set(keep_base)]

		# elif str(self.resname[0]).upper() in ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']:
		if np.any(self.atomname=='CA'):
			backbone_atoms = np.array(['N','H','CA','HCA','HCA1','HCA2','C','O','OXT']) ## backbone and gly

			an = np.array([ani.upper() for ani in self.atomname])
			keep_backbone = np.isin(an,backbone_atoms)
			keep_sidechain = np.bitwise_not(keep_backbone)

			return [self.get_set(keep_backbone),self.get_set(keep_sidechain)]

		return [self]

	def com(self):
		return self.xyz.mean(0)

	def info_atom(self,i):
		return [self.atomid[i],self.resid[i],self.resname[i],self.chain[i],self.element[i],self.atomname[i],self.conf[i],self.xyz[i],self.occupancy[i],self.bfactor[i],self.hetatom[i],self.modelnum[i],self.authresid[i]]

	def wrapcom_into_unitcell(self,grid):
		shift = np.zeros(3)
		l = grid.nxyz*grid.dxyz
		for i in range(3):
			com = self.com()
			while com[i] < grid.origin[i]:
				self.xyz[:,i] += l[i]
				shift[i] += l[i]
				com = self.com()
		# print(np.all(self.com()>grid.xyzmin())*np.all(self.com()<grid.xyzmax()))
		try:
			self.shift += shift
		except:
			self.shift = shift

def fake_mol(xyz):
	'''
	fake a molecule from xyz coordinate (natoms,3)
	this is enough to work in blobview. Not sure about other situations
	'''
	atomid = np.arange(xyz.shape[0])
	resid = np.zeros(atomid.size,dtype='int')
	resname = np.array(['   ' for _ in range(atomid.size)])
	atomname = np.array(['C' for _ in range(atomid.size)])
	chain = np.array([' ' for _ in range(atomid.size)])
	element = atomname.copy()
	conf = element.copy()
	occupancy = resid.copy()
	bfactor = occupancy.astype('double')
	hetatom = np.array([False for _ in range(atomid.size)])
	modelnum = 	atomname = np.array(['1' for _ in range(atomid.size)])
	return atomcollection(atomid,resid,resname,atomname,chain,element,conf,xyz,occupancy,bfactor,hetatom,modelnum)

def load(fname,only_polymers=False,firstmodel=True,authid=False):
	from . import io
	if fname.endswith('.mmcif') or fname.endswith('.cif') or fname.endswith('.cif.gz'):
		return io.load_mmcif(fname,only_polymers,firstmodel,authid)
