from ..molecule import atomcollection
import numpy as np

def load_pdb(fname):
	# http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
	## document is 1 indexed. ...
	with open(fname,'r') as f:

		atomid = []
		atomname = []
		resid = []
		resname = []
		chain = []
		occupancy = []
		bfactor = []
		element = []
		hetatom = []
		conf = []
		xyz = []
		for line in f:

			recordname = line[:6]
			if recordname == "ATOM  ":
				hetatom.append(False)
			elif recordname == "HETATM":
				hetatom.append(True)
			else:
				continue

			atomid.append(int(line[7-1:11])) ## serial
			atomname.append(line[13-1:16].replace(' ',''))
			# altloc.append(line[17-1])
			resname.append(line[18-1:20].replace(' ',''))
			chain.append(line[22-1])
			resid.append(int(line[23-1:26]))
			# icode.append(line[27-1])
			x = float(line[31-1:38])
			y = float(line[39-1:46])
			z = float(line[47-1:54])
			xyz.append([x,y,z])
			occupancy.append(float(line[55-1:60]))
			bfactor.append(float(line[61-1:66]))
			element.append(line[77-1:78])
			# charge.append(line[79-1:80])
			conf.append('.')

	return atomcollection(np.array(atomid),np.array(resid),np.array(resname),np.array(atomname),np.array(chain),np.array(element),np.array(conf),np.array(xyz),np.array(occupancy),np.array(bfactor),np.array(hetatom))


def write_pdb(fname,mol):
	with open(fname,'w') as f:
		for i in range(mol.natoms):
			line = ' '*80+'\n'
			start = 'HETATM' if mol.hetatom[i] else "ATOM  "

			line = "%6s%5d %4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"%(start,
				mol.atomid[i],
				mol.atomname[i],
				mol.resname[i],
				mol.chain[i],
				mol.resid[i],
				mol.xyz[i,0],
				mol.xyz[i,1],
				mol.xyz[i,2],
				mol.occupancy[i],
				mol.bfactor[i],
				mol.element[i]
			)

			f.write(line)
