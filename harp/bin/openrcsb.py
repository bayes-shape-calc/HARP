#!/usr/bin/env python

import webbrowser
import argparse
import sys

def openrcsb(pdbid):
	if len(pdbid) == 4:
		url = 'https://www.rcsb.org/structure/%s'%(pdbid.lower())
		print('opening %s'%(url))
		webbrowser.open(url)
	else:
		print('malformed PDB ID: %s'%(pdbid))

def main():
	parser = argparse.ArgumentParser(description="Open RCSB PDB page in webbrowser")
	parser.add_argument('pdbid', type=str, help='The ID of the structure in the RCSB PDB to load')
	args = parser.parse_args()

	openrcsb(args.pdbid)

if __name__ == '__main__':
	main()
