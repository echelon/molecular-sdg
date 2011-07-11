#!/usr/bin/env python2.6

import sys

from matrix import MolMatrix
from smiles import Smiles
from smiles import smiles_to_matrix

def main():
	"""Main function"""

	if len(sys.argv) < 2:
		print "Need to supply SMILES text as argument."
		sys.exit()
		
	#mol = MolMatrix(sys.argv[1])
	#smiles = Smiles(sys.argv[1])
	#print smiles.numAtoms()
	mol = smiles_to_matrix(sys.argv[1])


	mol.print_matrix()
	print "\n"
	mol.canonicalize()
	mol.print_matrix()
	print "\n"

	#mol.print_matrix()	

if __name__ == '__main__':
	main()
