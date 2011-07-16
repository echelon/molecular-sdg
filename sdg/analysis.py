# SDG Preassembly Analysis
# From [Helson 1999]
# XXX: This will be *extremely* messy until I have things sorted out
# in my mind as to how the SDG pipeline works. 

from perception.chain import *

def ring_perception(graph):
	# Ring perception greatly aids in the assembly phase. 
	# TODO: Need ring perception algorithm
	# Sources are [Balducci 1994] (recommended), [Figueras 1996]
	# Review [Downs 1989] (Found a copy!)
	pass

def build_datastructs(graph):
	# Helson p 326: Recommends building data structures to supplement
	# mere connection tables... 
	pass

def atom_prioritize(graph):
	# TODO: Won't work until ring perception works.
	# This prioritizes the drawing of congested atoms first
	pass

CHAIN_FIXED_ANGLE = 120

