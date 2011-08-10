#!/usr/bin/env python

"""
This represents the current work on the Molecular Graph and Structure
Diagram Generation (SDG) work that I am doing. There are several files
in this repository that are not relevant anymore, but contain an earlier
approach at drawing that I wish to retain for the time being.

At present, the code I am working on is in 'matrix.py' and 'smiles.py'.
"""

# Python libs
import sys
import random 
from cairo import *
from math import radians, sin, cos

# Parsing, misc.
from gui import Window
from examples import get_example
from smiles import Smiles
from smiles import smiles_to_molecule
from algo.path import *

# Data structures
from molecule import Molecule
from ring import partition_rings

# Perception
from perception.rings import *
from perception.chain import *

# Analysis Phase
from sdg.ring_analysis import *
from sdg.ring_construction import *

class Globals(object):
	"""Used as a Global Dictionary."""
	pass

def redraw():
	"""
	Update the image.
	"""
	# XXX: DO NOT CHANGE OTHER GUI COMPONENTS! ONLY IMAGE

	# Extract globals.
	ringGroups = Globals.ringGroups
	drawable = Globals.drawable

	if not drawable.window or not Globals.ringGroups:
		# Not ready
		return

	ctx = drawable.window.cairo_create()

	# CLEAR
	pat = SolidPattern(1.0, 1.0, 1.0, 0.9)
	ctx.rectangle(0,0, 500, 500)
	ctx.set_source(pat)
	ctx.fill()

	def draw_edge(ptA, ptB, color, xOff=100, yOff=100):
		"""
		JUST A TEST. More sophisticated later!
		"""
		ctx.set_source_rgb(color['r'], color['g'], color['b'])
		ctx.new_path()
		ctx.move_to(ptA.x + xOff, ptA.y + yOff)
		ctx.line_to(ptB.x + xOff, ptB.y + yOff)
		ctx.close_path()
		ctx.stroke()

	xOff = 50
	yOff = 50
	for group in ringGroups:
		xOff += 100
		yOff += 100
		for ring in group:
			print ring
			# XXX XXX XXX XXX COLOR AND RAND OFFSET HELP DEBUG
			color = {
				'r': random.uniform(0.0, 0.6),
				'g': random.uniform(0.0, 0.6),
				'b': random.uniform(0.0, 0.6),
			}
			randX = random.randint(0, 10)
			randY = random.randint(0, 10)
			for bond in ring.bonds:
				bond = list(bond)
				atomA = bond[0]
				atomB = bond[1]
				ptA = ring.pos[ring.index(atomA)]
				ptB = ring.pos[ring.index(atomB)]

				draw_edge(ptA, ptB, color, xOff + randX, yOff + randY)

		# Switch out the random spread for ring groups
		xOff, yOff = yOff, xOff

# TODO: Needs rename
def parse_smiles_text(smiles, informalName=None):
	"""
	Process SMILES text into molecular information and a structure
	diagram.
	"""
	# Extract globals
	window = Globals.window
	debugText = Globals.debugText

	mol = smiles_to_molecule(smiles)

	if mol == Globals.molecule:
		print "\n>>> Same molecule.\n"
		return

	#mol.print_matrix()

	# Perception algorithms. 
	rings = identify_rings(mol)
	chains = identify_chains(mol, rings)

	# Ring analysis
	ringGroups = partition_rings(rings)

	#print ringGroups

	ring_analysis(ringGroups, mol)

	#for rg in ringGroups:
	#	print rg.peelOrder

	for group in ringGroups:
		#print group
		construct_group(group)

	# Update global information.
	Globals.chains = chains
	Globals.ringGroups = ringGroups
	Globals.molecule = mol

	# Update gui with name, etc.
	window.setSmilesLabel(smiles)
	window.setInformalLabel(informalName)

def gui_draw_callback():
	redraw()

def gui_smiles_text_callback(smiles):
	"""
	Take form input and attempt a SMILES dictionary lookup.
	Otherwise considered actual SMILES input.
	"""
	informalName = None
	ex = get_example(smiles)
	if ex:
		informalName = ex[1]
		smiles = ex[2]
		print ">>> Using %s per input.\n" % informalName

	parse_smiles_text(smiles, informalName)

def main():
	"""Main function"""
	# Extract arguments or get random SMILES
	smiles = None
	informalName = None
	ex = get_example(None if len(sys.argv) < 2 else sys.argv[1])
	if not ex:
		smiles = sys.argv[1]
		print "\n>>> Using input `%s` as SMILES.\n" % smiles
	else:
		informalName = ex[1]
		smiles = ex[2]
		print "\n>>> Using %s per argument.\n" % informalName

	# Init GUI.
	win = Window('Chemical Structure Diagram Generation (WIP)')

	# Set globals.
	Globals.drawable = win.drawable
	Globals.debugText = win.debugText
	Globals.window = win
	Globals.ringGroups = None
	Globals.molecule = None

	# Set callbacks.
	win.drawCallback = gui_draw_callback
	win.smilesCallback = gui_smiles_text_callback

	# Process initial data.
	parse_smiles_text(smiles, informalName)

	# Run GUI.
	win.run()

if __name__ == '__main__':
	main()

