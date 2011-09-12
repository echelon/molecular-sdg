#!/usr/bin/env python

"""
This represents the current work on the Molecular Graph and Structure
Diagram Generation (SDG) work that I am doing. There are several files
in this repository that are not relevant anymore, but contain an earlier
approach at drawing that I wish to retain for the time being.

At present, the code I am working on is in 'molecule.py' and 'smiles.py'.
"""

# Python libs
import sys
import random
import cairo
from cairo import *
from glib import markup_escape_text
from math import radians, sin, cos, ceil 

# Parsing, misc.
from gui import Window
from examples import get_example
from smiles import Smiles
from smiles import smiles_to_molecule
from algo.path import *

# Data structures
from molecule import Molecule
from ring import partition_rings
from chain import *

# Perception
from perception.rings import *
from perception.chains import *

# Analysis Phase
from analysis.rings import *

# Assembly Phase
from assembly import * 

# XXX WORK IN PROGRESS.
# Should support everything--chains, rings, etc.
def redraw():
	"""
	Update the image.
	"""
	
	drawable = Globals.drawable
	mol = Globals.molecule
	window = Globals.window

	if not drawable.window or not Globals.molecule:
		# Not ready
		return

	ctx = drawable.window.cairo_create()

	# Size the window
	size = drawable.window.get_size()
	#center = get_average_position(ringGroups)
	
	# CLEAR
	pat = SolidPattern(1.0, 1.0, 1.0, 1.0)
	ctx.rectangle(-8000, -8000, 800000, 800000)
	ctx.set_source(pat)
	ctx.fill()

	# Draw in center of context
	#mat = cairo.Matrix(1, 0, 0, 1, 
	#		size[0]/2 + center[0]*-1, size[1]/2 + center[1]*-1)
	#ctx.transform(mat)

	# XXX XXX DEBUG ONLY
	mat = cairo.Matrix(1, 0, 0, 1,
			size[0]/2 - 100, size[1]/2 - 100)
	ctx.transform(mat)
	# XXX XXX DEBUG ONLY

	# Scale the image
	scale = window.getScale()
	scale *= 0.5 # XXX XXX DEBUG ONLY
	mat = cairo.Matrix(scale, 0, 0, scale, 0, 0)
	ctx.transform(mat)

	window.debugText.hide() # XXX XXX DEBUG ONLY
	window.debugText.set_size_request(0,0) # XXX XXX DEBUG ONLY
	window.scroll.hide()
	window.window.resize(300, 500)

	def draw_edge(ptA, ptB, color, xOff=0, yOff=0, labA=0, labB=0):
		"""
		JUST A TEST. More sophisticated later!
		"""
		ctx.set_source_rgb(color['r'], color['g'], color['b'])
		ctx.new_path()
		ctx.move_to(ptA.x + xOff, ptA.y + yOff)
		ctx.line_to(ptB.x + xOff, ptB.y + yOff)
		ctx.close_path()
		ctx.stroke()

		#print "%d->%d | %s to %s" % (labA, labB, ptA, ptB)

	#print "\n"

	drawnBonds = []

	# XXX DEBUG ONLY: DO NOT DRAW FIRST OR LAST ATOM:
	last = mol.size - 1
	for a in range(mol.size):
		#drawnBonds.append(frozenset([a, 0]))
		#drawnBonds.append(frozenset([a, last]))
		pass

	for atom in range(mol.size):
		# XXX XXX XXX XXX COLOR AND RAND OFFSET HELP DEBUG
		color = {
			'r': random.uniform(0.0, 0.6),
			'g': random.uniform(0.0, 0.6),
			'b': random.uniform(0.0, 0.6),
		}
		randX = 0
		randY = 0
		randX = random.randint(0, 10)
		randY = random.randint(0, 10)

		for nbr in mol.alphaAtoms[atom]:
			bond = frozenset([atom, nbr])
			bond = list(bond)
			atomA = bond[0]
			atomB = bond[1]
			ptA = mol.pos[atom]
			ptB = mol.pos[nbr]

			bond = frozenset(bond)
			if bond in drawnBonds:
				continue
			
			drawnBonds.append(bond)

			draw_edge(ptA, ptB, color, randX, randY, atom, nbr)

	print "..."


def redraw_OLD_BUT_WORKED():
	"""
	Update the image.
	"""
	# XXX: DO NOT CHANGE OTHER GUI COMPONENTS! ONLY IMAGE

	# Extract globals.
	ringGroups = Globals.ringGroups
	drawable = Globals.drawable
	window = Globals.window

	# TODO: Update for chains, etc.
	def get_molecule_dimensions(ringGroups):
		"""From atom data, get the canvas size, (width, height)"""
		maxX = -5000 # These should be out of bounds.
		maxY = -5000
		minX = 100000000
		minY = 100000000
		for group in ringGroups:
			for ring in group:
				for i in range(len(ring)):
					x = ring.pos[i].x
					y = ring.pos[i].y
					if x > maxX:
						maxX = x
					elif x < minX:
						minX = x
					if y > maxY:
						maxY = y
					elif y < minY:
						minY = y
		
		width = maxX - minX
		height = maxY - minY
		return (width, height)

	# TODO: Update for chains, etc.
	def get_average_position(ringGroups):
		xSum = 0
		ySum = 0
		atoms = []
		for group in ringGroups:
			for ring in group:
				for i in range(len(ring)):
					atom = ring[i]
					if atom in atoms:
						continue
					atoms.append(atom)
					x = ring.pos[i].x
					y = ring.pos[i].y
					if x != None:
						xSum += x
					if y != None:
						ySum += y

		xAvg = xSum / len(atoms)
		yAvg = ySum / len(atoms)
		return (int(xAvg), int(yAvg))

	if not drawable.window or not Globals.ringGroups:
		# Not ready
		return

	ctx = drawable.window.cairo_create()

	# Size the window
	size = drawable.window.get_size()
	center = get_average_position(ringGroups)
	
	# CLEAR
	pat = SolidPattern(1.0, 1.0, 1.0, 1.0)
	ctx.rectangle(-8000, -8000, 800000, 800000)
	ctx.set_source(pat)
	ctx.fill()

	# Draw in center of context
	mat = cairo.Matrix(1, 0, 0, 1, 
			size[0]/2 + center[0]*-1, size[1]/2 + center[1]*-1)
	ctx.transform(mat)

	# Scale the image
	scale = window.getScale()
	mat = cairo.Matrix(scale, 0, 0, scale, 0, 0)
	ctx.transform(mat)

	def draw_edge(ptA, ptB, color, xOff=0, yOff=0):
		"""
		JUST A TEST. More sophisticated later!
		"""
		ctx.set_source_rgb(color['r'], color['g'], color['b'])
		ctx.new_path()
		ctx.move_to(ptA.x + xOff, ptA.y + yOff)
		ctx.line_to(ptB.x + xOff, ptB.y + yOff)
		ctx.close_path()
		ctx.stroke()

	# TODO: DEBUG ONLY
	# XXX: Not very well documented..
	labels = {}
	def add_label(pt, atomNum=0):
		x = int(round(ceil(pt.x), -1))
		y = int(round(ceil(pt.y), -1))
		p = "%d,%d" % (x, y)
		if p in labels:
			if atomNum not in labels[p]['atoms']:
				labels[p]['atoms'].append(atomNum)
			return
		labels[p] = {'atoms': [atomNum], 'pos': pt}

	def draw_labels():
		def label_position(pt, label=''):
			ctx.move_to(pt.x, pt.y)
			ctx.set_source_rgb(0, 0, 0)
			ctx.show_text(label)

		for lab in labels.values():
			atoms = "%d" % lab['atoms'][0]
			for atom in lab['atoms'][1:]:
				atoms += ", %d" % atom
			pt = lab['pos']
			txt = "%s (%d, %d)" % (atoms, int(pt.x), int(pt.y))
			label_position(lab['pos'], txt)

	for group in ringGroups:
		for ring in group:
			# XXX XXX XXX XXX COLOR AND RAND OFFSET HELP DEBUG
			color = {
				'r': random.uniform(0.0, 0.6),
				'g': random.uniform(0.0, 0.6),
				'b': random.uniform(0.0, 0.6),
			}
			randX = 0
			randY = 0
			#randX = random.randint(0, 10)
			#randY = random.randint(0, 10)
			for bond in ring.bonds:
				bond = list(bond)
				atomA = bond[0]
				atomB = bond[1]
				ptA = ring.pos[ring.index(atomA)]
				ptB = ring.pos[ring.index(atomB)]

				draw_edge(ptA, ptB, color, randX, randY)

			for i in range(len(ring)):
				atom = ring[i]
				#label_position(ring.pos[i], atom)
				add_label(ring.pos[i], atom)

	draw_labels()

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
	mol.informalName = informalName

	if mol == Globals.molecule:
		print "\n>>> Same molecule.\n"
		return

	debugText = ""

	# Perception algorithms. 
	rings = identify_rings(mol)
	chains = identify_chains(mol, rings)

	mol.setRings(rings)
	mol.setChains(chains)

	# Ring analysis
	ringGroups = partition_rings(rings)
	ring_analysis(ringGroups)
	ring_construction(ringGroups)

	mol.setRingGroups(ringGroups)

	# Assembly phase (TODO: WORK IN PROGRESS)
	assemble(mol)

	# Pango markup for debug window
	debugText += "<b><u><big>Data Structures</big></u></b>\n\n"
	debugText += "<b>Ring Groups</b>: %d\n" % len(ringGroups)
	debugText += "<b>Rings</b>: %d" % len(rings)
	debugText += "\n<b>Chains</b>: %d" % len(chains)
	debugText += "\n\n<big><u><b>Molecule Report</b></u></big>\n\n"
	debugText += "<tt>%s</tt>" % markup_escape_text(str(mol))
	debugText += "\n\n<big><b><u>Ring Group Report</u></b></big>\n"

	for group in ringGroups:
		#construct_group(group)
		# XXX TEMP DEBUG
		debugText += "\n<i>Ring Group (in Peel Order)</i>\n"
		i = 0
		out = group.peelOrder[:]
		out.reverse()
		for ring in out:
			debugText += "\n<b><big>%d</big></b> . %s" % (i, str(ring))
			i += 1

	# Update global information.
	Globals.chains = chains
	Globals.ringGroups = ringGroups
	Globals.molecule = mol

	# Update gui with name, etc.
	window.setSmilesLabel(smiles)
	window.setInformalLabel(informalName)
	window.setDebugLabel(debugText)

class Globals(object):
	"""Used as a Global Dictionary."""
	pass

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

	try:
		# Run GUI.
		win.run()
	except KeyboardInterrupt:
		print ""	

if __name__ == '__main__':
	main()

