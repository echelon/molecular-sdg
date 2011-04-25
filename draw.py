from atom import Atom, create_graph
from cairo import *
import math


def draw_molecule(ctx, rootMol):
	"""Pass the cairo context and the molecule root"""

	# CLEAR
	pat = SolidPattern(1.0, 1.0, 1.0, 0.9)
	ctx.rectangle(0,0, 500, 500)
	ctx.set_source(pat)
	ctx.fill()

	# POSITION EVERYTHING
	queue = []
	queue.append(rootMol) # root

	while len(queue) > 0:
		mol = queue.pop(0)

		for ch in mol.children:
			queue.append(ch)

		# Root molecule positioning
		if not mol.parent:
			mol.x = 250
			mol.y = 250
			mol.orient = 1
			continue

		# Has a parent...
		parent = mol.parent

		if parent.available['left']:
			parent.positionLeft(mol)

		elif parent.available['up'] and parent.orient == 1:
			parent.positionUp(mol)
		
		elif parent.available['down'] and parent.orient == -1:
			parent.positionDown(mol)

	def draw_line(x1, y1, x2, y2):
		ctx.new_path()
		ctx.set_source_rgb(0.0, 0.0, 0.0)
		ctx.move_to(x1, y1)
		ctx.line_to(x2, y2)
		ctx.close_path()
		ctx.stroke()

	# ACTUAL DRAWING
	queue = []
	queue.append(rootMol)

	while len(queue) > 0:
		mol = queue.pop(0)

		for ch in mol.children:
			queue.append(ch)

		if not mol.parent:
			continue

		parent = mol.parent
		x1 = parent.x
		y1 = parent.y
		x2 = mol.x
		y2 = mol.y

		draw_line(x1, y1, x2, y2)


