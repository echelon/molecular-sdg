"""
Ring Construction picks up where Ring Analysis left off.
It is the process of laying out the ring coordinates.
"""

from math import *
from cairo import *

class Point(object):
	def __init__(self, x, y):
		self.x = x
		self.y = y
	def __str__(self):
		return "(%f, %f)" % (self.x, self.y)
	def __repr__(self):
		return str(self)

# XXX: 'num' is a hack
# TODO: Handle cw/ccw direction.
def draw_test(ctx, ptA, ptB, bondLen=30.0, direc='cw', num=5):
	# CLEAR
	pat = SolidPattern(1.0, 1.0, 1.0, 0.9)
	ctx.rectangle(0,0, 500, 500)
	ctx.set_source(pat)
	ctx.fill()

	def draw_line(ptA, ptB):
		ctx.new_path()
		ctx.set_source_rgb(0.0, 0.0, 0.0)
		ctx.move_to(ptA.x, ptA.y)
		ctx.line_to(ptB.x, ptB.y)
		ctx.close_path()
		ctx.stroke()

	def draw_spiral(center, num, phi, r):
		ctx.set_source_rgb(0.0, 0.0, 0.0)
		# TODO/TEST w/o first or second postions
		theta = 0
		for i in range(num):
			theta += phi
			ctx.new_path()
			px = ptA.x + cos(theta) * r
			py = ptA.y + sin(theta) * r
			ctx.move_to(ptA.x, ptA.y)
			ctx.line_to(px, py)
			ctx.close_path()
			ctx.stroke()

	def draw_spiral2(positions):
		"""
		Draw regular polygons. (WORK IN PROGRESS)
		Input: A list of each vertex position. 
		"""
		# TODO: Must use actual 1st and 2nd positions. 
		# TODO: Use matrix stacks to translate a local coord system.
		positions = positions[:]

		# In order to draw edge from last->first
		positions.append(positions[0])

		first = positions.pop(0)
		next_ = first
		last = 0.0

		ctx.set_source_rgb(0.0, 0.0, 0.0)
		ctx.new_path()
		while len(positions) > 0:
			last = next_
			next_ = positions.pop(0)
			ctx.move_to(last.x, last.y)
			ctx.line_to(next_.x, next_.y)

		ctx.close_path()
		ctx.stroke()

	#draw_line(ptA, ptB)
	#draw_line(pts['o'], pts['c'])

	size = int(num)
	#d = regular_polygon(size, ptA, ptB, bondLen, direc)
	#draw_spiral(d['o'], size, d['phi'], d['r'])

	positions = regular_polygon(size, ptA, ptB, bondLen, direc)
	draw_spiral2(positions)


def regular_polygon(size, ptA, ptB=None, bondLen=100.0, direc='cw'):
	"""
	Calculate the positions for a regular polygon. 

	Inputs:
		size - number of atoms in the ring.
		ptA - the first (x,y) posistion to start with
		ptB - the second (x,y) position. 
			  optional if alignment with x-axis is desired.
		bondLen - desired length of the bonds.
		direc - direction of the construction.

	Objectives: 
		1. Find center point O
		2. Find angle to first vertex
		3. Add or subtract phi once per vertex. (Easy).
	"""
	if direc != 'cw':
		direc = 'ccw'

	# TODO
	# If a second point is not specified, this is probably a core ring,
	# and we are free to align the ring with the coordinate system. 
	if ptB == None:
		if n % 2 == 0:
			# Even-numbered polygon: set ptB's x to same value. 
			ptB = (ptA[0], 0) # TODO: Other coord
		else: 
			# Odd-numbered polygon: set ptB's y to same value.
			ptB = (0, ptA[1]) # TODO: Other coord

	# Discern bond length.
	# TODO: What to do about bond length input?
	x = abs(ptA.x - ptB.x)
	y = abs(ptA.y - ptB.y)
	L = sqrt(x**2 + y**2) # Bond length

	# Characteristic angle; angle between two vertices
	phi = radians(360.0)/size

	# Distance from center to each vertex, eg O->A
	r = L/(2*sin(phi/2))

	# Distance from center to point in center of line AB
	z = sqrt(r**2 - 0.25*L**2) # Also: z = r*cos(phi/2)

	# Point between A and B
	cX = ptA.x + 0.5*(ptB.x - ptA.x)
	cY = ptA.y + 0.5*(ptB.y - ptA.y)
	ptC = Point(cX, cY)

	# Scaled Perpendicular direction vector from C to O
	# http://answers.google.com/answers/threadview/id/419874.html
	k = L/2.0
	u = ((ptA.y - ptC.y)/k,
		(ptC.x - ptA.x)/k)

	# Central point, calculated with direction vector. 
	oX = ptC.x + z* u[0]
	oY = ptC.y + z* u[1]
	ptO = Point(oX, oY)

	# Calculate the positions of each atom
	positions = []
	theta = 0
	for x in range(size):
		px = ptA.x + cos(theta) * r
		py = ptA.y + sin(theta) * r
		positions.append(Point(px, py))
		theta += phi

	return positions
	#return {'o':ptO, 'c': ptC, 'phi': phi, 'r': r}


def open_polygon():
	"""
	Before implementing this, try doing:

		Drawing:
			* Cairo

		Gui:
			* Gtk in /gui
			* Later: Do MVC pattern. Maybe in /mvc
			* Much later: GTK & QT /gui/gtk and /gui/qt, or handle w/ views

		Core datatypes: (LATER)
			* Ring, RingSystem, RingGroup.
			* Chain
			* Molecule
			* Matrix

	"""
	pass

