"""
Ring Construction picks up where Ring Analysis left off.
It is the process of laying out the ring coordinates.
"""

from math import *
from cairo import *
import sys

from ring import RING_TYPES

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
def draw_test(ctx, ptA, ptB=None, bondLen=30.0, direc='cw', num=5):
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


# ====================================================
# XXX: Actual work is here. Above is just debug/viz.
# ====================================================

def construct_group(ringGroup):
	"""
	Construct the coordinates, CFS, etc. for a ring group.
	Follows from [Helson] p~335
	"""

	remRings = list(ringGroup[:])

	# Take care of core ring
	core = remRings.pop()
	if core.type != RING_TYPES.CORE:
		print "TODO: RING IS NOT CORE!"
		sys.exit()

	print core[1]
	#sys.exit()
	#core.pos[0].x = foo
	#core.pos[0].y = foo

	# Assign points for core (using CW)

	print "\n\n"
	#points = regular_polygon(core, bondLen=100.0, direc='cw')
	#points = regular_polygon(6, Point(8000.0, -500.0), bondLen=100.0, direc='cw')
	points = regular_polygon(6, Point(0.0, 0.0), 
								Point(150.0, 150.0), 
								bondLen=100.0, direc='cw')

	print points
	print "\n\n"

	return
	# IN LOOP:
	# 1. Get fusion atoms
	# 2. May need to swap fusion atoms to ensure direction still CW
	# 3. Give positions. 

	# Work on remaining rings.
	while len(remRings) > 0:
		ring = remRings.pop()

		if len(done) == 0:
			pass


def attach_fused(ring):
	pass

# TODO for regular_polygon: handle cw/ccw. Still not sure how it works.
def regular_polygon(size, ptA, ptB=None, bondLen=100.0, direc='cw'):
	"""
	Calculate the coordinate positions for a regular polygon.

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
		3. Find vertex points by adding or subtracting phi 
		   once per vertex. (Easy).
	
	Returns: A list of Point object coordinates for each vertex.
	"""
	if direc != 'cw':
		direc = 'ccw'

	# If a second point is not specified, this is probably a core ring,
	# and we are free to align the ring with the coordinate system. 
	# TODO: Update to ensure proper cw/ccw handling
	if ptB == None:
		if size % 2 == 0:
			ptB = Point(ptA.x, ptA.y + bondLen)
		else: 
			ptB = Point(ptA.x + bondLen, ptA.y)

	# Discern bond length.
	# FIXME: Kind of inconsistent to recalculate.
	x = abs(ptA.x - ptB.x)
	y = abs(ptA.y - ptB.y)
	L = sqrt(x**2 + y**2) # Bond length

	# Characteristic angle; the angle between two vertices (from
	# polygon center 'O'), eg. angle <AOB
	phi = radians(360.0)/size

	# Distance from polygon center 'O' to each vertex, eg. O->A
	r = L/(2*sin(phi/2))

	# Point between the line joining A and B
	cX = ptA.x + 0.5*(ptB.x - ptA.x)
	cY = ptA.y + 0.5*(ptB.y - ptA.y)
	ptC = Point(cX, cY)

	# Distance from polygon center 'O' to the point in the 
	# center of line AB, 'C'. ie. O->C
	z = sqrt(r**2 - 0.25*L**2) # Also: z = r*cos(phi/2)

	# Scaled Perpendicular direction vector from C to O
	# http://answers.google.com/answers/threadview/id/419874.html
	k = L/2.0
	u = ((ptA.y - ptC.y)/k,
		(ptC.x - ptA.x)/k)

	# Polygon central point, calculated with the direction vector.
	oX = ptC.x + z* u[0]
	oY = ptC.y + z* u[1]
	ptO = Point(oX, oY)

	# Find angle to point A. 
	adj = ptA.x - ptO.x
	opp = ptA.y - ptO.y

	theta = 0
	if adj == 0.0:
		if opp >= 0:
			theta = pi * 1/2
		else:
			theta = pi * 3/2
	else:
		theta = atan(opp/adj)
	
	# Calculate the positions of each atom
	positions = []
	for x in range(size):
		px = ptO.x + cos(theta) * r
		py = ptO.y + sin(theta) * r
		positions.append(Point(px, py))
		theta += phi # CW

	print ""
	print positions
	return positions

"""
# Entirely TODO
def open_polygon():
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

	pass
"""

