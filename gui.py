import random 
import math
from cairo import *
import pygtk
pygtk.require('2.0')
import gtk

from math import radians, sin, cos
from sdg.ring_construction import *

class Window(object):

	def __init__(self, title='Untitled Window'):
		self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
		self.window.set_title(title)

		self.hbox= gtk.HBox(False, 10)
		self.window.add(self.hbox)

		self.vbox = gtk.VBox(False, 0)
		self.hbox.add(self.vbox)

		self.label = gtk.Label()
		self.vbox.add(self.label)

		self.molEntry = gtk.Entry(max=0)
		self.vbox.add(self.molEntry)

		vbox2 = gtk.VBox(False, 10)
		vbox2.show()
		self.hbox.add(vbox2)

		lab = gtk.Label('Debug Text')
		lab.set_padding(0, 0)
		lab.show()
		vbox2.add(lab)

		self.debugText = gtk.Label()
		self.debugText.set_padding(20, 10)
		vbox2.add(self.debugText)
		
		self.drawable = gtk.DrawingArea()
		self.drawable.set_size_request(500, 500)
		self.vbox.add(self.drawable)

		adj = gtk.Adjustment(0.0, 0.0, 360.0, 1.0, 10.0, 1.0)
		self.scale = gtk.HScale(adj)
		self.vbox.add(self.scale)

		poly = gtk.Adjustment(0, 3, 50, 1, 1, 0)
		spin = gtk.SpinButton(poly, 1.0, 0)
		spin.show()
		self.vbox.add(spin)

		# Callbacks
		self.window.connect('destroy', lambda w: gtk.main_quit())
		self.window.connect('delete_event', lambda x, y: gtk.main_quit())
		self.drawable.connect('expose-event', self.expose)
		self.molEntry.connect('activate', self.textActivate)
		self.scale.connect('value-changed', self.scaleChanged)
		spin.connect('value-changed', self.spinChanged)

		self.atom = None

	def setText(self, text):
		self.label.set_text(text)
		self.molEntry.set_text(text)

	def run(self):
		self.drawable.show()
		self.label.show()
		self.molEntry.show()
		self.debugText.show()
		self.vbox.show()
		self.hbox.show()
		self.scale.show()
		self.window.show()
		gtk.main() 

	def expose(self, widget, event):
		cr = widget.window.cairo_create()
		self.drawLines()

	def scaleChanged(self, scale):
		angle = scale.get_value()
		angle = round(angle)
		self.drawLines(angle=angle)

	def spinChanged(self, spin):
		size = spin.get_value()
		self.drawLines(size=size)

	def textActivate(self, widget):
		text = self.molEntry.get_text()
		self.label.set_text(text)
		self.drawLines()

	def drawLines(self, pt1=None, pt2=None, angle=None, size=None):

		# XXX: Draw whatever is passed from main.py

		ctx = self.drawable.window.cairo_create()

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
		for group in self.ringGroups:
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

			xOff, yOff = yOff, xOff

		# TODO: Cairo Canvas Matrix Transform 
		# Matrix stack: save() restore()
		#ang = radians(180)
		#m = Matrix(cos(ang), sin(ang), -sin(ang), cos(ang), 0, 0) # ROT
		#ctx.transform(m)

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

	# XXX: Deprecated after positions internalized in Ring
	#positions = regular_polygon(size, ptA, ptB, bondLen, direc)
	#draw_spiral2(positions)


