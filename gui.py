import math
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

		L = 10
		pt1 = Point(250, 260)
		pt2 = Point(120, 180+L)

		if angle != None:
			angle = radians(angle)
			x = pt1.x + cos(angle) * L
			y = pt1.y + sin(angle) * L
			pt2 = Point(x, y)

		if not size:
			size = 5

		# TODO: Matrix Transform 
		# Matrix stack: save() restore()
		#ang = radians(180)
		#m = Matrix(cos(ang), sin(ang), -sin(ang), cos(ang), 0, 0) # ROT
		#ctx.transform(m)

		ctx = self.drawable.window.cairo_create()
		draw_test(ctx, pt1, pt2, num=size)
		#draw_test(ctx, pt1, num=size)


