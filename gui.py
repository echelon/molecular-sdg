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

		self.informalNameLabel = gtk.Label()
		self.setInformalLabel(None)
		vbox2.add(self.informalNameLabel)

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
		self.molEntry.connect('activate', self.smilesChanged)
		self.scale.connect('value-changed', self.scaleChanged)
		spin.connect('value-changed', self.spinChanged)

		# Must set these callbacks
		self.drawCallback = None
		self.smilesCallback = None

	def setSmilesLabel(self, smiles, todoRemove=None): # TODO
		self.label.set_text(smiles)
		self.molEntry.set_text(smiles)

	def setInformalLabel(self, informalName):
		self.informalNameLabel.set_text("Informal name:\n%s" 
				% str(informalName))

	def run(self):
		self.drawable.show()
		self.label.show()
		self.informalNameLabel.show()
		self.molEntry.show()
		self.debugText.show()
		self.vbox.show()
		self.hbox.show()
		self.scale.show()
		self.window.show()
		gtk.main() 

	def expose(self, widget, event):
		self.draw()

	def scaleChanged(self, scale):
		angle = scale.get_value()
		angle = round(angle)
		pass # TODO

	def spinChanged(self, spin):
		size = spin.get_value()
		pass # TODO

	def smilesChanged(self, widget):
		text = self.molEntry.get_text()
		
		# Callback.
		if not self.smilesCallback:
			print "No smiles callback set."
			return
		self.smilesCallback(text)

		self.draw()

	def draw(self):
		if not self.drawCallback:
			print "No draw callback set."
			return
		self.drawCallback()

