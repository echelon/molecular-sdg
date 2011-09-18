import random 
import math
import pygtk
pygtk.require('2.0')
import gtk

from math import radians, sin, cos

class Window(object):

	def __init__(self, title='Untitled Window'):
		self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
		self.window.set_title(title)

		self.hbox= gtk.HBox(False, 0)
		self.hbox.set_spacing(0)
		self.window.add(self.hbox)

		self.vbox = gtk.VBox(False, 0)
		self.hbox.add(self.vbox)


		self.label = gtk.Label()
		self.vbox.add(self.label)

		self.molEntry = gtk.Entry(max=0)
		self.vbox.add(self.molEntry)
		
		scroll = gtk.ScrolledWindow()
		scroll.show()
		scroll.set_border_width(10)
		self.scroll = scroll
		self.hbox.add(scroll)

		self.informalNameLabel = gtk.Label() # TODO: Remove.

		self.debugText = gtk.Label()
		self.debugText.set_padding(0, 0)
		self.debugText.set_line_wrap(True)
		scroll.add_with_viewport(self.debugText)
		scroll.set_size_request(450, 100) # Min size
		#print scroll.get_placement()
		#scroll.set_placement(gtk.CORNER_BOTTOM_RIGHT)

		self.drawable = gtk.DrawingArea()
		self.drawable.set_size_request(500, 500)
		self.vbox.add(self.drawable)

		adj = gtk.Adjustment(0.5, 0.2, 1.5)
		self.scale = gtk.HScale(adj)
		self.vbox.add(self.scale)

		#poly = gtk.Adjustment(0, 3, 50, 1, 1, 0)
		#spin = gtk.SpinButton(poly, 1.0, 0)
		#spin.show()
		#self.vbox.add(spin)

		def exitkey(ev):
			# ESC 
			if ev.keyval == gtk.keysyms.Escape:
				gtk.main_quit()

			# CTRL-C / CTRL-D
			keyname = gtk.gdk.keyval_name(ev.keyval)
			ekey = ['C', 'D']
			if ev.state == gtk.gdk.CONTROL_MASK and keyname.upper() in ekey:
				gtk.main_quit()

			return False # Bubble event

		# Callbacks
		self.window.connect('destroy', lambda w: gtk.main_quit())
		self.window.connect('delete_event', lambda w, e: gtk.main_quit())
		self.window.connect('key-press-event', lambda w, e: exitkey(e))
		self.drawable.connect('expose-event', self.expose)
		self.molEntry.connect('activate', self.smilesChanged)
		self.scale.connect('value-changed', self.scaleChanged)
		#spin.connect('value-changed', self.spinChanged)

		# Must set these callbacks
		self.drawCallback = None
		self.smilesCallback = None

	def setSmilesLabel(self, smiles, todoRemove=None): # TODO
		if not self.molEntry.get_text():
			self.molEntry.set_text(smiles)
		self.label.set_text(smiles)

	def setInformalLabel(self, informalName):
		self.informalNameLabel.set_markup("<b>Informal name:</b>\n%s" 
				% str(informalName))

	def setDebugLabel(self, debugText):
		self.debugText.set_markup(str(debugText))

	def getScale(self):
		return self.scale.get_value()

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
		self.draw()

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

