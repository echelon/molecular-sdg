import math
import pygtk
pygtk.require('2.0')
import gtk

from atom import Atom, create_graph
from draw import draw_molecule
from smiles import Smiles

class Window(object):

	def __init__(self, title):
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

		# Callbacks
		self.window.connect('destroy', lambda w: gtk.main_quit())
		self.window.connect('delete_event', lambda x, y: gtk.main_quit())
		self.drawable.connect('expose-event', self.expose)
		self.molEntry.connect('activate', self.textActivate)

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
		self.window.show()
		gtk.main() 

	def expose(self, widget, event):

		cr = widget.window.cairo_create()
		
		#cr.set_line_width(9)
		#cr.set_source_rgb(0.7, 0.2, 0.0)
		
		#w = self.window.allocation.width
		#h = self.window.allocation.height
		
		#cr.translate(w/2, h/2)
		#cr.arc(0, 0, 50, 0, 2*math.pi)
		#cr.stroke_preserve()
		#cr.set_source_rgb(0.3, 0.4, 0.6)
		#cr.fill()


	def textActivate(self, widget):
		text = self.molEntry.get_text()
		self.label.set_text(text)

		dbg = ""

		smiles = Smiles(text)
		numAtoms = smiles.num_atoms()
		dbg += "\n\nNum Atoms: %d" % numAtoms

		self.atom = create_graph(text)	
		dbg += "\n\nGraph\n===============\n"
		dbg += str(self.atom)

		self.debugText.set_text(dbg)

		ctx = self.drawable.window.cairo_create()
		draw_molecule(ctx, self.atom)


