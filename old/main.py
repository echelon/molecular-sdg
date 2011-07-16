#!/usr/bin/env python
from atom import Atom, create_graph
from cairo import *
import math
from gui import Window

# Examples 
hexane = "CCCCCC"
hexene = "C=CCCCC"						# 1-hexene
hexyne = "C#CCCCC"						# 1-hexyne
isohex = "CC(C)CCC"						# isohexane or 2-methylpentane
adenine = "n1c(c2c(nc1)ncn2)N"
acetic = "CC(=O)O" 						# Acetic Acid
pbromo = "C1=CC(=CC=C1Cl)Br" 			# p-BromoChloro Benzene 
vanillin = "O=Cc1ccc(O)c(OC)c1"
naph = "c1cccc2c1cccc2" 				# Naphthalene

"""mol = create_graph(hexane)
create_graph(hexene)
create_graph(hexyne)
create_graph(isohex)

surface = ImageSurface(FORMAT_ARGB32, 500, 500)
ctx = Context(surface)

# Background
height = 500
pat = SolidPattern(1.0, 1.0, 1.0, 0.9)
ctx.rectangle (0,0, 500, 500)
ctx.set_source (pat)
ctx.fill ()
"""


win = Window("Testing")

win.setText(hexane)
win.run()




def draw_line(x1, y1, x2, y2):
	ctx.new_path()
	ctx.set_source_rgb(0.0, 0.0, 0.0)
	ctx.move_to(x1, y1)
	ctx.line_to(x2, y2)
	ctx.close_path()
	ctx.stroke()

def write_text(x, y, text):
	ctx.set_source_rgb(0.0, 0.0, 0.0)
	ctx.select_font_face ("Sans", FONT_SLANT_NORMAL, FONT_WEIGHT_BOLD)
	ctx.set_font_size(21)
	x_bearing, y_bearing, width, height, \
			x_advance, y_advance = ctx.text_extents(text)

	ctx.move_to(x, y)
	ctx.show_text(text)
	ctx.fill()

#write_text(400, 200, 'N')

def test(x, y, direc='1', length=30):
	ANG = 180 - (109/2 + 90)
	direc = 1 if direc > 0 else -1
	endX = x + length
	endY = y + math.tan(ANG) * direc
	print ANG
	print math.tan(ANG)
	draw_line(x, y, endX, endY)
	return (endX, endY)

#xy = test(100, 100)
#xy = test(xy[0], xy[1], -1)
#xy = test(xy[0], xy[1], 1)
#xy = test(xy[0], xy[1], -1)
#xy = test(xy[0], xy[1], 1)

def position_molecule(root, x=100, y=100):

	def position_atom(root, x, y, direc):
		root.x = x
		root.y = y
		for child in root.children:
			position_atom(root)

	position_atom(root, x, y, 'l')

#draw_molecule(mol)
#surface.write_to_png("output.png")

