#! /usr/bin/env python3
# coding=utf-8

# run without arguments for usage info

import re
import sys
import os
from optparse import OptionParser
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

descr = "This Python script allows to visualise a 3D Pareto curve generated by PRISM."
parser = OptionParser(usage="usage: %prog input_file",version="alpha1 (2012-03-04)", description=descr)
(options, args) = parser.parse_args()

if len(args) != 1:
	parser.print_help()
	sys.exit(1)

#open input file
try :
	file = open(args[0])
except IOError as e:
	print "Cannot read the input file '" + e.filename + "': " + e.strerror >> sys.stderr
	sys.exit(1)

col = []
maxX = 0;
maxY = 0;
maxZ = 0;
for line in file:
	tile = []
        coords = line.split(';')
	for pointStr in coords:
		point = pointStr.split(',')
		x= float(point[0])
		y= float(point[1])
		maxX = x if (maxX < x) else maxX  
		maxY = y if (maxY < y) else maxY  
		if len(point) > 2:
			z = float(point[2])
			maxZ = z if (maxZ < z) else maxZ  
			vert = (x,y,z)
		else:
			vert = (x,y)
		tile.append(vert)
	col.append(tile)

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
fig.clf()
dim = len(col[0][0])
if dim == 3:
	ax = fig.add_subplot(111, projection='3d')
	poly = Poly3DCollection(col)
	ax.add_collection3d(poly)
	ax.set_zlabel("objective 3")
	ax.set_zlim(0, maxZ * 1.1)
else:
	ax = fig.add_subplot(111)
	poly = PolyCollection(col)
	ax.add_collection(poly)

ax.set_xlim(0, maxX * 1.1)
ax.set_ylim(0, maxY * 1.1)
ax.set_xlabel("objective 1")
ax.set_ylabel("objective 2")

plt.show()
plt.draw()
