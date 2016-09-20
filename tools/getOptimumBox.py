import sys
from math import *



numberOfPointByA = 2
offset = 10 # A

def lineToTab(line):
	line=line.strip()
	line=line.replace("\t", " ")
	while "  " in line:
		line = line.replace("  ", " ")
	return line.split(" ")


fileName = sys.argv[1]


xs=[]
ys=[]
zs=[]
with open(fileName, 'r') as fileIn:
	lines = fileIn.readlines()
	for line in lines[3:]:
		tab = lineToTab(line)
		xs.append(float(tab[4]))
		ys.append(float(tab[5]))
		zs.append(float(tab[6]))

dx = max(xs)-min(xs)
dy = max(ys)-min(ys)
dz = max(zs)-min(zs)

nx=ceil((dx+2*offset)*numberOfPointByA)
ny=ceil((dy+2*offset)*numberOfPointByA)
nz=ceil((dz+2*offset)*numberOfPointByA)

lx=nx/numberOfPointByA
ly=ny/numberOfPointByA
lz=nz/numberOfPointByA


print "boxnod = "+str(nx)+" "+str(ny)+" "+str(nz)
print "boxlen = "+str(lx)+" "+str(ly)+" "+str(lz)
