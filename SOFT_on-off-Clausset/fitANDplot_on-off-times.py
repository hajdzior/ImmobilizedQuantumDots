from math import *
from random import *
import pylab as pl
import numpy as np
import sys, getopt

# Import files
import genPL as gPL
import fitPL as fPL
import plotPL as pPL
import plotPL1 as pPL1

def main(argv):
	inputfile = ''
	#outputfile = ''
	try:
		#opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
		opts, args = getopt.getopt(argv,"hi:o:",["ifile="])
	except getopt.GetoptError:
		#print 'test.py -i <inputfile> -o <outputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			#print 'test.py -i <inputfile> -o <outputfile>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		#elif opt in ("-o", "--ofile"):
			#outputfile = arg

	# data
	f = open (inputfile, "rt")
	#print inputfile
	l = f.readline()
	on = []
	off = []
	while l:
		va = l.split()
		on.append(float(va[0]))
		off.append(float(va[1]))
		l = f.readline()

	## simulate
	#x = gPL.randht(nsim,'powerlaw',alpha);

	# fit
	[alphaON, xminON, LON] = fPL.plfit(on)
	[alphaOFF, xminOFF, LOFF] = fPL.plfit(off)

	# print
	pl.figure(1)
	pl.subplot(1,2,1)
	pl.plot(on,'g*')
	pl.subplot(1,2,2)
	pl.hist(on, 100)
	pl.savefig('../tmp/Raw-on.png')

	pl.figure(2)
	pl.subplot(1,2,1)
	h1 = pPL1.plplot(on,xminON,alphaON,'on time')
	textstr = '$a =%.2f$\n $xmin=%.2f$\n $L=%.2f$'%(alphaON, xminON, LON)
	pl.text(0.05, 0.95, textstr, fontsize=14, verticalalignment='top')
	pl.subplot(1,2,2)
	h = pPL.plplot(on,xminON,alphaON,'on time','ro', 'k--')
	pl.savefig('../tmp/Fits-on.png')
	
	pl.figure(3)
	pl.subplot(1,2,1)
	pl.plot(off,'g*')
	pl.subplot(1,2,2)
	pl.hist(off, 100)
	pl.savefig('../tmp/Raw-off.png')

	pl.figure(4)
	pl.subplot(1,2,1)
	h1 = pPL1.plplot(off,xminOFF,alphaOFF,'off time')
	textstr = '$a =%.2f$\n $xmin=%.2f$\n $L=%.2f$'%(alphaOFF, xminOFF, LOFF)
	pl.text(0.05, 0.95, textstr, fontsize=14, verticalalignment='top')
	pl.subplot(1,2,2)
	h = pPL.plplot(off,xminOFF,alphaOFF,'off time','g*', 'k--')
	pl.savefig('../tmp/Fits-off.png')
	#pl.show()

	pl.figure(5)
	h = pPL.plplot(on,xminON,alphaON,'times','r-', 'k--')
	h = pPL.plplot(off,xminOFF,alphaOFF,'times','g-', 'k--')

	pl.savefig('../tmp/comp.png')
if __name__ == "__main__":
	main(sys.argv[1:])
