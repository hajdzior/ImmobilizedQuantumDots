import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import sys, getopt

def main(argv):
	inputfile = ''
	outputfile = ''
	try:
			opts, args = getopt.getopt(argv,"hi:j:o:",["ifile=","ifile2","ofile="])
	except getopt.GetoptError:
		print 'test.py -i <inputfile> -j <inputfile2> -o <outputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'test.py -i <inputfile> -j <inputfile2> -o <outputfile>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-j", "--ifile2"):
			inputfile2 = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
	#print 'Input file is ', inputfile
	#print 'Input file is ', inputfile2
	#print 'Output file is ', outputfile

	# data traj
	is_traj = 0
	f = open (inputfile, "rt")
	l = f.readline()
	x = []
	y = []
	y1 = []
	
	while l:
		va = l.split()
		if (len(va)==5 and va[0]=="record#"):
			is_traj = 1  
			l = f.readline()
			va = l.split()
		if is_traj == 1:
			x.append(float(va[0]))
			y.append(float(va[1]))
		l = f.readline()
	f.close() #close the data file
	
	# data thresh
	f2 = open (inputfile2, "rt")
	l2 = f2.readline()
	th_u = []
	th_d = []

	while l2:
		va = l2.split()
		th_u.append(float(va[0]))
		th_d.append(float(va[1]))
		l2 = f2.readline()
		
	#print th_u
	#print th_d
	nullfmt   = NullFormatter()         # no labels

	# definitions for the axes
	left, width = 0.05, 0.75
	bottom, height = 0.1, 0.8
	bottom_h = left_h = left+width+0.02

	rect_scatter = [left, bottom, width, height]
	rect_histy = [left_h, bottom, 0.15, height]

	## start with a rectangular Figure
	plt.figure(1, figsize=(15,3))

	axScatter = plt.axes(rect_scatter)
	axHisty = plt.axes(rect_histy)

	## no labels
	axHisty.yaxis.set_major_formatter(nullfmt)

	# the scatter plot:
	x_out = x[::10]   # NOTE Bierzemy co 10-ty punkt z trajektorii!
	y_out = y[::10]
	if len(x_out) < 100000:
		axScatter.plot(x_out, y_out,'-')
	else:
		axScatter.plot(x_out[:10000], y_out[:10000],'-')
	axScatter.plot([np.min(x_out),np.max(x_out)],[th_u[0],th_d[0]],'r--')
	axScatter.plot([np.min(x_out),np.max(x_out)],[th_u[1],th_u[1]],'g--')
	axScatter.plot([np.min(x_out),np.max(x_out)],[th_d[1],th_d[1]],'g--')
	axScatter.plot([np.min(x_out),np.max(x_out)],[th_u[2],th_u[2]],'y--')
	axScatter.plot([np.min(x_out),np.max(x_out)],[th_d[2],th_d[2]],'y--')

	# now determine nice limits by hand:
	binwidth = 1
	xmax = np.max(np.fabs(x_out))
	ymax = np.max(np.fabs(y))

	axScatter.set_xlim( (0, xmax) )
	axScatter.set_ylim( (0, ymax) )

	bins = np.arange(0, ymax + binwidth, binwidth)
	hist_val = axHisty.hist(y, bins=bins, orientation='horizontal')
	axHisty.plot([0,np.max(hist_val[0])],[th_u[0],th_d[0]],'r--')
	axHisty.plot([0,np.max(hist_val[0])],[th_u[1],th_u[1]],'g--')
	axHisty.plot([0,np.max(hist_val[0])],[th_d[1],th_d[1]],'g--')
	axHisty.plot([0,np.max(hist_val[0])],[th_u[2],th_u[2]],'y--')
	axHisty.plot([0,np.max(hist_val[0])],[th_d[2],th_d[2]],'y--')
	axHisty.set_ylim( axScatter.get_ylim() )
	
	plt.savefig('../trajsPLOT/TOTAL-' + outputfile)
	#plt.show()
	
	
if __name__ == "__main__":
   main(sys.argv[1:])
