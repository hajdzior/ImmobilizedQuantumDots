import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import sys, getopt

PLOT_SIZE = 5000
PLOTS_PER_PAGE = 9

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

	total_plots = len(x) / PLOT_SIZE
	
	for page in xrange(total_plots / PLOTS_PER_PAGE + 1):
	    plt.figure(page, figsize=(19,28))

	    if total_plots >= PLOTS_PER_PAGE:
	      num_plots = PLOTS_PER_PAGE
	    else:
	      num_plots = total_plots

	    for p in xrange(num_plots):
		drawSubPlot(page, p, x, y, th_u, th_d)
		#drawHistogram(page, p, x, y, th_u, th_d)
	
	    total_plots -= PLOTS_PER_PAGE
	
	    plt.savefig('../trajsPLOT/%d-' % page + outputfile)
	    #plt.show()
	
def createNewSubPlot(p, histogram = False):
  	# definitions for the axes
	if histogram:
	    left, width = 0.05, 0.75 # with histogram
	else:
	    left, width = 0.05, 0.9  # no hist
	bottom, height = 0.1 * (PLOTS_PER_PAGE - p - 1) + 0.03, 0.083
	left_h = left + width + 0.02

	rect_scatter = [left, bottom, width, height]
	rect_histy = [left_h, bottom, 0.15, height]

	axScatter = plt.axes(rect_scatter)
	if histogram:
	  axHisty = plt.axes(rect_histy)
	else:
	  axHisty = None
	
	return axScatter, axHisty

def drawSubPlot(page, p, x, y, th_u, th_d):
	axScatter, axHisty = createNewSubPlot(p)

	# the scatter plot:
	plot_num = p + page * PLOTS_PER_PAGE
	start = plot_num * PLOT_SIZE
	end = start + PLOT_SIZE
	axScatter.plot(x[start:end], y[start:end],'-')
	
	axScatter.plot([x[start],x[end]],[th_u[0],th_d[0]],'r--')
	axScatter.plot([x[start],x[end]],[th_u[1],th_u[1]],'g--')
	axScatter.plot([x[start],x[end]],[th_d[1],th_d[1]],'g--')
	axScatter.plot([x[start],x[end]],[th_u[2],th_u[2]],'y--')
	axScatter.plot([x[start],x[end]],[th_d[2],th_d[2]],'y--')
	axScatter.set_xlim([x[start],x[end]])
	axScatter.set_ylim([0, np.max(np.fabs(y))])

def drawHistogram(page, p, x, y, th_u, th_d):
	axScatter, axHisty = createNewSubPlot( p)
	nullfmt   = NullFormatter()         # no labels
	axHisty.yaxis.set_major_formatter(nullfmt)
	
	# now determine nice limits by hand:
	binwidth = 1
	#xmax = np.max(np.fabs(x))
	ymax = np.max(np.fabs(y))

	bins = np.arange(0, ymax + binwidth, binwidth)
	hist_val = axHisty.hist(y, bins=bins, orientation='horizontal')
	axHisty.plot([0,np.max(hist_val[0])],[th_u[0],th_d[0]],'r--')
	axHisty.plot([0,np.max(hist_val[0])],[th_u[1],th_u[1]],'g--')
	axHisty.plot([0,np.max(hist_val[0])],[th_d[1],th_d[1]],'g--')
	axHisty.plot([0,np.max(hist_val[0])],[th_u[2],th_u[2]],'y--')
	axHisty.plot([0,np.max(hist_val[0])],[th_d[2],th_d[2]],'y--')
	axHisty.set_ylim( axScatter.get_ylim() )
	
if __name__ == "__main__":
   main(sys.argv[1:])
