import numpy as np
import matplotlib.pyplot as plt
import sys, getopt 
 
def main(argv):
	inputfile = ''
	outputfile = ''
	try:
			opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
		print 'test.py -i <inputfile> -o <outputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'test.py -i <inputfile> -o <outputfile>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
	#print 'Input file is "', inputfile
	#print 'Output file is "', outputfile
 
 
	# Read data
	f = open (inputfile, "rt")
	l = f.readline()
	x = []
	y = []
	y1 = []
	is_traj = 0
	
	while l:
		va = l.split()
		if (len(va)==5 and va[0]=="record#"):
			is_traj = 1  
			l = f.readline()
			va = l.split()
		if is_traj == 1:
			x.append(float(va[1]))
		l = f.readline()
	f.close() #close the data file

	# Plot data
	fig1 = plt.figure()
	plt.plot(x[0:-2],x[1:-1],'.r')
	plt.xlabel('x')
	plt.ylabel('y')

	# Estimate the 2D histogram
	nbins = 200
	H, xedges, yedges = np.histogram2d(x[0:-2],x[1:-1],bins=nbins)

	# H needs to be rotated and flipped
	H = np.rot90(H)
	H = np.flipud(H)

	# Mask zeros
	Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero

	# Plot 2D histogram using pcolor
	fig2 = plt.figure()
	plt.pcolormesh(xedges,yedges,Hmasked)
	plt.xlabel('$I_{n}$')
	plt.ylabel('$I_{n+1}$')
	cbar = plt.colorbar()
	cbar.ax.set_ylabel('Counts')
	plt.savefig('../tmp/int2dHist.png')
	#plt.show()
	
if __name__ == "__main__":
   main(sys.argv[1:])