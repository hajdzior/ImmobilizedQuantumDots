# -*- coding: utf-8 -*-

from math import *
from random import *
import pylab as pl
import numpy as np
import sys, getopt
import powerlaw

# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Set up 

# <codecell>


pl.rcParams['xtick.major.pad']='8'
pl.rcParams['ytick.major.pad']='8'
#pylab.rcParams['font.sans-serif']='Arial'


from matplotlib import rc
rc('font', family='sans-serif')
rc('font', size=10.0)
rc('text', usetex=False)

from matplotlib.font_manager import FontProperties

panel_label_font = FontProperties().copy()
panel_label_font.set_weight("bold")
panel_label_font.set_size(12.0)
panel_label_font.set_family("sans-serif")

# <codecell>


def main(argv):
	inputfile = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["ifile="])
	except getopt.GetoptError:
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg


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

#_________ON___________________________
	pl.clf()
	data = on
		
	fit = powerlaw.Fit(data, xmin = min(data))
	label_now = r"Data, $x_{min}=data_{min}=$"+str(fit.power_law.xmin)
	figON = fit.power_law.plot_ccdf(color='b', linestyle='--', label=label_now)
	label_now = r"Fit="+'%.4f' % fit.power_law.alpha+"$\pm$"+'%.4f' % fit.power_law.sigma
	fit.plot_ccdf(color='b', linewidth=2, ax=figON, label=label_now)

	
	fit = powerlaw.Fit(data)
	label_now = r"Data, $x_{min}=$"+str(fit.power_law.xmin)
	fit.power_law.plot_ccdf(color='r', linestyle='--', ax=figON, label=label_now)
	label_now = r"Fit="+'%.4f' % fit.power_law.alpha+"$\pm$"+'%.4f' % fit.power_law.sigma
	fit.plot_ccdf(color='r', linewidth=2, ax=figON, label=label_now)

	fit = powerlaw.Fit(data, xmin = min(data), xmax = fit.power_law.xmin)
	label_now = r"Data, $x_{min}=$"+str(fit.power_law.xmin)
	fit.power_law.plot_ccdf(color='r', linestyle='--', ax=figON, label=label_now)
	label_now = r"Fitkkkkkk="+'%.4f' % fit.power_law.alpha+"$\pm$"+'%.4f' % fit.power_law.sigma
	fit.plot_ccdf(color='g', linewidth=2, ax=figON, label=label_now)
	
	handles, labels = figON.get_legend_handles_labels()
	leg = figON.legend(handles, labels, loc=3)
	leg.draw_frame(True)
	figON.set_ylabel(u"p(X≥x)")
	figON.set_xlabel(r"On times")

	figname = '../tmp/onPlenz'
	pl.savefig(figname+'.png', bbox_inches='tight')


#_________OFF___________________________
	pl.clf()
	data = off
	#fit = powerlaw.Fit(data, discrete=True)
	fit = powerlaw.Fit(data, xmin = min(data))
	label_now = r"Data, $x_{min}=data_{min}=$"+str(fit.power_law.xmin)
	figON = fit.power_law.plot_ccdf(color='b', linestyle='--', label=label_now)
	label_now = r"Fit="+'%.4f' % fit.power_law.alpha+"$\pm$"+'%.4f' % fit.power_law.sigma
	fit.plot_ccdf(color='b', linewidth=2, ax=figON, label=label_now)

	
	fit = powerlaw.Fit(data)
	label_now = r"Data, $x_{min}=$"+str(fit.power_law.xmin)
	fit.power_law.plot_ccdf(color='r', linestyle='--', ax=figON, label=label_now)
	label_now = r"Fit="+'%.4f' % fit.power_law.alpha+"$\pm$"+'%.4f' % fit.power_law.sigma
	fit.plot_ccdf(color='r', linewidth=2, ax=figON, label=label_now)

	
	handles, labels = figON.get_legend_handles_labels()
	leg = figON.legend(handles, labels, loc=3)
	leg.draw_frame(True)
	figON.set_ylabel(u"p(X≥x)")
	figON.set_xlabel(r"Off times")

	figname = '../tmp/offPlenz'
	pl.savefig(figname+'.png', bbox_inches='tight')
									
if __name__ == "__main__":
	main(sys.argv[1:])