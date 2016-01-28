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

##_________ON___________________________
	pl.clf()
	data = on
	#fit = powerlaw.Fit(data, discrete=True)
	fit = powerlaw.Fit(data, xmin = min(data))
	label_now = r"Data: $x_{min}=data_{min}=$"+str(fit.truncated_power_law.xmin)
	figON = fit.plot_ccdf(linewidth=2, label=label_now)
	label_now = r"Fit, $\alpha$ = "+'%.4f' % fit.truncated_power_law.alpha+"; $\lambda = $"+'%.4f' % fit.truncated_power_law.Lambda
	fit.truncated_power_law.plot_ccdf(color='b', linestyle='--', ax=figON, label=label_now)
	

	
	fit = powerlaw.Fit(data)
	label_now = r"Data: $x_{min}=$"+str(fit.truncated_power_law.xmin)
	fit.plot_ccdf(color='r', linewidth=2, ax=figON, label=label_now)
	label_now = r"Fits, $\alpha$ = " + '%.4f' % fit.truncated_power_law.alpha+"; $\lambda = $"+'%.4f' % fit.truncated_power_law.Lambda
	fit.truncated_power_law.plot_ccdf(color='r', linestyle='--', ax=figON, label=label_now)
	
	handles, labels = figON.get_legend_handles_labels()
	leg = figON.legend(handles, labels, loc=3)
	leg.draw_frame(True)
	figON.set_ylabel(u"p(X≥x)")
	figON.set_xlabel(r"On times")

	figname = '../tmp/onPlenz'
	pl.savefig(figname+'.png', bbox_inches='tight')

	##_________OFF___________________________
	pl.clf()
	data = off
	#fit = powerlaw.Fit(data, discrete=True)
	fit = powerlaw.Fit(data, xmin = min(data))
	label_now = r"Data: $x_{min}=data_{min}=$"+str(fit.truncated_power_law.xmin)
	figOFF = fit.plot_ccdf(linewidth=2, label=label_now)
	label_now = r"Fit, $\alpha$ = "+'%.4f' % fit.truncated_power_law.alpha+"; $\lambda = $"+'%.4f' % fit.truncated_power_law.Lambda
	fit.truncated_power_law.plot_ccdf(color='b', linestyle='--', ax=figOFF, label=label_now)
		
	fit = powerlaw.Fit(data)
	label_now = r"Data: $x_{min}=$"+str(fit.truncated_power_law.xmin)
	fit.plot_ccdf(color='r', linewidth=2, ax=figOFF, label=label_now)
	label_now = r"Fits, $\alpha$ = " + '%.4f' % fit.truncated_power_law.alpha+"; $\lambda = $"+'%.4f' % fit.truncated_power_law.Lambda
	fit.truncated_power_law.plot_ccdf(color='r', linestyle='--', ax=figOFF, label=label_now)
	
	handles, labels = figOFF.get_legend_handles_labels()
	leg = figOFF.legend(handles, labels, loc=3)
	leg.draw_frame(True)
	figOFF.set_ylabel(u"p(X≥x)")
	figOFF.set_xlabel(r"Off times")

	figname = '../tmp/offPlenz'
	pl.savefig(figname+'.png', bbox_inches='tight')

	
	
	## COMPARE ----------------------------------------------------------------
	#___________ON______________________________
	pl.clf()
	data = on
	fit = powerlaw.Fit(data)
	####
	fit.distribution_compare('power_law', 'truncated_power_law')
	fig = fit.plot_ccdf(linewidth=3, label='Empirical Data')
	fit.power_law.plot_ccdf(ax=fig, color='r', linestyle='--', label='Power law fit')
	fit.truncated_power_law.plot_ccdf(ax=fig, color='g', linestyle='--', label='Truncated power law fit')
	####
	fig.set_ylabel(u"p(X≥x)")
	fig.set_xlabel("On times")
	handles, labels = fig.get_legend_handles_labels()
	fig.legend(handles, labels, loc=3)

	
	figname = '../tmp/Fig_on_PL_truncPL'
	
	R_n,p_n = fit.distribution_compare('power_law', 'truncated_power_law', normalized_ratio=True)
	R,p =  fit.distribution_compare('power_law', 'truncated_power_law')
	
	textstr = '$R$ =' + '%.4f' % R + '\n $p$ = ' + '%.4f' % p +'\n $R_{norm}$ =' + '%.4f' % R_n + "\n $p_{norm}$ = " + '%.4f' % p_n
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	fig.text(0.05, 0.45, textstr, transform=fig.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
	pl.savefig(figname+'.png', bbox_inches='tight')
	
	# ______OFF______________
	pl.clf()
	data = off
	fit = powerlaw.Fit(data)
	####
	fit.distribution_compare('power_law', 'truncated_power_law')
	fig = fit.plot_ccdf(linewidth=3, label='Empirical Data')
	fit.power_law.plot_ccdf(ax=fig, color='r', linestyle='--', label='Power law fit')
	fit.truncated_power_law.plot_ccdf(ax=fig, color='g', linestyle='--', label='Truncated power law fit')
	####
	fig.set_ylabel(u"p(X≥x)")
	fig.set_xlabel("Off times")
	handles, labels = fig.get_legend_handles_labels()
	fig.legend(handles, labels, loc=3)

	
	figname = '../tmp/Fig_off_PL_truncPL'
	
	R_n,p_n = fit.distribution_compare('power_law', 'truncated_power_law', normalized_ratio=True)
	R,p =  fit.distribution_compare('power_law', 'truncated_power_law')
	
	textstr = '$R$ =' + '%.4f' % R + '\n $p$ = ' + '%.4f' % p +'\n $R_{norm}$ =' + '%.4f' % R_n + "\n $p_{norm}$ = " + '%.4f' % p_n
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	fig.text(0.05, 0.45, textstr, transform=fig.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
	pl.savefig(figname+'.png', bbox_inches='tight')
	
	#pl.show()
	#------------------------------------------------------------------
if __name__ == "__main__":
	main(sys.argv[1:])