import pylab as pl
import numpy as np
import sys, getopt
import matplotlib.mlab as mlab
from gaussfitter import multigaussfit
import math as mt
import os
import re

def n_gausses(pars=None,a=None,dx=None,sigma=None):
    """
    Returns a function that sums over N gaussians, where N is the length of
    a,dx,sigma *OR* N = len(pars) / 3

    The background "height" is assumed to be zero (you must "baseline" your
    spectrum before fitting)

    pars  - a list with len(pars) = 3n, assuming a,dx,sigma repeated
    dx    - offset (velocity center) values
    sigma - line widths
    a     - amplitudes
    chii: f.norm -->   The value of the summed squared residuals for the returned parameter
                values.
    """
    if len(pars) % 3 == 0:
        a = [pars[ii] for ii in xrange(0,len(pars),3)]
        dx = [pars[ii] for ii in xrange(1,len(pars),3)]
        sigma = [pars[ii] for ii in xrange(2,len(pars),3)]
    elif not(len(dx) == len(sigma) == len(a)):
        raise ValueError("Wrong array lengths! dx: %i  sigma: %i  a: %i" % (len(dx),len(sigma),len(a)))

    def g(x):
        v = np.zeros((len(dx), len(x)))
        for i in range(len(dx)):
            v[i] = a[i] * np.exp( - ( x - dx[i] )**2 / (2.0*sigma[i]**2) )
        return v
    return g

def plot_fits(hist_x, hist_y, sum_gauss, pars, ngauss_now):
  	# plot histogram and fits
	# histogrm
	pl.plot(hist_x[:-1], hist_y)
	# multigauss fit line
	lab = "n$_{Gauss}$ = " + str(ngauss_now)
	pl.plot(hist_x[:-1],sum_gauss, label = lab)
	pl.legend(ncol=3, loc='best', 
        columnspacing=1.0, labelspacing=0.2,
        handletextpad=0.0, handlelength=1.5,
        fancybox=True, shadow=True, fontsize = 'small')# bbox_to_anchor=[0.5, 1.1]
	# Gauss components
	comp = n_gausses(pars=pars)(hist_x[:-1])
	
	if ngauss_now !=2 and ngauss_now !=6:
		pl.xticks([])
		pl.yticks([])
	for i in range(len(comp)):
	  pl.plot(hist_x[:-1], comp[i],'r--')
    
def calc_AIC(chi, npar, N):
	#chi -> sum of sqr res
	#npar-> nubmer of params in the model
	#N -> sample size
	AIC1 = -N*mt.log10(N) + N*mt.log10(chi)+2*npar
	AIC2 =  N*mt.log10(chi)+2*npar
	return AIC1, AIC2

def scan_dir(dir):
	list_dat = []
	for name in os.listdir(dir):
		path = os.path.join(dir, name)
		if os.path.isdir(path):
			if re.search('.*-trajs$', path): # $:end-line mark
				path_name = re.search('.*-trajs$', path).group(0)
				list_dat.append(path_name)
			else:
				list_dat = list_dat + scan_dir(path)
	return list_dat

def write_results_line(name, index, ngauss_now, pars):
	f3.write(name + ' ')
	for c in range(9):
		if c < ngauss_now:
			f3.write('& %.3f' % pars[c*3 + index])
		else:
			f3.write('&')
	f3.write(' \\\\\n')
	

def fit(ngauss_now):
	print "%d Gauss fit" % ngauss_now

	(pars,sum_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y,ngauss = ngauss_now, params=params, limitedmin=limitedmin_set*ngauss_now, quiet=False)
	chi_norm = chii/len(hist_y)
	
	#-----------------------------------------------------------------
	#!!!!!!!!!!!!!!!!!! I can also add errors!!!!!!!!!!!!!!!!!!!!!!!!!
	#-----------------------------------------------------------------
	
	write_results_line('amplitude', 0, ngauss_now, pars)
	write_results_line('offset', 1, ngauss_now, pars)
	write_results_line('width', 2, ngauss_now, pars)
	f3.write('chi&'+ "%.3f" %(chii) + '& chi norm&'+ "%.3f" %(chi_norm) + '&&&&& \\\\\\hline\n')
	
	
	#AIC model selection
	AIC1.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[0])
	AIC2.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[1])

	pl.subplot(2,4,ngauss_now-1) # WARNING zakladamy ze ngauss = (2,3,4,..,9)
	plot_fits(hist_x, hist_y, sum_gauss, pars, ngauss_now)

	#print "pars", pars, "sum_gauss", sum_gauss, "errors", errors, "chi_norm", chi_norm
	return (pars, sum_gauss, errors, chii)


# --------------- MAIN -------------

dirs_to_analyze = scan_dir('../')
print dirs_to_analyze

for dir_now in dirs_to_analyze:
	print "================================================================================="
	print "I'm working on darta from: " + dir_now
	files_to_analyze = os.listdir(dir_now)
	#print files_to_analyze 
	
	if not os.path.exists(dir_now[:-7] + '/GaussFit'):
		os.mkdir(dir_now[:-7] + '/GaussFit')
		os.mkdir(dir_now[:-7] + '/GaussFit/PCH_gaussFit')
		os.mkdir(dir_now[:-7] + '/GaussFit/PCH_gaussFitPars')
		os.mkdir(dir_now[:-7] + '/GaussFit/PCH_modelSel')
	
	
	for file_now in files_to_analyze:
		print "================================================================================="
		print "I'm working on data from: " + dir_now
		print "---------------------------------------------------------------------------------------"
		print file_now

		inputfile = dir_now + '/' + file_now
		f = open (inputfile, "rt")

		l = f.readline()
		time = []
		counts = []
		is_traj = 1

		while l:
			va = l.split()
			if (len(va)==5 and va[0]=="record#"):
				is_traj = 1  
				l = f.readline()
				va = l.split()
			if is_traj == 1:
				time.append(float(va[0]))
				counts.append(float(va[1]))
			l = f.readline()
		f.close() #close the data file

		# data thresh
		inputfile2 = dir_now[:-7] + 'thresh/' + file_now[:-13] + '-threshVAL.txt'
		f2 = open (inputfile2, "rt")
		l2 = f2.readline()
		th_u = []
		th_d = []
		tmp_f3 = dir_now[:-7] + 'GaussFit/PCH_gaussFitPars/gaussFitTEX_%s.txt' % file_now[0:-4] 
		f3 = open(tmp_f3,'w')
		while l2:
			va = l2.split()
			th_u.append(float(va[0]))
			th_d.append(float(va[1]))
			l2 = f2.readline()

		## calculate histogram
		binwidth = 1
		xmax = np.max(np.fabs(time))
		ymax = np.max(np.fabs(counts))

		bins = np.arange(-0.5, ymax + binwidth-0.5, binwidth)
		hist_y, hist_x = np.histogram(counts, bins=bins)#, normed = True
		hist_x = bins + 0.5

		#find init pars - look for maxima in the left and right parts of the histogram
		first_x = hist_x[0:int(th_u[0])]
		first_y = hist_y[0:int(th_u[0])]
		second_x = hist_x[int(th_u[0]):]
		second_y = hist_y[int(th_u[0]):]

		# Start by fitting 2 gauss functions
		a_f = max(first_y) #amplitude for first max
		o_f = np.where(first_y==max(first_y))[0][0] + 1 # offset for first max
		w_f = mt.sqrt(o_f) #width for first max

		a_s = max(second_y)#amplitude for second max
		o_s = np.where(hist_y==max(second_y))[0][-1]# offset for second max
		w_s = mt.sqrt(o_s)#width  for first max

		#Do we need limitations? [amplitude, offset, width] - default bottom limit is 0
		limitedmin_set=[True, False, True]
		
		AIC1 = []
		AIC2 = []
		
		params = [a_f,o_f,w_f, a_s,o_s,w_s]

		(pars, sum_gauss, errors, chii) = fit(2)
		
		# For gauss number 3 -> 9 add starting offsets for the extra Gauss functions automatically
		for gauss_number in range(3,10):
			differences = [hist_y[i] - sum_gauss[i] for i in range(len(hist_y))]
			
			new_offset = differences.index(max(differences))   #find the place with max difference between fitted and actual values
			new_amp = hist_y[new_offset]
			new_width = mt.sqrt(new_offset)
			
			params.extend((new_amp, new_offset, new_width))
			
			(pars, sum_gauss, errors, chii) = fit(gauss_number)

		outputfile = dir_now[:-7] + '/GaussFit/PCH_gaussFit/PCH_gaussFit' + file_now[0:-4] + '.png'
		pl.savefig(outputfile)
		f3.close()
		pl.clf()

		# Save AIC's to file
		minAIC1 = min(AIC1)
		for l in range(len(AIC1)):
			if (AIC1[l]==minAIC1):
				AIC1[l] = '\\textbf{' + "%.2f" % (AIC1[l])+'}'
			else:
				AIC1[l] = "%.2f" % (AIC1[l])

		AIC1_f = open(dir_now[:-7] + '/GaussFit/PCH_modelSel/aic1Gauss_TEX.txt','a')
		
		stringified_aic1 = ' & '.join([str(i) for i in AIC1])
		AIC1_f.write('$' + inputfile[12:] +'$&'+ stringified_aic1 +'\\\\\n')
		AIC1_f.close()

		minAIC2 = min(AIC2)
		for l in range(len(AIC2)):
			if (AIC2[l]==minAIC2):
				AIC2[l] = '\\textbf{' + "%.2f" % (AIC2[l])+'}'
			else:
				AIC2[l] = "%.2f" % (AIC2[l])
		AIC2_f = open(dir_now[:-7] + '/GaussFit/PCH_modelSel/aic2Gauss_TEX.txt','a')
		AIC2_f.write('$' + inputfile[12:]+'$&'+AIC2[0]+'&'+AIC2[1]+'&'+AIC2[2]+'&'+AIC2[3]+'&'+AIC2[4]+'&'+AIC2[5]+'&'+AIC2[6]+'&'+AIC2[7]+'\\\\\n')#+'&'+AIC2[5]
		AIC2_f.close()
