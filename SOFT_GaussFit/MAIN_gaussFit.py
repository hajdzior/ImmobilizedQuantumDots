import pylab as pl
import numpy as np
import sys, getopt
import matplotlib.mlab as mlab
from gaussfitter import multigaussfit
from gaussfitter import n_gaussian
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

def plot_fits(hist_x, hist_y, n_gauss, pars):
  	# plot histogram and fits
	# histogrm
	pl.plot(hist_x[:-1], hist_y)
	# multigauss fit line
	lab = "n$_{Gauss}$ = " + str(ngauss_now)
	pl.plot(hist_x[:-1],n_gauss, label = lab)
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

#def main(argv):
    #inputfile = ''
    #outputfile = ''
    #try:
		    #opts, args = getopt.getopt(argv,"hi:j:o:",["ifile=","ifile2","ofile="])
    #except getopt.GetoptError:
	    #print 'test.py -i <inputfile> -j <inputfile2> -o <outputfile>'
	    #sys.exit(2)
    #for opt, arg in opts:
	    #if opt == '-h':
		    #print 'test.py -i <inputfile> -j <inputfile2> -o <outputfile>'
		    #sys.exit()
	    #elif opt in ("-i", "--ifile"):
		    #inputfile = arg
	    #elif opt in ("-j", "--ifile2"):
		    #inputfile2 = arg
	    #elif opt in ("-o", "--ofile"):
		    #outputfile = arg

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

dirs_to_analyze = scan_dir('../')
print dirs_to_analyze

for dir_now in dirs_to_analyze:
	print "================================================================================="
	print "I'm working on darta from: " + dir_now
	files_to_analyze = os.listdir(dir_now)
	#print files_to_analyze 
	
	if not os.path.exists(dir_now[:-7] + '/tmp'):
		os.mkdir(dir_now[:-7] + '/tmp')
	if not os.path.exists(dir_now[:-7] + '/GaussFit'):
		os.mkdir(dir_now[:-7] + '/GaussFit')
		os.mkdir(dir_now[:-7] + '/GaussFit/PCH_gaussFit')
		os.mkdir(dir_now[:-7] + '/GaussFit/PCH_gaussFitPars')
		os.mkdir(dir_now[:-7] + '/GaussFit/PCH_modelSel')
	
	
	for file_now in files_to_analyze:
		print "================================================================================="
		print "I'm working on darta from: " + dir_now
		print "---------------------------------------------------------------------------------------"
		print file_now
		inputfile = dir_now + '/' + file_now
		# data
		f = open (inputfile, "rt")
		#print inputfile
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
		tmp_f3 = dir_now[:-7] + '/tmp/gaussFitTEX.txt'
		f3 = open(tmp_f3,'w')
		while l2:
			va = l2.split()
			th_u.append(float(va[0]))
			th_d.append(float(va[1]))
			l2 = f2.readline()

		## calculate histogram
		## now determine nice limits by hand:
		binwidth = 1
		xmax = np.max(np.fabs(time))
		ymax = np.max(np.fabs(counts))

		bins = np.arange(-0.5, ymax + binwidth-0.5, binwidth)
		hist_y, hist_x = np.histogram(counts, bins=bins)#, normed = True
		hist_x = bins + 0.5



		#find init pars
		first_x = hist_x[0:int(th_u[0])]
		first_y = hist_y[0:int(th_u[0])]
		second_x = hist_x[int(th_u[0]):]
		second_y = hist_y[int(th_u[0]):]

		#print th_u[0]
		#print first_x[:]
		#print second_x[0]
		#print max(first_y)
		#print max(second_y)
		##print hist_x.index(20.0)
		#print np.where(first_y==max(first_y))
		#print np.where(hist_y==max(second_y))

		a_f = max(first_y) #amplitude for first max
		o_f = np.where(first_y==max(first_y))[0][0] + 1 # offset for first max
		w_f = mt.sqrt(o_f) #width for first max
		#print array[tmp_where[0][0]]

		a_s = max(second_y)#amplitude for second max
		o_s = np.where(hist_y==max(second_y))[0][-1]# offset for second max
		w_s = mt.sqrt(o_s)#width  for first max

		#print a_f
		#print o_f
		#print w_f
		#print a_s
		#print o_s
		#print w_s

		# parameters: ampl, offset, width
		#AM	
		#limits = [0.,False, False]
		
		#Do we need lilitattions? [amplitude, offset, width]
		limitedmin_set=[True, False, True]
		
		AIC1 = []
		AIC2 = []
		print "2 Gauss fit"
		ngauss_now = 2
		#(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y,ngauss = ngauss_now, params=[a_f,o_f,w_f, a_s,o_s,w_s])
		(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y,ngauss = ngauss_now, params=[a_f,o_f,w_f, a_s,o_s,w_s], limitedmin=limitedmin_set*ngauss_now, quiet=False)
		chi_norm = chii/len(hist_y)
		#-----------------------------------------------------------------
		#!!!!!!!!!!!!!!!!!! I can also add errors!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		#-----------------------------------------------------------------
		f3.write('amplitude & ' + "%.3f" %(pars[0]) + '&' + "%.3f" %(pars[3]) + '&&&&&&& \\\\\n')
		f3.write('offset & ' + "%.3f" %(pars[1]) + '&' + "%.3f" %(pars[4]) + '&&&&&&& \\\\\n')
		f3.write('width & ' + "%.3f" %(pars[2]) + '&' + "%.3f" %(pars[5]) + '&&&&&&& \\\\\n')
		f3.write('chi&'+ "%.3f" %(chii) + '& chi norm&'+ "%.3f" %(chi_norm) + '&&&&& \\\\\\hline\n')
		#AIC model selection
		AIC1.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[0])
		AIC2.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[1])

		pl.subplot(2,4,1)
		plot_fits(hist_x, hist_y, n_gauss, pars)
		#f3.write(pars)
		#print "pars"
		#print pars
		#print "n_gauss"
		#print n_gauss
		#print "errors" 
		#print errors 
		#print "chii"
		#print chii/len(hist_y)


		
		print "3 Gauss fit"
		ngauss_now = 3
		o_m = (o_f+o_s)/2
		a_m = hist_y[int(o_m)]
		w_m = mt.sqrt(o_m)
		#print o_m
		#print a_m
		(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y,ngauss = ngauss_now, params=[a_f,o_f,w_f, a_m,o_m,w_m, a_s,o_s,w_s], limitedmin=limitedmin_set*ngauss_now, quiet=False)
		#(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y,ngauss = ngauss_now, params=[a_f,o_f,w_f, a_m,o_m,w_m, a_s,o_s,w_s], limitedmin=limits*ngauss_now)
		chi_norm = chii/len(hist_y)
		f3.write('amplitude & ' + "%.3f" %(pars[0]) + '&' + "%.3f" %(pars[3])  + '&' + "%.3f" %(pars[6]) + '&&&&&& \\\\\n')
		f3.write('offset & ' + "%.3f" %(pars[1]) + '&' + "%.3f" %(pars[4])  + '&' + "%.3f" %(pars[7]) + '&&&&&& \\\\\n')
		f3.write('width & ' + "%.3f" %(pars[2]) + '&' + "%.3f" %(pars[5])  + '&' + "%.3f" %(pars[8]) + '&&&&&& \\\\\n')
		f3.write('chi&'+ "%.3f" %(chii) + '& chi norm&'+ "%.3f" %(chi_norm) + '&&&&& \\\\\\hline\n')

		#AIC model selection
		AIC1.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[0])
		AIC2.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[1])

		pl.subplot(2,4,2)
		plot_fits(hist_x, hist_y, n_gauss, pars)


		#print  chii/len(counts)
		#print "pars"
		#print pars
		#print "n_gauss"
		#print n_gauss
		#print "errors" 
		#print errors 
		#print "chii"
		#print chii/len(hist_y)


		print "4 Gauss fit"
		ngauss_now = 4
		o_1 = int((o_f+o_s)/4)
		a_1 = hist_y[o_1]
		w_1 = mt.sqrt(o_1)
		o_2 = int(3*(o_f+o_s)/4)
		a_2 = hist_y[o_2]
		w_2 = mt.sqrt(o_2)
		#(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y, ngauss = ngauss_now, params=[a_f,o_f,w_f, a_1,o_1,w_1, a_2,o_2,w_2, a_s,o_s,w_s])
		(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y, ngauss = ngauss_now, params=[a_f,o_f,w_f, a_1,o_1,w_1, a_2,o_2,w_2, a_s,o_s,w_s], limitedmin=limitedmin_set*ngauss_now, quiet=False)
		chi_norm = chii/len(hist_y)
		f3.write('amplitude & ' + "%.3f" %(pars[0]) + '&' + "%.3f" %(pars[3])  + '&' + "%.3f" %(pars[6]) +  '&' + "%.3f" %(pars[9]) + '&&&&& \\\\\n')
		f3.write('offset & ' + "%.3f" %(pars[1]) + '&' + "%.3f" %(pars[4])  + '&' + "%.3f" %(pars[7]) +  '&' + "%.3f" %(pars[10]) + '&&&&& \\\\\n')
		f3.write('width & ' + "%.3f" %(pars[2]) + '&' + "%.3f" %(pars[5])  + '&' + "%.3f" %(pars[8]) +  '&' + "%.3f" %(pars[11]) + '&&&&& \\\\\n')
		f3.write('chi&'+ "%.3f" %(chii) + '& chi norm&'+ "%.3f" %(chi_norm) + '&&&&& \\\\\\hline\n')

		#AIC model selection
		AIC1.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[0])
		AIC2.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[1])

		pl.subplot(2,4,3)
		plot_fits(hist_x, hist_y, n_gauss, pars)
		#print  chii/len(counts)

		print "5 Gauss fit"
		ngauss_now = 5
		#minpars=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,]    
		#(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y,ngauss=ngauss_now, params=[a_f,o_f,w_f, a_1,o_1,w_1, a_m,o_m,w_m, a_2,o_2,w_2, a_s,o_s,w_s])
		(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y,ngauss=ngauss_now, params=[a_f,o_f,w_f, a_1,o_1,w_1, a_m,o_m,w_m, a_2,o_2,w_2, a_s,o_s,w_s], limitedmin=limitedmin_set*ngauss_now, quiet=False)
		chi_norm = chii/len(hist_y)
		f3.write('amplitude & ' + "%.3f" %(pars[0]) + '&' + "%.3f" %(pars[3])  + '&' + "%.3f" %(pars[6]) +  '&' + "%.3f" %(pars[9]) +  '&' + "%.3f" %(pars[12]) + '&&&& \\\\\n')
		f3.write('offset & ' + "%.3f" %(pars[1]) + '&' + "%.3f" %(pars[4])  + '&' + "%.3f" %(pars[7]) +  '&' + "%.3f" %(pars[10]) +  '&' + "%.3f" %(pars[3]) + '&&&& \\\\\n')
		f3.write('width & ' + "%.3f" %(pars[2]) + '&' + "%.3f" %(pars[5])  + '&' + "%.3f" %(pars[8]) +  '&' + "%.3f" %(pars[11]) +  '&' + "%.3f" %(pars[14]) + '&&&& \\\\\n')
		f3.write('chi&'+ "%.3f" %(chii) + '& chi norm&'+ "%.3f" %(chi_norm) + '&&&&&\\\\\\hline\n')

		#AIC model selection
		AIC1.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[0])
		AIC2.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[1])

		pl.subplot(2,4,4)
		plot_fits(hist_x, hist_y, n_gauss, pars)
		#print  chii/len(counts)

		print "6 Gauss fit"
		ngauss_now = 6
		#o_0 = int(o_f*0.8)
		o_0 = o_f*0.8    
		a_0 = hist_y[int(o_0)]
		w_0 = mt.sqrt(o_0)
		(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y,ngauss=ngauss_now, params=[a_0,o_0,w_0, a_f,o_f,w_f, a_1,o_1,w_1, a_m,o_m,w_m, a_2,o_2,w_2, a_s,o_s,w_s], limitedmin=limitedmin_set*ngauss_now, quiet=False)
		#(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y,ngauss=ngauss_now, params=[a_0,o_0,w_0, a_f,o_f,w_f, a_1,o_1,w_1, a_m,o_m,w_m, a_2,o_2,w_2, a_s,o_s,w_s], limitedmin=limits*ngauss_now)
		chi_norm = chii/len(hist_y)

		f3.write('amplitude & ' + "%.3f" %(pars[0]) + '&' + "%.3f" %(pars[3])  + '&' + "%.3f" %(pars[6]) +  '&' + "%.3f" %(pars[9]) +  '&' + "%.3f" %(pars[12]) +  '&' + "%.3f" %(pars[15]) + '&&&\\\\\n')	
		f3.write('offset & ' + "%.3f" %(pars[1]) + '&' + "%.3f" %(pars[4])  + '&' + "%.3f" %(pars[7]) +  '&' + "%.3f" %(pars[10]) +  '&' + "%.3f" %(pars[3]) +  '&' + "%.3f" %(pars[16]) + '&&&\\\\\n')
		f3.write('width & ' + "%.3f" %(pars[2]) + '&' + "%.3f" %(pars[5])  + '&' + "%.3f" %(pars[8]) +  '&' + "%.3f" %(pars[11]) +  '&' + "%.3f" %(pars[14]) +  '&' + "%.3f" %(pars[17]) + '&&&\\\\\n')
		f3.write('chi&'+ "%.3f" %(chii) + '& chi norm&'+ "%.3f" %(chi_norm) + '&&&&& \\\\\\hline\n')

		#AIC model selection
		AIC1.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[0])
		AIC2.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[1])

		pl.subplot(2,4,5)
		plot_fits(hist_x, hist_y, n_gauss, pars)

		print "7 Gauss fit"
		ngauss_now = 7
		#o_n = int(o_s*1.2)
		o_n = o_s*1.2
		a_n = hist_y[int(o_n)]
		w_n = mt.sqrt(o_n)
		(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y,ngauss=ngauss_now, params=[a_0,o_0,w_0, a_f,o_f,w_f, a_1,o_1,w_1, a_m,o_m,w_m, a_2,o_2,w_2, a_s,o_s,w_s, a_n,o_n,w_n], limitedmin=limitedmin_set*ngauss_now, quiet=False)
		#(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y,ngauss=ngauss_now, params=[a_0,o_0,w_0, a_f,o_f,w_f, a_1,o_1,w_1, a_m,o_m,w_m, a_2,o_2,w_2, a_s,o_s,w_s], limitedmin=limits*ngauss_now)
		chi_norm = chii/len(hist_y)

		f3.write('amplitude & ' + "%.3f" %(pars[0]) + '&' + "%.3f" %(pars[3])  + '&' + "%.3f" %(pars[6]) +  '&' + "%.3f" %(pars[9]) +  '&' + "%.3f" %(pars[12]) +  '&' + "%.3f" %(pars[15]) +  '&' + "%.3f" %(pars[18]) + '&&\\\\\n')	
		f3.write('offset & ' + "%.3f" %(pars[1]) + '&' + "%.3f" %(pars[4])  + '&' + "%.3f" %(pars[7]) +  '&' + "%.3f" %(pars[10]) +  '&' + "%.3f" %(pars[3]) +  '&' + "%.3f" %(pars[16]) +  '&' + "%.3f" %(pars[19]) + '&&\\\\\n')
		f3.write('width & ' + "%.3f" %(pars[2]) + '&' + "%.3f" %(pars[5])  + '&' + "%.3f" %(pars[8]) +  '&' + "%.3f" %(pars[11]) +  '&' + "%.3f" %(pars[14]) +  '&' + "%.3f" %(pars[17]) +  '&' + "%.3f" %(pars[20]) + '&&\\\\\n')
		f3.write('chi&'+ "%.3f" %(chii) + '& chi norm&'+ "%.3f" %(chi_norm) + '&&&&& \\\\\\hline\n')

		#AIC model selection
		AIC1.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[0])
		AIC2.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[1])

		pl.subplot(2,4,6)
		plot_fits(hist_x, hist_y, n_gauss, pars)

		print "8 Gauss fit"
		ngauss_now = 8
		#o_5_1 = int((o_f+o_s)/5)
		o_5_1 = (o_f+o_s)/5.0
		a_5_1 = hist_y[int(o_5_1)]
		w_5_1 = mt.sqrt(o_5_1)
		#o_5_2 = int(2*(o_f+o_s)/5)
		o_5_2 = 2.0*(o_f+o_s)/5
		a_5_2 = hist_y[int(o_5_2)]
		w_5_2 = mt.sqrt(o_5_2)
		#o_5_3 = int(3*(o_f+o_s)/5)
		o_5_3 = 3.0*(o_f+o_s)/5
		a_5_3 = hist_y[int(o_5_3)]
		w_5_3 = mt.sqrt(o_5_3)
		#o_5_4 = int(4*(o_f+o_s)/5)
		o_5_4 = 4.0*(o_f+o_s)/5
		a_5_4 = hist_y[int(o_5_4)]
		w_5_4 = mt.sqrt(o_5_4)
		(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y,ngauss=ngauss_now, params=[a_0,o_0,w_0, a_f,o_f,w_f, a_5_1,o_5_1,w_5_1, a_5_2,o_5_2,w_5_2, a_5_3,o_5_3,w_5_3, a_5_4,o_5_4,w_5_4, a_s,o_s,w_s, a_n,o_n,w_n], limitedmin=limitedmin_set*ngauss_now, quiet=False)
		#(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y,ngauss=ngauss_now, params=[a_0,o_0,w_0, a_f,o_f,w_f, a_1,o_1,w_1, a_m,o_m,w_m, a_2,o_2,w_2, a_s,o_s,w_s], limitedmin=limits*ngauss_now)
		chi_norm = chii/len(hist_y)

		f3.write('amplitude & ' + "%.3f" %(pars[0]) + '&' + "%.3f" %(pars[3])  + '&' + "%.3f" %(pars[6]) +  '&' + "%.3f" %(pars[9]) +  '&' + "%.3f" %(pars[12]) +  '&' + "%.3f" %(pars[15]) +  '&' + "%.3f" %(pars[18]) +  '&' + "%.3f" %(pars[21]) + '& \\\\\n')	
		f3.write('offset & ' + "%.3f" %(pars[1]) + '&' + "%.3f" %(pars[4])  + '&' + "%.3f" %(pars[7]) +  '&' + "%.3f" %(pars[10]) +  '&' + "%.3f" %(pars[3]) +  '&' + "%.3f" %(pars[16]) +  '&' + "%.3f" %(pars[19]) +  '&' + "%.3f" %(pars[22]) + '& \\\\\n')
		f3.write('width & ' + "%.3f" %(pars[2]) + '&' + "%.3f" %(pars[5])  + '&' + "%.3f" %(pars[8]) +  '&' + "%.3f" %(pars[11]) +  '&' + "%.3f" %(pars[14]) +  '&' + "%.3f" %(pars[17]) +  '&' + "%.3f" %(pars[20]) +  '&' + "%.3f" %(pars[23]) + '& \\\\\n')
		f3.write('chi&'+ "%.3f" %(chii) + '& chi norm&'+ "%.3f" %(chi_norm) + '&&&&& \\\\\\hline\n')

		#AIC model selection
		AIC1.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[0])
		AIC2.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[1])

		pl.subplot(2,4,7)
		plot_fits(hist_x, hist_y, n_gauss, pars)

		print "9 Gauss fit"
		ngauss_now = 9
		o_s2 = o_s*1.1
		a_s2 = hist_y[int(o_s2)]
		w_s2 = mt.sqrt(o_s2)
		(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y,ngauss=ngauss_now, params=[a_0,o_0,w_0, a_f,o_f,w_f, a_5_1,o_5_1,w_5_1, a_5_2,o_5_2,w_5_2, a_5_3,o_5_3,w_5_3, a_5_4,o_5_4,w_5_4, a_s,o_s,w_s, a_s2,o_s2,w_s2, a_n,o_n,w_n], limitedmin=limitedmin_set*ngauss_now, quiet=False)
		#(pars,n_gauss,errors, chii) =  multigaussfit(hist_x[:-1],hist_y,ngauss=ngauss_now, params=[a_0,o_0,w_0, a_f,o_f,w_f, a_1,o_1,w_1, a_m,o_m,w_m, a_2,o_2,w_2, a_s,o_s,w_s], limitedmin=limits*ngauss_now)
		chi_norm = chii/len(hist_y)

		f3.write('amplitude & ' + "%.3f" %(pars[0]) + '&' + "%.3f" %(pars[3])  + '&' + "%.3f" %(pars[6]) +  '&' + "%.3f" %(pars[9]) +  '&' + "%.3f" %(pars[12]) +  '&' + "%.3f" %(pars[15]) +  '&' + "%.3f" %(pars[18]) +  '&' + "%.3f" %(pars[21]) +  '&' + "%.3f" %(pars[24]) + '\\\\\n')	
		f3.write('offset & ' + "%.3f" %(pars[1]) + '&' + "%.3f" %(pars[4])  + '&' + "%.3f" %(pars[7]) +  '&' + "%.3f" %(pars[10]) +  '&' + "%.3f" %(pars[3]) +  '&' + "%.3f" %(pars[16]) +  '&' + "%.3f" %(pars[19]) +  '&' + "%.3f" %(pars[22]) +  '&' + "%.3f" %(pars[25]) + '\\\\\n')
		f3.write('width & ' + "%.3f" %(pars[2]) + '&' + "%.3f" %(pars[5])  + '&' + "%.3f" %(pars[8]) +  '&' + "%.3f" %(pars[11]) +  '&' + "%.3f" %(pars[14]) +  '&' + "%.3f" %(pars[17]) +  '&' + "%.3f" %(pars[20]) +  '&' + "%.3f" %(pars[23]) +  '&' + "%.3f" %(pars[26]) + '\\\\\n')
		f3.write('chi&'+ "%.3f" %(chii) + '& chi norm&'+ "%.3f" %(chi_norm) + '&&&&& \\\\\\hline\n')

		#AIC model selection
		AIC1.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[0])
		AIC2.append(calc_AIC(chii, ngauss_now*3, len(hist_y) )[1])

		pl.subplot(2,4,8)
		plot_fits(hist_x, hist_y, n_gauss, pars)

		#print  chii/len(counts)
		minAIC1 = min(AIC1)
		for l in range(len(AIC1)):
			if (AIC1[l]==minAIC1):
				AIC1[l] = '\\textbf{' + "%.2f" % (AIC1[l])+'}'
			else:
				AIC1[l] = "%.2f" % (AIC1[l])

		AIC1_f = open(dir_now[:-7] + '/GaussFit/PCH_modelSel/aic1Gauss_TEX.txt','a')
		AIC1_f.write('$' + inputfile[12:] +'$&'+AIC1[0]+'&'+AIC1[1]+'&'+AIC1[2]+'&'+AIC1[3]+'&'+AIC1[4]+'&'+AIC1[5]+'&'+AIC1[6]+'&'+AIC1[7]+'\\\\\n')#
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
		
		outputfile = dir_now[:-7] + '/GaussFit/PCH_gaussFit/PCH_gaussFit' + file_now[0:-4] + '.png'
		pl.savefig(outputfile)
		f3.closed
		pl.clf()
