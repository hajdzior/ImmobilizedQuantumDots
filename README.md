# README - software
## Data analysis (c) 2016

# RUN_calcTrajs_on-off-times.sh

### Files in: SOFT_calcTrajs-onOffTimes
	* ReadIATandCalcTrajs_v0_2	->	C++ programme
		# translates the binary file into ASCII and calculates the trajectory with a given bin size (in seconds)
		# input:	ptu/result_file.ptu
		# output:	trajsTXT/*.txt
	
	* ON_OFF_times_v2_1.sce		->	Scilab script 
		# calculates thresholds and on and off times
		# input: 	trajsTXT/*.txt
		# output:	thresh/*.txt, on_off_times/*.txt

	* plotTraj-summary.py			->	Python script
		# plots complete trajectory (low resolution), photon counting histograms and thresholds
		# input:	trajsTXT/*.txt, thresh/*.txt
		# output:	trajsPLOT/*.png

	* plotTraj-details.py			->	Python script
		# plots detailed trajectories and thresholds
		# input:	trajsTXT/*.txt, thresh/*.txt
		# output:	trajsPLOT/*.png


	* 2dHistInt.py				->	Python script
		# plots 2D intensity histograms
		# input:	trajsTXT/*.txt
		# output:	2D-int_hist/*.png


**************************************************
   RUN_FitPlot-on-off-times_Clausset_v0_1.sh
**************************************************
Files in:	SOFT_on-off-Clausset


**************************************************
   RUN_FitPlot-on-off-times_Plenz_v0_1.sh
**************************************************
Files in:	SOFT_on-off-Plenz

	* powerlaw.py				->	the powerlaw Package 
		# used in on-off-Plenz.py and on-off-Plenz_v3truncPL.py
			> 'import powerlaw' in line 8

	* on-off-Plenz.py			->	Python script
		# on-off analysis, power law fits
		# input:	on_off_times/*.txt
		# output:	on_off_analysis/on_off_plots_Plenz/powerLaw/*.png

	* on-off-Plenz_v3truncPL.py	->	Python script
		# on-off analysis, truncated power law fits + model selection (power-law and truncated power law)
		# input:	on_off_times/*.txt
		# output:	on_off_analysis/on_off_plots_Plenz/powerLaw_trunc/*_trunc.png -> fit
					on_off_analysis/on_off_plots_Plenz/powerLaw_trunc/*__PL_truncPL.png -> model selection
	

**************************************************
                 RUN_gaussFit-v2.sh
**************************************************
Files in:	SOFT_GaussPoissFit
	Fits 2, 3, 4, 5, 6, ... Gauss distributions to the PCH
	* mpfit.py					->	Package/class 
		# mpfit used in gaussfitter (line 18 of gaussfitter_v0_1.py

	* gaussfitter_v0_1.py		->	Function definitions
		# n_gaussian and multigaussfit used in fit_gauss_v2-1autom.py

	* fit_gauss_v2-1autom.py	->	Python script
		# Fits gauss distributions and plots the result
		# input:	trajsTXT/*.txt, thresh/*.txt
		# output:	GaussFit/PCH_gaussFitPars/*.png, PCH_gaussFitPars/*.txt, PCH_modelSel/*.txt
		
**************************************************
                 RUN_gaussPoissFit-v2.sh
**************************************************
Files in:	SOFT_GaussPoissFit


