#!/bin/bash

#------------------------
# This script prepres the plots and fits of the ON and OFF times distributions
# As the parmeterr the input file (ON or OFF times) are given

# input folder: 'on_off_times'
# output folder: 'on_off_plots's
#------------------------

i=0
cd on_off_times
for filename in *0-001*
do
	file=${filename:0:-4}
	((i = $i+1))
	echo "I am working on file #" $i ": " $filename
	echo "Analysis: power-law" 
	python ../SOFT_on-off-Plenz/on-off-Plenz.py -i ../on_off_times/$filename
	cp ../tmp/onPlenz.png ../on_off_analysis/on_off_plots_Plenz/powerLaw/$file'-Plenz-on.png'
	cp ../tmp/offPlenz.png ../on_off_analysis/on_off_plots_Plenz/powerLaw/$file'-Plenz-off.png'
	echo "Analysis: truncated power-law" 
	python ../SOFT_on-off-Plenz/on-off-Plenz_v3truncPL.py -i ../on_off_times/$filename
	cp ../tmp/onPlenz.png ../on_off_analysis/on_off_plots_Plenz/powerLaw_trunc/$file'_Plenz_on_trunc.png'
	cp ../tmp/offPlenz.png ../on_off_analysis/on_off_plots_Plenz/powerLaw_trunc/$file'_Plenz_off_trunc.png'
	cp ../tmp/Fig_on_PL_truncPL.png ../on_off_analysis/on_off_plots_Plenz/powerLaw_trunc/$file'_Plenz_on_PL_truncPL.png'
	cp ../tmp/Fig_off_PL_truncPL.png ../on_off_analysis/on_off_plots_Plenz/powerLaw_trunc/$file'_Plenz_off_PL_truncPL.png'
done

cd ..
rm tmp/*

echo "Done, are you glad? Have a nice day!"

