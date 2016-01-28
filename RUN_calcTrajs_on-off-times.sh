#!/bin/bash

#---------------------------
# Last mod: 11.07.14

# This script clcultes the trjectories with the given bin size (in seconds)

# Executes the C++ progrmme "calcTrjs", as a parameter the trajectory bin size has to be given!
# The C++ programme is based on the one downloded from PicoQunt, in the output file some aditional ('string') 
#   information are given, the file hs to be "clened" to remove the text from the pure trajectory

# The final file is two-column:
# time & counts
#input folder: 'ptu'
#output folder: 'trjsTXT'

# At the end the python script ('plot.pl') is executed to drw the fluorescence trjectory and save it to the file
# As parameters, the name of the trajectory file (input) nd the name of the plot file (output) are given. 
#input folder: 'trjsTXT'
#output folder: 'plots'

# MW 11.07.2014
#----------------------------

# GET CURRENT SOFTWARE VERSION

SOFT_VERSION=`git show --format=%ci | head -n 1`
echo "Calculating with software from $SOFT_VERSION."
  

set -e #errors
array=(0.001)
names=("0-001")

# array=(0.1 0.01 0.001)
# names=("0-1" "0-01" "0-001")

# array=(0.0001 0.00001 )
# names=("0-0001" "0-00001")

#----------------------------------------
# Calculate trajectories from .ptu
#----------------------------------------
cd ptu
numfiles=(*)
i=0
echo "..Calculating the trajectories"
for filename in *.ptu
do
	file=${filename:0:-4}
	((i = $i+1))
	echo ".. I am working on file #" $i "/" ${#numfiles[@]}":" $filename
	for ix in ${!array[*]}
	do
		printf "  ..%s\n" "${array[$ix]}"
		../SOFT_calcTrajs-onOffTimes/ReadIATandCalcTrajs $file.ptu ../trajsTXT/$file--"${names[$ix]}".txt ${array[$ix]}
	done
done	
cd ..

#--------------------------
# Calculate on-off times
#--------------------------
cd trajsTXT
numfiles=(*)
i=0
echo "..Calculating on-off times"
for trajectory in *.txt
do
	file=${trajectory:0:-4}
	((i = $i+1))
	echo ".. I am working on file #" $i "/" ${#numfiles[@]}":" $trajectory
	cp $trajectory ../tmp/trajectory.txt
	
	scilab -nwni -f ../SOFT_calcTrajs-onOffTimes/ON_OFF_times.sce
	cp ../tmp/on-off_times_thresh-1.txt ../on_off_times/$file-on-off_times_thresh-1.txt
	cp ../tmp/on-off_times_thresh-2--1-3.txt ../on_off_times/$file-on-off_times_thresh-2--1-3.txt
	cp ../tmp/on-off_times_thresh-2--1-5.txt ../on_off_times/$file-on-off_times_thresh-2--1-5.txt
	cp ../tmp/thresholds.txt ../thresh/$file-threshTEX.txt
	cp ../tmp/thresholds_val.txt ../thresh/$file-threshVAL.txt
done
cd ..

#--------------------------
# Plot trajs
#--------------------------

./RUN_plot_traj.sh

#--------------------------
# 2D hist of intensities
#--------------------------
cd trajsTXT
numfiles=(*)
i=0
echo "..Calculating 2d intensity histograms "
for filename in *.txt
do
	file=${filename:0:-4}
	((i = $i+1))
	echo ".. I am working on file #" $i "/" ${#numfiles[@]}":" $trajectory
	python ../SOFT_calcTrajs-onOffTimes/2dHistInt.py -i $filename
	cp ../tmp/int2dHist.png ../2D-int_hist/$file-int2dHist.png
done
cd ..

rm tmp/*

echo $SOFT_VERSION > trajsPLOT/soft_version.txt
git rev-parse HEAD >> trajsPLOT/soft_version.txt

echo "Done, are you glad? Have a nice day!"
