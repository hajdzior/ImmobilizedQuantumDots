#!/bin/bash

#------------------------
# This script prepres the plots and fits of the ON and OFF times distributions
# As the parmeterr the input file (ON or OFF times) are given

# input folder: 'on_off_times'
# output folder: 'on_off_plots's
#------------------------

i=0

cd trajsTXT
for filename in *0-001*
do
	file=${filename:0:-4}
	((i = $i+1))
	echo "I am working on file #" $i ": " $filename
	python ../SOFT_GaussFit/fit_gauss_v2-1autom.py -i ../trajsTXT/$filename -j ../thresh/$file"-threshVAL.txt" -o '../GaussFit/PCH_gaussFit/'$file'_PCHgaussFit.png'
	cp ../tmp/gaussFitTEX.txt ../GaussFit/PCH_gaussFitPars/$file'-GaussFitTEX.txt'
done

cd ..
# rm tmp/*

echo "Done, are you glad? Have a nice day!"

