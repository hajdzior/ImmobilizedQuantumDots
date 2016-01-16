#!/bin/bash

cd trajsTXT
numfiles=(*)
i=0
echo "..Plotting trajectories"
for trajectory in *.txt
do
        file=${trajectory:0:${#trajectory}-4}
        ((i = $i+1))
        echo ".. I am working on file #" $i "/" ${#numfiles[@]}":" $trajectory
        python ../SOFT_calcTrajs-onOffTimes/plotTraj-details.py -i $trajectory -j ../thresh/$file"-threshVAL.txt" -o $file".png" # dokl. traj podzielona na strony
        python ../SOFT_calcTrajs-onOffTimes/plotTraj-summary.py -i $trajectory -j ../thresh/$file"-threshVAL.txt" -o $file".png" # cala traj + hist
done
cd ..
