#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 target_dir"
fi

TARGET=$1
CURRENT=`pwd`

LINKS=".git SOFT_calcTrajs-onOffTimes RUN_calcTrajs_on-off-times.sh RUN_plot_traj.sh"
DIRS="trajsTXT trajsPLOT thresh tmp ptu on_off_analysis on_off_analysis/on_off_plots_Clausset on_off_analysis/on_off_plots_Plenz on_off_times GaussFit 2D-int_hist"

mkdir $TARGET

for directory in $DIRS
do
  mkdir $TARGET/$directory
done

for link in $LINKS
do
  ln -s $CURRENT/$link $TARGET/$link
done
