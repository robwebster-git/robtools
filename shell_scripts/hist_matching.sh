#!/bin/bash

# First argument $1 is text file containing full paths to reference files in reference year
if [ -z "$1" ]; then
   echo   
   echo "Please supply a text file of input reference TIFs"
   echo
   exit 
fi

if [ -z "$2" ]; then
   echo   
   echo "Please supply a reference year in 2 digits eg 19 for 2019"
   echo
   exit
fi

if [ -z "$3" ]; then
   echo
   echo "Please supply a target year in 2 digits eg 16 for 2016"
   echo
   exit
fi


# Histograms to use
reference_year=$2

# Rebalance the histograms from this year with those from the reference year above
target_year=$3

echo "Using reference year 20${reference_year}"
echo
echo "Matching files for target year 20${target_year} with histogram values from reference year"
echo
sleep 2

while read f
do
   Histogram ${f} RefHistogram.txt
   fn=$(basename $f)
   fn2=${fn/$reference_year/$target_year}
   NewYearFile="/home/nx06/nx06/shared/RF/S2AWS/OA_GE${target_year}_36MZE/Fmask_Fourier/${fn2}"
   CorrNewYearFile="/home/nx06/nx06/shared/RF/S2AWS/OA_GE${target_year}_36MZE/Fmask_Fourier/hist_corrected_${fn2}"
   Histogram $NewYearFile NewHistogram.txt RefHistogram.txt $CorrNewYearFile
done < $1
