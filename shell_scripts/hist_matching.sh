#!/bin/bash

# First argument $1 is text file containing full paths to reference files in reference year
if [ -z "$1" ]; then
   echo   
   echo "USAGE: bash hist_matching.sh  [TEXT_FILE_OF_INPUT_FOURIER_FILES]  [REFERENCE_YEAR_YY]  [TARGET_YEAR_YY]  [AOI_NAME]  [TILE_ID]"
   echo
   exit 
fi


# Histograms to use
reference_year=$2

# Rebalance the histograms from this year with those from the reference year above
target_year=$3

aoi=$4
tile=$5

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
   NewYearFile="/home/nx06/nx06/shared/RF/S2AWS/${aoi}${target_year}_${tile}/Fmask_Fourier/${fn2}"
   CorrNewYearFile="/home/nx06/nx06/shared/RF/S2AWS/${aoi}${target_year}_${tile}/Fmask_Fourier/hist_matched_${fn2}"
   Histogram $NewYearFile NewHistogram.txt RefHistogram.txt $CorrNewYearFile
done < $1
