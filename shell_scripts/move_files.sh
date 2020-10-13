#!/bin/bash
CHUNK=$1

#  The rest of the arguments should be S2 tile IDs eg 17SMS 17SMR 18SUE etc

for TILE in ${@:2}
do
  SOURCE=$MAINDIR/Projects/FC_CL_${TILE}/Scores/
  DEST=$MAINDIR/Projects/FC_CL_${CHUNK}/
  mv $SOURCE/*.tif $DEST 
done
