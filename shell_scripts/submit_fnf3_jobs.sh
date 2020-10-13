#!/bin/bash

#  Set AOI as first commandline argument

AOI=$1

#  The rest of the arguments should be S2 tile IDs eg 17SMS 17SMR 18SUE etc

for TILE in ${@:2}
do
  MAXHEIGHTDIR=$MAINDIR/Projects/${AOI}_${TILE}/Scores/tree_max_height
  CCDIR=$MAINDIR/Projects/${AOI}_${TILE}/Scores/tree_canopy_cover
  FNF2DIR=$MAINDIR/Projects/${AOI}/FNF/FNF2/Per_Tile  
  OUTDIR=$MAINDIR/Projects/${AOI}/FNF/FNF3/Per_Tile
  #OUTDIR=~/FNF3
  CODEDIR=/home/nx06/nx06/rwgsi/git/ms-core/src/moonshine/util
  python $CODEDIR/create_fnf3.py --tree_height $MAXHEIGHTDIR --canopy_cover $CCDIR --fnf2 $FNF2DIR --outdir $OUTDIR
done
