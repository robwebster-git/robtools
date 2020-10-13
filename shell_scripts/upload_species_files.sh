#!/bin/bash

AOI=$1

for TILE in ${@:2}
do
   for SPECIES in 111 691 721 802 820  
     do
       FILEPATH=$MAINDIR/Projects/${AOI}_${TILE}/Scores/*${SPECIES}_pc/
       aws s3 cp $FILEPATH s3://gsi-test-storage/Projects/FC_CL/species/local_species_pertile/${SPECIES}/ --recursive --exclude "*.*" --include "*pc_adj.*.tif"
     done
   echo "Finished Tile $TILE"
done
echo "Finished all"
