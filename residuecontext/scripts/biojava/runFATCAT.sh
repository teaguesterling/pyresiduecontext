#!/bin/bash

# examples:

#  bash runFATCAT.sh -pdb1 4hhb.A -pdb2 4hhb.B -pdbFilePath /tmp/ -autoFetch -show3d 
#  bash runFATCAT.sh -pdb1 4hhb.A -pdb2 4hhb.B -pdbFilePath /tmp/ -printXML 
#  bash runFATCAT.sh -pdb1 4hhb.A -pdb2 4hhb.B -pdbFilePath /tmp/ -printFatCat
# 
#  run a flexible alignment:
#  bash runFATCAT.sh -pdb1 4hhb.A -pdb2 4hhb.B -pdbFilePath /tmp/ -printFatCat -flexible 
# 
# for more examples see runCE.sh
# jCE works almost exactly the same...

# send the arguments to the java app
args="$*"

function scriptdir {
    cd "$(dirname "$1")"
    cd "$(dirname "$(readlink "$1" 2>/dev/null || basename "$1" )")"
    pwd 
}
DIR="$(scriptdir "$0" )"

java -Xmx500M -cp "$DIR/jars/*" org.biojava.bio.structure.align.fatcat.FatCat $args
