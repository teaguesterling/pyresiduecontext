#!/bin/bash

REDUCEDIR="$( dirname "$0" )"
PDB="$1"
NOH="${PDB}.noh"
PRE="${PDB}.prepdb"
TMP="${PDB}.crgtmp"
CRG="${2-${PDB}.crg}"

# DEFAULT DOCK37 Reduce args
RHSDICT="${REDUCEDIR}/reduce_wwPDB_het_dict.txt"
RARGS=( -OH -HIS -ALLALT -ROTNH3 -Keep -METALBump1.5 -NONMETALBump-5.0 )

grep -v '^HETATM' "${PDB}" > "${NOH}"
${REDUCEDIR}/cleanpdb.py "${NOH}" "${PRE}" inital
${REDUCEDIR}/reduce -db "${RHSDICT}" "${RARGS[@]}" "${PRE}" > "${TMP}"
${REDUCEDIR}/cleanpdb.py "${TMP}" "${CRG}" final

