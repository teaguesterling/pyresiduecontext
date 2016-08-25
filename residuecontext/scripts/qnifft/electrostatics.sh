#!/bin/bash
USAGE="$0 <PDB_FILE> [<OUTPUT_FILE>]"

PDB="$1"
PDB_NAME="$( basename "$PDB" .pdb )"
PARAM_FILE="${PDB_NAME}-qnifft.parm"

SIZE=193
QDIR="$( dirname "$0" )"
DELDEF="${QDIR}/delphi.def"
CHARGE="${QDIR}/amb.crg.oxt"
RADIUS="${QDIR}/default.siz"

QNIFFT="${QDIR}/qnifft22_${SIZE}"

# Delphi Defaults
cat "$DELDEF" > "${PARAM_FILE}"

# DOCK37 Parameters
echo "grid=${SIZE}" >> "${PARAM_FILE}"
echo "charge=${CHARGE}" >> "${PARAM_FILE}"
echo "radius=${RADIUS}" >> "${PARAM_FILE}"
echo "pdb_input=${PDB}" >> "${PARAM_FILE}"
echo "pdb_output_file=${PDB_NAME}.atm" >> "${PARAM_FILE}"
echo "phi_output_file=${PDB_NAME}.phi" >> "${PARAM_FILE}"
echo "border=15" >> "${PARAM_FILE}"
echo "sizing=border border_solvent=10." >> "${PARAM_FILE}"

${QNIFFT} "${PARAM_FILE}"
