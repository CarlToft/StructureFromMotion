#!/bin/bash
#Usage: ./start_matlab <path_to_lib_folder>
SFM_LIB_PATH=${1:-"$HOME"}
SFM_ABS_LIB_PATH="$( readlink -f ${SFM_LIB_PATH} )"
MOSEK_PATHS=( ${SFM_ABS_LIB_PATH}/mosek/*/tools/platform/linux64x86/bin )
MOSEK_FOLDER="${MOSEK_PATHS[0]}"

echo \
"Library path: ${SFM_ABS_LIB_PATH}
Mosek path: ${MOSEK_FOLDER}"

SFM_LIB_PATH=${SFM_ABS_LIB_PATH} LD_LIBRARY_PATH=${MOSEK_FOLDER} matlab
