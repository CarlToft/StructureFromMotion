#!/bin/bash
#Usage: ./start_matlab <path_to_lib_folder>
set -x
SFM_LIB_PATH=${1:-"$HOME"}
MOSEK_PATH="${SFM_LIB_PATH}/mosek/8/tools/platform/linux64x86/bin"
SFM_LIB_PATH=${SFM_LIB_PATH} LD_LIBRARY_PATH=${MOSEK_PATH} matlab

