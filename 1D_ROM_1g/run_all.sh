#!/usr/bin/env bash 

# Script to run all the examples in 1D_ROM_td.prm
# Must be run for FEMFFUSION MAIN dir as
# sh 1D_ROM_1g/run_all.sh

./femffusion.exe -f 1D_ROM_1g/1D_ROM_td.prm
./femffusion.exe -f 1D_ROM_1g/1D_ROM_pod_dyn3.prm
./femffusion.exe -f 1D_ROM_1g/1D_ROM_pod_dyn5.prm
./femffusion.exe -f 1D_ROM_1g/1D_ROM_pod_dyn3_LUPOD.prm
./femffusion.exe -f 1D_ROM_1g/1D_ROM_pod_dyn5_LUPOD.prm
./femffusion.exe -f 1D_ROM_1g/1D_ROM_pod_dyn11_LUPOD.prm


