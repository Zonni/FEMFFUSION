#!/usr/bin/env bash 

# Script to run all the examples in 1D_ROM_td.prm
# Must be run for FEMFFUSION MAIN dir as
# sh 1D_ROM_1g/run_all.sh

./femffusion.exe -f 1D_ROM/1D_ROM_td.prm
./femffusion.exe -f 1D_ROM/1D_ROM_pod10_dyn.prm
./femffusion.exe -f 1D_ROM/1D_ROM_pod_sta10.prm
./femffusion.exe -f 1D_ROM/1D_ROM_pod_sta20.prm
./femffusion.exe -f 1D_ROM/1D_ROM_pod_sta10_LUPOD.prm
./femffusion.exe -f 1D_ROM/1D_ROM_pod_sta20_LUPOD.prm
./femffusion.exe -f 1D_ROM/1D_ROM_rampmat12_10.prm
