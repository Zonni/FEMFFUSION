#!/bin/bash
set -e  # Stop at first error


./femffusion.exe -f 2D_C5G7_DSPN/c5g7_SDP1.prm
./femffusion.exe -f 2D_C5G7_DSPN/c5g7_SDP2.prm 
./femffusion.exe -f 2D_C5G7_DSPN/c5g7_SDP3.prm 

./femffusion.exe -f 2D_C5G7_DSPN/c5g7_SP1.prm
./femffusion.exe -f 2D_C5G7_DSPN/c5g7_SP3.prm 
./femffusion.exe -f 2D_C5G7_DSPN/c5g7_SP5.prm 
./femffusion.exe -f 2D_C5G7_DSPN/c5g7_SP7.prm 

./femffusion.exe -f 2D_C5G7_DSPN/c5g7_SDP1_Nazari.prm
./femffusion.exe -f 2D_C5G7_DSPN/c5g7_SDP2_Nazari.prm 
./femffusion.exe -f 2D_C5G7_DSPN/c5g7_SDP3_Nazari.prm 