#!/bin/bash
#Run from FEMFFUSION_DIR
set -e  # Stop at first error

./femffusion.exe -f 3D_ACP100/ACP100_FE2_td_FOM.prm
./femffusion.exe -f 3D_ACP100/ACP100_FE2_td_ROM5.prm
./femffusion.exe -f 3D_ACP100/ACP100_FE2_td_ROM10.prm
./femffusion.exe -f 3D_ACP100/ACP100_FE2_td_ROM25.prm
./femffusion.exe -f 3D_ACP100/ACP100_FE2_td_ROM50.prm
./femffusion.exe -f 3D_ACP100/ACP100_FE2_td_ROM100.prm