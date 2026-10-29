#!/bin/bash
set -e  # Stop at first error

echo 'DYNAMIC'
./femffusion.exe -f 3D_Langenbuch/3D_Langenbuch_fe3_rom_dynamic.prm 
echo 'static_bank_swept.'
./femffusion.exe -f 3D_Langenbuch/3D_Langenbuch_fe3_rom_static_bank_swept.prm
echo 'static_time_sampled'
./femffusion.exe -f 3D_Langenbuch/3D_Langenbuch_fe3_rom_static_time_sampled.prm
echo '_bank_ihs.'
./femffusion.exe -f 3D_Langenbuch/3D_Langenbuch_fe3_rom_static_bank_ihs.prm
