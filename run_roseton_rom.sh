#!/bin/bash
set -e  # Stop at first error

echo 'DYNAMIC'
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_dyn_ns6.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_dyn_ns10.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_dyn_ns16.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_dyn_ns20.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_dyn_ns40.prm
echo 'Static time-sampled'
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ts_ns6.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ts_ns10.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ts_ns16.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ts_ns20.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ts_ns40.prm
echo 'Static bank-sweptv'
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_bs_ns6.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_bs_ns10.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_bs_ns16.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_bs_ns20.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_bs_ns40.prm
echo 'Static IHS-sampled'
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ihs_ns6.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ihs_ns10.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ihs_ns16.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ihs_ns20.prm
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ihs_ns40.prm

