#!/bin/bash
set -e  # Stop at first error

echo 'FOM'
./femffusion.exe -f 3D_Roseton/roseton_FE3_60p.prm 
echo 'DYNAMIC'
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_dyn_ns6.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_dyn_ns6.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_dyn_ns10.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_dyn_ns10.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_dyn_ns16.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_dyn_ns16.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_dyn_ns20.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_dyn_ns20.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_dyn_ns40.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_dyn_ns40.out
echo 'Static time-sampled'
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ts_ns6.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_ts_ns6.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ts_ns10.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_ts_ns10.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ts_ns16.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_ts_ns16.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ts_ns20.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_ts_ns20.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ts_ns40.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_ts_ns40.out
echo 'Static bank-sweptv'
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_bs_ns6.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_bs_ns6.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_bs_ns10.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_bs_ns10.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_bs_ns16.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_bs_ns16.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_bs_ns20.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_bs_ns20.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_bs_ns40.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_bs_ns40.out
echo 'Static IHS-sampled'
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ihs_ns6.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_ihs_ns6.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ihs_ns10.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_ihs_ns10.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ihs_ns16.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_ihs_ns16.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ihs_ns20.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_ihs_ns20.out
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rom_ihs_ns40.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rom_ihs_ns40.out
echo 'RPOD'
./femffusion.exe -f 3D_Roseton/roseton_FE2_60p_rrom.prm -fe_degree 3 -out_file 3D_Roseton/roseton_FE3_60p_rrom.out

