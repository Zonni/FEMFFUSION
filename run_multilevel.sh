#!/bin/bash

./femffusion.exe -f 3D_VVER440_NOISE/VVER440_EFPD_007_fe2.prm -pc_complex gauss_seidel_2x2 -ksp_monitor
./femffusion.exe -f 3D_VVER440_NOISE/VVER440_EFPD_007_fe2.prm -pc_complex gauss_seidel -ksp_monitor
./femffusion.exe -f 3D_VVER440_NOISE/VVER440_EFPD_007_fe2.prm -pc_complex multilevel -ksp_monitor -n_fe_coarse 1 -pc_coarse gs-cgilu -matrixfree_type_A non_diagonal
./femffusion.exe -f 3D_VVER440_NOISE/VVER440_EFPD_007_fe2.prm -pc_complex multilevel -ksp_monitor -n_fe_coarse 1 -pc_coarse gs-cgilu -matrixfree_type_A full_matrixfree
./femffusion.exe -f 3D_VVER440_NOISE/VVER440_EFPD_007_fe2.prm -pc_complex multilevel -ksp_monitor -n_fe_coarse 1  -matrixfree_type_A non_diagonal
./femffusion.exe -f 3D_VVER440_NOISE/VVER440_EFPD_007_fe2.prm -pc_complex multilevel -ksp_monitor -n_fe_coarse 1  -matrixfree_type_A full_matrixfree

./femffusion.exe -f 3D_VVER440_NOISE/VVER440_EFPD_007_fe3.prm -pc_complex gauss_seidel_2x2 -ksp_monitor
./femffusion.exe -f 3D_VVER440_NOISE/VVER440_EFPD_007_fe3.prm -pc_complex gauss_seidel -ksp_monitor
./femffusion.exe -f 3D_VVER440_NOISE/VVER440_EFPD_007_fe3.prm -pc_complex multilevel -ksp_monitor -n_fe_coarse 1 -pc_coarse gs-cgilu -matrixfree_type_A non_diagonal
./femffusion.exe -f 3D_VVER440_NOISE/VVER440_EFPD_007_fe3.prm -pc_complex multilevel -ksp_monitor -n_fe_coarse 1 -pc_coarse gs-cgilu -matrixfree_type_A full_matrixfree
./femffusion.exe -f 3D_VVER440_NOISE/VVER440_EFPD_007_fe3.prm -pc_complex multilevel -ksp_monitor -n_fe_coarse 1  -matrixfree_type_A non_diagonal
./femffusion.exe -f 3D_VVER440_NOISE/VVER440_EFPD_007_fe3.prm -pc_complex multilevel -ksp_monitor -n_fe_coarse 1  -matrixfree_type_A full_matrixfree