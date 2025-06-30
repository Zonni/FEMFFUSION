#!/bin/bash
#Run it from FEMFFUSION_DIR!
set -e  # Stop at first error

## List of snapshot value
param_file="3D_VVER440_static_rom/3D_VVER440_POD_FE1.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 10 -out_file '3D_VVER440_static_rom/3D_VVER440_POD10_FE1.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 20 -out_file '3D_VVER440_static_rom/3D_VVER440_POD20_FE1.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 40 -out_file '3D_VVER440_static_rom/3D_VVER440_POD40_FE1.out'

param_file="3D_VVER440_static_rom/3D_VVER440_LUPODext_FE1.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 10 -out_file '3D_VVER440_static_rom/3D_VVER440_LUPODext10_FE1.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 20 -out_file '3D_VVER440_static_rom/3D_VVER440_LUPODext20_FE1.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 40 -out_file '3D_VVER440_static_rom/3D_VVER440_LUPODext40_FE1.out'

param_file="3D_VVER440_static_rom/3D_VVER440_LUPODext_FE2.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 10 -out_file '3D_VVER440_static_rom/3D_VVER440_LUPODext10_FE2.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 20 -out_file '3D_VVER440_static_rom/3D_VVER440_LUPODext20_FE2.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 40 -out_file '3D_VVER440_static_rom/3D_VVER440_LUPODext40_FE2.out'

param_file="3D_VVER440_static_rom/3D_VVER440_POD_FE2.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 10 -out_file '3D_VVER440_static_rom/3D_VVER440_POD10_FE2.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 20 -out_file '3D_VVER440_static_rom/3D_VVER440_POD20_FE2.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 40 -out_file '3D_VVER440_static_rom/3D_VVER440_POD40_FE2.out'

param_file="3D_VVER440_static_rom/3D_VVER440_POD_FE3.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_test 10 -n_snap 10 -out_file '3D_VVER440_static_rom/3D_VVER440_POD10_FE3.out'
mpirun -n 1 femffusion.exe -f $param_file -n_test 10 -n_snap 20 -out_file '3D_VVER440_static_rom/3D_VVER440_POD20_FE3.out'
mpirun -n 1 femffusion.exe -f $param_file -n_test 10 -n_snap 40 -out_file '3D_VVER440_static_rom/3D_VVER440_POD40_FE3.out'

param_file="3D_VVER440_static_rom/3D_VVER440_LUPODext_FE3.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_test 10 -n_snap 10 -out_file '3D_VVER440_static_rom/3D_VVER440_LUPODext10_FE3.out'
mpirun -n 1 femffusion.exe -f $param_file -n_test 10 -n_snap 20 -out_file '3D_VVER440_static_rom/3D_VVER440_LUPODext20_FE3.out'
mpirun -n 1 femffusion.exe -f $param_file -n_test 10 -n_snap 40 -out_file '3D_VVER440_static_rom/3D_VVER440_LUPODext40_FE3.out'
