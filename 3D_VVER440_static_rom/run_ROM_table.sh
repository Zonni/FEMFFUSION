#!/bin/bash
#Run it from FEMFFUSION_DIR!
set -e  # Stop at first error

# List of snapshot value


param_file="3D_VVER440_static_rom/3D_VVER440_POD.prm"

#./femffusion.exe -f $param_file -n_snap 10 -out_file  '3D_VVER440_static_rom/3D_VVER440_POD10.out'
#./femffusion.exe -f $param_file -n_snap 20 -out_file  '3D_VVER440_static_rom/3D_VVER440_POD20.out'
#./femffusion.exe -f $param_file -n_snap 40 -out_file  '3D_VVER440_static_rom/3D_VVER440_POD40.out'

#param_file="3D_VVER440_static_rom/3D_VVER440_LUPOD.prm"
#./femffusion.exe -f $param_file -n_snap 10 -out_file  '3D_VVER440_static_rom/3D_VVER440_LUPOD10.out'
#./femffusion.exe -f $param_file -n_snap 20 -out_file  '3D_VVER440_static_rom/3D_VVER440_LUPOD20.out'
#./femffusion.exe -f $param_file -n_snap 40 -out_file  '3D_VVER440_static_rom/3D_VVER440_LUPOD40.out'

param_file="3D_VVER440_static_rom/3D_VVER440_LUPODext.prm"
./femffusion.exe -f $param_file -n_snap 10  -n_lupod_points 52030 -out_file '3D_VVER440_static_rom/3D_VVER440_LUPODext10_group_wise.out' 
./femffusion.exe -f $param_file -n_snap 20  -n_lupod_points 52030 -out_file '3D_VVER440_static_rom/3D_VVER440_LUPODext20_group_wise.out'
./femffusion.exe -f $param_file -n_snap 40  -n_lupod_points 52030 -out_file '3D_VVER440_static_rom/3D_VVER440_LUPODext40_group_wise.out'
