#!/bin/bash
#Run from FEMFFUSION_DIR
set -e  # Stop at first error

# List of snapshot value
param_file="3D_Langenbuch_static/3D_Langenbuch_POD_group_wise.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 10 -out_file '3D_Langenbuch_static/3D_Langenbuch_POD10_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 20 -out_file '3D_Langenbuch_static/3D_Langenbuch_POD20_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 40 -out_file '3D_Langenbuch_static/3D_Langenbuch_POD40_group_wise.out'


param_file="3D_Langenbuch_static/3D_Langenbuch_LUPOD_group_wise.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 10 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPOD10_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 20 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPOD20_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 40 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPOD40_group_wise.out'

param_file="3D_Langenbuch_static/3D_Langenbuch_LUPODext_group_wise.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 10  -n_lupod_points 4309 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPODext10_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 20  -n_lupod_points 4309 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPODext20_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 40  -n_lupod_points 4309 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPODext40_group_wise.out'


param_file="3D_Langenbuch_static/3D_Langenbuch_Random_group_wise.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 10  -n_lupod_points 4309 -out_file '3D_Langenbuch_static/3D_Langenbuch_Random10_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 20  -n_lupod_points 4309 -out_file '3D_Langenbuch_static/3D_Langenbuch_Random20_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 40  -n_lupod_points 4309 -out_file '3D_Langenbuch_static/3D_Langenbuch_Random40_group_wise.out'


param_file="3D_Langenbuch_static/3D_Langenbuch_FEM1_group_wise.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 10  -out_file '3D_Langenbuch_static/3D_Langenbuch_FEM1_10_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 20  -out_file '3D_Langenbuch_static/3D_Langenbuch_FEM1_20_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 40  -out_file '3D_Langenbuch_static/3D_Langenbuch_FEM1_40_group_wise.out'


param_file="3D_Langenbuch_static/3D_Langenbuch_POD_mono.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 10 -out_file '3D_Langenbuch_static/3D_Langenbuch_POD10_mono.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 20 -out_file '3D_Langenbuch_static/3D_Langenbuch_POD20_mono.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 40 -out_file '3D_Langenbuch_static/3D_Langenbuch_POD40_mono.out'


param_file="3D_Langenbuch_static/3D_Langenbuch_LUPOD_mono.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 10 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPOD10_mono.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 20 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPOD20_mono.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 40 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPOD40_mono.out'

param_file="3D_Langenbuch_static/3D_Langenbuch_LUPODext_mono.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 10 -n_lupod_points 8618 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPODext10_mono.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 20 -n_lupod_points 8618 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPODext20_mono.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 40 -n_lupod_points 8618 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPODext40_mono.out'



