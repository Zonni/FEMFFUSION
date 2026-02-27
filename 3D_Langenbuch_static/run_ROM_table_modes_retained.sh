#!/bin/bash
#Run from FEMFFUSION_DIR
set -e  # Stop at first error

# List of snapshot value
param_file="3D_Langenbuch_static/3D_Langenbuch_POD_group_wise.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 50 -m_req 10 -out_file '3D_Langenbuch_static/3D_Langenbuch_POD50_m10_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 50 -m_req 20 -out_file '3D_Langenbuch_static/3D_Langenbuch_POD50_m20_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 50 -m_req 40 -out_file '3D_Langenbuch_static/3D_Langenbuch_POD50_m40_group_wise.out'


param_file="3D_Langenbuch_static/3D_Langenbuch_LUPOD_group_wise.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 50 -m_req 10 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPOD50_m10_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 50 -m_req 20 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPOD50_m20_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 50 -m_req 40 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPOD50_m40_group_wise.out'

param_file="3D_Langenbuch_static/3D_Langenbuch_LUPODext_group_wise.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 50 -m_req 10 -n_lupod_points 4309 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPODext50_m10_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 50 -m_req 20 -n_lupod_points 4309 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPODext50_m20_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 50 -m_req 40 -n_lupod_points 4309 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPODext50_m40_group_wise.out'


param_file="3D_Langenbuch_static/3D_Langenbuch_LUPODext_group_wise.prm"
mpirun -n 1 femffusion.exe -f $param_file -n_snap 50 -m_req 10 -n_lupod_points 10 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPOD50_m10_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 50 -m_req 20 -n_lupod_points 20 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPOD50_m20_group_wise.out'
mpirun -n 1 femffusion.exe -f $param_file -n_snap 50 -m_req 40 -n_lupod_points 40 -out_file '3D_Langenbuch_static/3D_Langenbuch_LUPOD50_m40_group_wise.out'

