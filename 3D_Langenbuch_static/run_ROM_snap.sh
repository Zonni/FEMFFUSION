#!/bin/bash
#Run from FEMFFUSION_DIR
set -e  # Stop at first error

# List of snapshot values
#snap_values=(2 3 5 10 15 20 25 35 40 45 50 60 70 80 90 100 125 150 175 200 300 400 500)
snap_values=(40)

# Base command
executable="./femffusion.exe"
param_file="3D_Langenbuch_static/3D_Langenbuch_POD_group_wise.prm"

# Loop through the snapshot values
for snap in "${snap_values[@]}"; do
    echo 'N_snapshots: ' $snap
    out_file="3D_Langenbuch_static/POD_n_snapshots/3D_Langenbuch_POD${snap}_group_wise.out"
    $executable -f $param_file -n_snap $snap -out_file "$out_file"
done