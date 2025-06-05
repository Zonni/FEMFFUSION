#!/bin/bash
set -e  # Stop at first error

# Generate integer snapshot values
n_points=()
for i in {1..1}; do
    n_points+=($(echo "(10773 * 0.02 * $i) / 1" | bc)) 
done

# Base command
executable="./femffusion.exe"
param_file="3D_Langenbuch_static/3D_Langenbuch_LUPODext_group_wise.prm"

# Loop through the snapshot values
for p in "${n_points[@]}"; do
    out_file="3D_Langenbuch_static/LUPOD_points/3D_Langenbuch_LUPODext${p}_group_wise.out"
    echo 'N LUPODext POINTS ' $p 'at' $out_file
    $executable -f $param_file -n_snap 20 -out_file "$out_file" -n_lupod_points "$p"
done