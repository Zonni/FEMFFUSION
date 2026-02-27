#!/bin/bash
# Run from FEMFFUSION_DIR
set -e  # Stop at first error

# Setup directory
mkdir -p '3D_VVER440_static_rom/POD_n_snapshots_FE2'

# List of snapshot values
n_snaps=(4 5 10 15 20 25 30 40  50 60 70 80 90 100)
param_file="3D_VVER440_static_rom/3D_VVER440_LUPODext_FE2.prm"

echo "Launching ${#n_snaps[@]} simulations ..."

# Loop through the array one by one
for ns in "${n_snaps[@]}"; do
    out_file="3D_VVER440_static_rom/POD_n_snapshots_FE2/3D_VVER440_LUPODext${ns}_group_wise.out"
    log_file="./3D_VVER440_static_rom/POD_n_snapshots_FE2/3D_VVER440_LUPODext${ns}_group_wise.log"

    echo "Running LUPOD EXT group_wise with $ns snapshots..."
    # Run the executable
    ./femffusion.exe -f "$param_file" -n_snap "$ns" -out_file "$out_file" 
done

echo "All simulations completed."