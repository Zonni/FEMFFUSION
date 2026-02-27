#!/bin/bash
# Run from FEMFFUSION_DIR
set -e  # Stop at first error

# Setup directory
mkdir -p '3D_VVER440_static_rom/POD_n_snapshots_FE1'

# List of snapshot values
n_snaps=(2 3 4 5 10 15 20 25 30 40  50 60 70 80 90 100)
param_file="3D_VVER440_static_rom/3D_VVER440_POD_FE1.prm"

echo "Launching ${#n_snaps[@]} simulations ..."

# Loop through the array one by one
for ns in "${n_snaps[@]}"; do
    out_file="3D_VVER440_static_rom/POD_n_snapshots_FE1/3D_VVER440_POD${ns}_group_wise.out"
    log_file="./3D_VVER440_static_rom/POD_n_snapshots_FE1/3D_VVER440_POD${ns}_group_wise.log"

    echo "Running POD_group_wise with $ns snapshots..."
    
    # Run the executable
    ./femffusion.exe -f "$param_file" -n_snap "$ns" -out_file "$out_file" > "$log_file" 2>&1
done

echo "All simulations completed."