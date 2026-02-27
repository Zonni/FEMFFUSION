#!/bin/bash
#Run from FEMFFUSION_DIR
set -e  # Stop at first error

mkdir -p '3D_Langenbuch_static/POD_n_snapshots'
# List of snapshot values
#n_snaps=(2 3 4 5 10 15 20 25 35 40 45 50 60 70 80 90 100 125 150 175 200 225 250 275 300 325 350 375 400 425 450 475 500)
n_snaps=(250 275 325 350 375 425 450 475)
# Función para ejecutar una simulación (definida para usar en parallel)
run_simulation() {
    local ns=$1
    local out_file="3D_Langenbuch_static/POD_n_snapshots/3D_Langenbuch_POD${ns}_group_wise.out"
    local log_file="./3D_Langenbuch_static/POD_n_snapshots/3D_Langenbuch_POD${ns}_group_wise.log"
    local param_file="3D_Langenbuch_static/3D_Langenbuch_POD_group_wise.prm"

    echo "Running POD_group_wise with $ns snapshots..."
    ./femffusion.exe -f "$param_file" -n_snap $ns -out_file "$out_file" > "$log_file" 2>&1
}
export -f run_simulation

echo "n_snaps: ${n_snaps[@]}"
echo "Launching ${#n_snaps[@]} simulations in parallel..."
# Usar parallel con los valores de n_snaps
printf "%s\n" "${n_snaps[@]}" | parallel -j 2 run_simulation
