#!/bin/bash
set -e  # Stop at first error
mkdir -p '3D_Langenbuch_static/LUPOD_points'

# Generate integer snapshot values
n_points=()
for i in {1..50}; do
    n_points+=($(echo "(10773 * 0.02 * $i) / 1" | bc))
done

# Función para ejecutar una simulación (definida para usar en parallel)
run_simulation() {
    local p=$1
    local out_file="./3D_Langenbuch_static/LUPOD_points/3D_Langenbuch_LUPODext${p}_group_wise.out"
    local log_file="./3D_Langenbuch_static/LUPOD_points/3D_Langenbuch_LUPODext${p}_group_wise.log"

    # Base command
    executable="./femffusion.exe"
    param_file="3D_Langenbuch_static/3D_Langenbuch_LUPODext_group_wise.prm"

    echo "Running LUPODext with $p points..."
    "$executable" -f "$param_file" -n_snap 20 -out_file "$out_file" -n_lupod_points "$p" > "$log_file" 2>&1
}
export -f run_simulation

echo "n_points: ${n_points[@]}"
echo "Launching ${#n_points[@]} simulations in parallel..."
# Usar parallel con los valores de n_points
printf "%s\n" "${n_points[@]}" | parallel -j 38 run_simulation
