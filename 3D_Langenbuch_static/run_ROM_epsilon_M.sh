#!/bin/bash
# Run this script from FEMFFUSION_DIR
set -e  # Stop at first error

# Create output directory
mkdir -p 3D_Langenbuch_static/epsilonM/

# Path to parameter file
param_file="3D_Langenbuch_static/3D_Langenbuch_POD_group_wise.prm"

# Define epsilon values
eps=(
  3e-1 1e-1
  3e-2 1e-2
  3e-3 1e-3
  3e-4 1e-4
  3e-5 1e-5
  3e-6 1e-6
  3e-7 1e-7
  3e-8 1e-8
  3e-9 1e-9
  3e-10 1e-10
  3e-11 1e-11
  3e-12 1e-12
  3e-13 1e-13
  3e-14 1e-14
  3e-15 1e-15
  3e-16 1e-16
)
#eps=(0.9)

export param_file  # For use in parallel subshells

# Function to run one simulation
run_femffusion() {
    e="$1"
    mantissa="${e%%e-*}"
    exponent="${e##*e-}"
    exponent_padded=$(printf "%02d" "$exponent")
    e_label="${mantissa}em${exponent_padded}"

    out_file="3D_Langenbuch_static/epsilonM/3D_Langenbuch_eps_${e_label}.out"
    log_file="3D_Langenbuch_static/epsilonM/3D_Langenbuch_eps_${e_label}.log"

    echo "   Epsilon_M $e → $out_file  (log: $log_file)"
    ./femffusion.exe -f "$param_file" -n_snap 100 -epsilon_M "$e" -out_file "$out_file" 2>&1 | tee "$log_file"
}

export -f run_femffusion

# Run all simulations in parallel using up to 20 jobs
parallel -j32 run_femffusion ::: "${eps[@]}"
