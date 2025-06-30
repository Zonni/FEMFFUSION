#!/bin/bash
# Run this script from FEMFFUSION_DIR
set -e  # Stop at first error

# Create output directory
mkdir -p 3D_Langenbuch_static/mreq/

# Path to parameter file
param_file="3D_Langenbuch_static/3D_Langenbuch_POD_group_wise.prm"

# Define epsilon values
mreq=()
for ((i=2; i<=100; i+=2)); do
  mreq+=($i)
done


export param_file  # For use in parallel subshells

# Function to run one simulation
run_femffusion() {
    m="$1"

    out_file="3D_Langenbuch_static/mreq/3D_Langenbuch_mreq${m}.out"
    log_file="3D_Langenbuch_static/mreq/3D_Langenbuch_mreq_${m}.log"
    echo "   N Snaps Retained $m at $out_file"
    ./femffusion.exe -f "$param_file" -n_snap 100 -m_req $m -out_file "$out_file" 2>&1 | tee "$log_file"
}

export -f run_femffusion

# Run all simulations in parallel using up to 20 jobs
parallel -j32 run_femffusion ::: "${mreq[@]}"
