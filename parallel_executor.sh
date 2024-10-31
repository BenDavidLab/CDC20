#!/bin/bash

# Constants
NUM_WORKERS=4  # Set the number of workers (W)
OUTPUT_DIR="results/PRISM_23_correlations"
R_PROFILE_PATH=""
R_SCRIPT_COMMAND="Rscript --no-save --no-restore drug_sensitivity_gene_exp_corr.R"

mkdir -p "$OUTPUT_DIR"

# Function to run the R script with a given index
run_rscript() {
    local i=$1
    export R_PROFILE="$R_PROFILE_PATH"
    $R_SCRIPT_COMMAND $i $OUTPUT_DIR
}


# Export the function to be available in subshells
export -f run_rscript
export R_SCRIPT_COMMAND OUTPUT_DIR R_PROFILE_PATH

echo "Starting parallel execution with $NUM_WORKERS workers"

seq 1 6525 | xargs -P $NUM_WORKERS -I {} bash -c 'run_rscript "$@"' _ {}

echo "All tasks are done!"
