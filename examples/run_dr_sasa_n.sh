#!/bin/bash

# Configuration
TEST_PDB="/home/alessio/dr_sasa_python/tests/data/pred.pdb"
DR_SASA="/home/alessio/dr_sasa_python/dr_sasa_n/build/dr_sasa"

OUTPUT_DIR="test_results"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Create results directory
RESULTS_DIR="${OUTPUT_DIR}/${TIMESTAMP}"
mkdir -p "$RESULTS_DIR"

# Log file
LOG_FILE="${RESULTS_DIR}/test_log.txt"

# Function to log messages
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Function to run test and save output
run_test() {
    local mode=$1
    local description=$2
    local args=$3
    local output_base="${RESULTS_DIR}/mode${mode}"
    
    log_message "Running test: $description"
    log_message "Command: $DR_SASA $args -i $TEST_PDB"
    
    # Run dr_sasa and capture output
    $DR_SASA $args -i "$TEST_PDB" > "${output_base}.out" 2>&1
    
    # Save specific outputs based on mode
    if [ $? -eq 0 ]; then
        log_message "Test completed successfully"
    else
        log_message "Test failed! Check ${output_base}.out for details"
    fi
}

log_message "Starting DR_SASA tests with structure: $TEST_PDB"

# Mode 0 (Simple SASA)
log_message "Testing Mode 0 - Simple SASA solver"
run_test 0 "Basic SASA" "-m 0"
run_test 0 "SASA probe 1.2" "-m 0 -r 1.2"
run_test 0 "SASA chain A" "-m 0 -chain A"

# Mode 1 (Delta SASA)
log_message "Testing Mode 1 - Delta SASA"
run_test 1 "Auto mode" "-m 1"
run_test 1 "Two chains" "-m 1 -chain AB -chain CD"
run_test 1 "No matrix" "-m 1 -chain AB -nomatrix"

# Mode 2 (Residue dSASA)
log_message "Testing Mode 2 - Residue dSASA"
run_test 2 "Single chain" "-m 2 -chain A"
run_test 2 "Multiple chains" "-m 2 -chain ABCD"

# Mode 3 (Atom dSASA)
log_message "Testing Mode 3 - Atom dSASA"
run_test 3 "Single chain" "-m 3 -chain A"
run_test 3 "Multiple chains" "-m 3 -chain ABCD"

# Mode 4 (Contact surface)
log_message "Testing Mode 4 - Contact surface"
run_test 4 "Auto mode" "-m 4"
run_test 4 "Specific chains" "-m 4 -chain AB"

# Mode 100 (Relative SASA)
log_message "Testing Mode 100 - Relative SASA"
run_test 100 "Basic relative" "-m 100"
run_test 100 "Chain selection" "-m 100 -chain A"

# Mode 106 (ASA by atom type)
log_message "Testing Mode 106 - ASA by atom type"
run_test 106 "Atom type analysis" "-m 106 -chain A"

# Additional variations
log_message "Testing additional flag combinations"
run_test 0 "Molecular surface" "-m 0 -r 0.0"
run_test 0 "Force reorder" "-m 0 -force_reorder"
run_test 0 "No HETATM" "-m 0 -nohetatm"

# Create summary file
cat > "${RESULTS_DIR}/summary.txt" << EOL
DR_SASA Test Results Summary
===========================
Date: $(date)
Structure: ${TEST_PDB}

Tests Run:
1. Mode 0 (Simple SASA):
   - Basic calculation
   - Probe radius 1.2
   - Chain selection

2. Mode 1 (Delta SASA):
   - Auto mode
   - Chain groups
   - No matrix

3. Mode 2 (Residue dSASA):
   - Single/Multiple chains

4. Mode 3 (Atom dSASA):
   - Single/Multiple chains

5. Mode 4 (Contact surface):
   - Auto/Chain modes

6. Mode 100 (Relative SASA):
   - Basic/Chain selection

7. Mode 106 (ASA by atom type)

Additional Tests:
- Zero probe radius
- Force reorder
- No HETATM

Results Location: ${RESULTS_DIR}
Log File: ${LOG_FILE}
EOL

log_message "Testing completed. Results saved in ${RESULTS_DIR}"
log_message "See summary.txt for overview"