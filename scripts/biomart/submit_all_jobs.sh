#!/bin/bash

# Submit all Biomart coordinate fetching jobs
# Usage: ./submit_all_jobs.sh [all|current|legacy]

# Default is to submit only current genome builds
MODE="${1:-current}"

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "=========================================="
echo "Submitting Biomart Coordinate Fetch Jobs"
echo "Mode: $MODE"
echo "Script directory: $SCRIPT_DIR"
echo "=========================================="
echo ""

# Function to submit a job and print info
submit_job() {
    local script=$1
    local description=$2
    local full_path="${SCRIPT_DIR}/$script"

    if [ -f "$full_path" ]; then
        echo "Submitting: $description"
        JOB_ID=$(sbatch "$full_path" | awk '{print $4}')
        echo "  Job ID: $JOB_ID"
        echo "  Script: $script"
        echo ""
        return 0
    else
        echo "WARNING: Script not found: $full_path"
        return 1
    fi
}

# Track job submissions
SUBMITTED=0
FAILED=0

# Submit current genome builds
if [ "$MODE" == "current" ] || [ "$MODE" == "all" ]; then
    echo "Current genome builds:"
    echo "----------------------"

    if submit_job "fetch_mouse_mm39.slurm" "Mouse mm39 (GRCm39) coordinates"; then
        ((SUBMITTED++))
    else
        ((FAILED++))
    fi

    if submit_job "fetch_human_hg38.slurm" "Human hg38 (GRCh38) coordinates"; then
        ((SUBMITTED++))
    else
        ((FAILED++))
    fi
fi

# Submit legacy genome builds
if [ "$MODE" == "legacy" ] || [ "$MODE" == "all" ]; then
    echo "Legacy genome builds:"
    echo "---------------------"

    if submit_job "fetch_mouse_mm10.slurm" "Mouse mm10 (GRCm38) coordinates"; then
        ((SUBMITTED++))
    else
        ((FAILED++))
    fi

    if submit_job "fetch_human_hg19.slurm" "Human hg19 (GRCh37) coordinates"; then
        ((SUBMITTED++))
    else
        ((FAILED++))
    fi
fi

# Print summary
echo "=========================================="
echo "Summary:"
echo "  Jobs submitted: $SUBMITTED"
if [ $FAILED -gt 0 ]; then
    echo "  Jobs failed to submit: $FAILED"
fi
echo ""
echo "Monitor job status with: squeue -u $USER"
echo "Check output logs in: scripts/biomart/logs/fetch_*_coords_*.out"
echo "Check error logs in: scripts/biomart/logs/fetch_*_coords_*.err"
echo "=========================================="

# Return non-zero if any submissions failed
exit $FAILED