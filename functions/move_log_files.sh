#!/bin/bash

move_log_files() {
    # Moves log files from job submission directory to designated output
    # directory.
    original_name=$1

    mkdir -p "${LOG_DIR}/"
    mv "${SLURM_SUBMIT_DIR}/${original_name}${SLURM_JOB_ID}.log" \
       "${LOG_DIR}/${original_name}${SLURM_JOB_ID}.log"
    mv "${SLURM_SUBMIT_DIR}/${original_name}${SLURM_JOB_ID}.err" \
       "${LOG_DIR}/${original_name}${SLURM_JOB_ID}.err"
}
