#!/bin/bash

move_log_files() {
   # SLURM workload manager is incapable of creating directories or using 
   # variables when generating log/error files. As a result, scripts create
   # temporary log/error files and this function moves them instead.
   original_name=$1

   mkdir -p "${LOG_DIR}/"
   mv "${SLURM_SUBMIT_DIR}/${original_name}${SLURM_JOB_ID}.log" \
      "${LOG_DIR}/${original_name}${SLURM_JOB_ID}.log"
   mv "${SLURM_SUBMIT_DIR}/${original_name}${SLURM_JOB_ID}.err" \
      "${LOG_DIR}/${original_name}${SLURM_JOB_ID}.err"
}
