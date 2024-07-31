#!/bin/bash

move_log_files() {
   # SLURM workload manager is incapable of creating directories or using 
   # variables when generating log/error files. As a result, scripts create
   # temporary log/error files and this function moves them instead.
   original_name=$1
   datetime=$(date +%d-%h~%H-%M)

   mkdir -p "${LOG_DIR}/${USER}/${original_name}"
   mv "${SLURM_SUBMIT_DIR}/${original_name}${SLURM_JOB_ID}.log" \
      "${LOG_DIR}/${USER}/${original_name}/${datetime}.log"
   mv "${SLURM_SUBMIT_DIR}/${original_name}${SLURM_JOB_ID}.err" \
      "${LOG_DIR}/${USER}/${original_name}/${datetime}.err"
}
