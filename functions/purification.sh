#!/bin/bash

purification_extractSitesWithHighMethylation() {
  output_directory=$1

  awk -v percent_threshold="${reference_percentage_threshold_h}" \
    -v read_threshold="${reference_read_depth_threshold_h}" \
    'function convert_to_percent(reads, total_reads) {
       return (int(reads / total_reads * 10000) / 100)
     }
     $5 >= read_threshold && convert_to_percent($4,$5) >= percent_threshold {print $5","convert_to_percent($4,$5)}' \
    "${oxBS_bed_file_location}" > \
      "${output_directory}/methylated.csv"
}

purification_extractSitesWithLowMethylation() {
  output_directory=$1

  awk -v percent_threshold=$((100 - ${reference_percentage_threshold_h})) \
    -v read_threshold="${reference_read_depth_threshold_h}" \
    'function convert_to_percent(reads, total_reads) {
       return (int(reads / total_reads * 10000) / 100)
     }
     $5 >= read_threshold && convert_to_percent($4,$5) >= percent_threshold {print $5","convert_to_percent($4,$5)}' \
    "${oxBS_bed_file_location}" > \
      "${output_directory}/unmethylated.csv"
}

purification_filterOutLowReadDepthSites() {
  output_directory=$1

  awk -v read_threshold="${minimum_read_depth}" \
    'function convert_to_percent(reads, total_reads) {
       return (int(reads / total_reads * 10000) / 100)
     }
    {OFS="\t"} 
    $5 >= read_threshold {print $1,$2,$3,"h",$5,"+",convert_to_percent($4,$5)}' \
    "${oxBS_bed_file_location}" > \
      "${output_directory}/filtered_reads.bed"
}

purification_calculateSiteMethylationProbability() {
  output_directory=$1

  module purge
  module load R/4.2.1-foss-2022a

  Rscript "${RSCRIPT_DIR}/binom.R" "${output_directory}"

  module purge
}

purification_removeDeterminedUnmethylatedSites() {
  output_directory=$1

  awk -v threshold="${binomial_threshold}" \
    '$9 < threshold' \
    "${output_directory}/processed_reads.bed" > \
    "${output_directory}/purified_reads.bed"
}
