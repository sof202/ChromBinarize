#!/bin/bash

purification_convertBSBedToMethylBedFormat() {
  # The default format required for BS-Seq files is not the same as ONT. This
  # is because ONT bed files are usually created with modkit and our lab uses
  # wgbs_tools for BS-Seq bed file creation. These do not output the same 
  # format of bed file, as such this function is here to normalise the format.
  output_file_name=$1
  input_bed_file=$2
  mark_name=$3

  awk -v mark_name="$mark_name" \
    'function convert_to_percent(reads, total_reads) {
      return (int(reads / total_reads * 10000) / 100)
    }
    {OFS="\t"}
    {print $1,$2,$3,mark_name,$5,"+",convert_to_percent($4,$5)}' \
    "${input_bed_file}" > \
      "${output_file_name}"
}

purification_extractSitesWithLowMethylation() {
  output_directory=$1
  input_bed_file=$2
  mark_name=$3

  if [[ "${mark_name}" == "m" ]]; then
    reference_percent_threshold="${reference_percentage_threshold_m:=95}"
    reference_read_depth_threshold="${reference_read_depth_threshold_m:=500}"
  elif [[ "${mark_name}" == "h" ]]; then
    reference_percent_threshold="${reference_percentage_threshold_h:=95}"
    reference_read_depth_threshold="${reference_read_depth_threshold_h:=50}"
  fi

  awk -v percent_threshold=$((100 - reference_percent_threshold)) \
    -v read_threshold="${reference_read_depth_threshold}" \
    -v mark_name="$mark_name" \
    '{OFS="\t"}
    $4 == mark_name && $5 >= read_threshold && $7 <= percent_threshold {print $5","$7}' \
    "${input_bed_file}" > \
      "${output_directory}/unmethylated_reads.bed"
}

purification_filterOutLowReadDepthSites() {
  output_directory=$1
  input_bed_file=$2
  mark_name=$3

  awk -v read_threshold="${minimum_read_depth}" \
    -v mark_name="$mark_name" \
    '{OFS="\t"} 
    $4 == mark_name && $5 >= read_threshold' \
    "${input_bed_file}" > \
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
    '$8 < threshold' \
    "${output_directory}/processed_reads.bed" > \
      "${output_directory}/purified_reads.bed"
}
