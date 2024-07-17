#!/bin/bash

purification_convertBSBedToMethylBedFormat() {
  # The default format required for BS-Seq files is not the same as ONT. This
  # is because ONT bed files are usually created with modkit and our lab uses
  # wgbs_tools for BS-Seq bed file creation. These do not output the same 
  # format of bed file, as such this function is here to normalise the format.
  mark_name=$1
  input_bed_file=$2
  output_file_name=$3

logs "${DEBUG_MODE:0}" \
"Converting ${input_bed_file} into methylbed format."

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
  mark_name=$1
  input_bed_file=$2
  output_file_name=$3

  if [[ "${mark_name}" =~ "m" ]]; then
    reference_percent_threshold="${reference_percentage_threshold_m:=5}"
    reference_read_depth_threshold="${reference_read_depth_threshold_m:=500}"
  elif [[ "${mark_name}" == "h" ]]; then
    reference_percent_threshold="${reference_percentage_threshold_h:=5}"
    reference_read_depth_threshold="${reference_read_depth_threshold_h:=50}"
  fi

logs "${DEBUG_MODE:0}" \
"Creating 'good reference set' from ${input_bed_file} where:
Percent methylation for sites are at most ${reference_percent_threshold},
Number of reads for sites are at least ${reference_read_depth_threshold}."

  awk -v percent_threshold="${reference_percent_threshold}" \
    -v read_threshold="${reference_read_depth_threshold}" \
    -v mark_name="$mark_name" \
    '{OFS="\t"}
    $4 == mark_name && $5 >= read_threshold && $7 <= percent_threshold {print $5,$7}' \
    "${input_bed_file}" > \
      "${output_file_name}"

  if [[ ! -s "${output_file_name}" ]]; then
errors "${output_file_name} is empty.
Your thresholds in the config file are likely too high for your dataset."
  fi
}

purification_filterOutLowReadDepthSites() {
  mark_name=$1
  input_bed_file=$2
  output_file_name=$3

logs "${DEBUG_MODE:0}" \
"Filtering sites from ${input_bed_file} where \
the number of reads is at least ${reference_read_depth_threshold}."

  awk -v read_threshold="${minimum_read_depth}" \
    -v mark_name="$mark_name" \
    '{OFS="\t"} 
    $4 == mark_name && $5 >= read_threshold' \
    "${input_bed_file}" > \
      "${output_file_name}"

  if [[ ! -s "${output_file_name}" ]]; then
errors "${output_file_name} is empty.
Your read threshold in the config file is likely too high for your dataset."
  fi
}

purification_calculateSiteMethylationProbability() {
  processing_directory=$1
  reference_set=$2
  input_file=$3
  output_file=$4

logs "${DEBUG_MODE:0}" \
"Calculating the probability that methylated reads are due to sequencing or \
basecalling errors."

  module purge
  module load R/4.2.1-foss-2022a

  Rscript "${RSCRIPT_DIR}/binom.R" \
    "${processing_directory}" \
    "${reference_set}" \
    "${input_file}" \
    "${output_file}"

  module purge
}

purification_removeDeterminedUnmethylatedSites() {
  input_file=$1
  output_file=$1

logs "${DEBUG_MODE:0}" \
"Removing sites that are deemed unmethylated (non-significantly methylated)"

  awk -v threshold="${binomial_threshold}" \
    '$8 < threshold' \
    "${input_file}" > \
    "${output_file}"

  if [[ ! -s "${output_file}" ]]; then
errors "${output_file} is empty.
Your binomial threshold in the config file is likely too low.
Alternatively, your read/percent thresholds may be too lenient."
  fi
}
