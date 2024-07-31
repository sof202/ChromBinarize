#!/bin/bash

binarization_createDirectories() {
  output_directory=$1   

  rm -rf "${output_directory}/split" \
    "${output_directory}/blanks" \
    "${output_directory}/bin_counts" \
    "${output_directory}/binarized"

  mkdir -p "${output_directory}/split" \
    "${output_directory}/blanks" \
    "${output_directory}/bin_counts" \
    "${output_directory}/binarized/dense" \
    "${output_directory}/binarized/sparse"
}

binarization_splitIntoChromosomes() {
  # ChromHMM requires a binary file for each chromosome, whilst bed files
  # usually cover the whole genome. Hence, splitting of the bed files is
  # required
  processing_directory=$1
  input_file_name=$2

logs "${DEBUG_MODE:0}" \
"Splitting ${input_file_name} into chromosomes 1-22 and X..."

  for chromosome in {1..22} X; do
    awk \
      -v chromosome="$chromosome" \
      '$1 == "chr"chromosome' \
      "${processing_directory}/${input_file_name}" > \
      "${processing_directory}/split/purified_chr${chromosome}.bed"
    done
}

binarization_createBlankBins() {
  output_directory=$1
  bin_size=$2

logs "${DEBUG_MODE:0}" \
"Creating blank bed files for chromosomes 1-22 and X..."

  if [[ ! -f "${chromosome_sizes}" ]]; then
errors "Could not find your chromosome sizes file at:
${chromosomes_sizes}.
Make sure the variable in the config file is pointing to the correct location."
  fi

  if [[ $(awk 'NR == 1 {print NF}' "${chromosome_sizes}") -ne 2 ]]; then
errors "Your chromosome sizes file at:
${chromosome_sizes}
does not have 2 columns. This is likely to cause the script to fail.
This file should be of the form:
chromosome    size"
  fi
  

  Rscript "${RSCRIPT_DIR}/create_blank_bed_files.R" \
    "${REPO_DIR}" \
    "${chromosome_sizes}" \
    "${bin_size}" \
    "${output_directory}/blanks"
}

binarization_countSignalIntersectionWithBins() {
  output_directory=$1   

logs "${DEBUG_MODE:0}" \
"Calculating the number of methylated sites in each ${bin_size}bp bin..."

  for chromosome in {1..22} X; do
    bedtools intersect \
      -wa \
      -c \
      -a "${output_directory}/blanks/chromosome${chromosome}.bed" \
      -b "${output_directory}/split/purified_chr${chromosome}.bed" > \
      "${output_directory}/bin_counts/chromosome${chromosome}.bed"
    done
}

binarization_createChromhmmBinaryFiles() {
  input_directory=$1
  output_directory=$2   
  mark_name=$3

logs "${DEBUG_MODE:0}" \
"Generating ChromHMM binary files...
Sparse files will be placed in: ${output_directory}-sparse
Dense files will be placed in: ${output_directory}-dense"

  rm -rf "${output_directory}-dense" "${output_directory}-sparse"
  mkdir -p "${output_directory}-dense" "${output_directory}-sparse"

  for chromosome in {1..22} X; do
    dense_file="${output_directory}-dense/${cell_type}_chr${chromosome}_binary.txt"
    echo -e "${cell_type}\tchr${chromosome}" > "${dense_file}" 
    echo "${mark_name}_dense" >> "${dense_file}"

    sparse_file="${output_directory}-sparse/${cell_type}_chr${chromosome}_binary.txt"
    echo -e "${cell_type}\tchr${chromosome}" > "${sparse_file}" 
    echo "${mark_name}_sparse" >> "${sparse_file}"

    Rscript "$RSCRIPT_DIR/binarize.R" \
      "${REPO_DIR}" \
      "${input_directory}/bin_counts/chromosome${chromosome}.bed" \
      "${dense_file}" \
      "${sparse_file}" \
      "${bin_size}" \
      "${beta_threshold}"

    number_of_dense_signatures=$(awk 'NR>2 && $1>0' "${dense_file}" | wc -l)
    number_of_sparse_signatures=$(awk 'NR>2 && $1>0' "${sparse_file}" | wc -l)

logs 1 \
"chromosome ${chromosome} has:
${number_of_dense_signatures} dense signatures,
${number_of_sparse_signatures} sparse signatures."

    if [[ "${number_of_dense_signatures}" -eq 0 ]]; then
errors "${chromosome}'s binary file has no true/1 entries. 
Either this chromosome's signatures are too sparse, \
or your 'binomial threshold' in the config file is too strict."
    fi

    gzip "${dense_file}" "${sparse_file}"
  done
}
