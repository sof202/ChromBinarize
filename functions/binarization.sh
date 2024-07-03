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
   output_directory=$1

   for chromosome in {1..22} X; do
     awk \
       -v chromosome="$chromosome" \
       '$1 == "chr"chromosome' \
       "${output_directory}/purified_reads.bed" > \
       "${output_directory}/split/purified_chr${chromosome}.bed"
   done
}

binarization_createBlankBins() {
   output_directory=$1

   module purge
   module load R/4.2.1-foss-2022a

   Rscript "${RSCRIPT_DIR}/create_blank_bed_files.R" \
     "$chromosome_sizes" \
     "$bin_size" \
     "${output_directory}/blanks"

   module purge
}

binarization_countSignalIntersectionWithBins() {
   output_directory=$1   

   module purge
   module load BEDTools

   for chromosome in {1..22} X; do
     bedtools intersect \
       -wa \
       -c \
       -a "${output_directory}/blanks/chromosome${chromosome}.bed" \
       -b "${output_directory}/split/purified_chr${chromosome}.bed" > \
       "${output_directory}/bin_counts/chromosome${chromosome}.bed"
   done

   module purge
}

binarization_createChromhmmBinaryFiles() {
   input_directory=$1
   output_directory=$2   
   mark_name=$3

   module purge
   module load R/4.2.1-foss-2022a

   for chromosome in {1..22} X; do
     dense_file="${output_directory}/dense/${cell_type}_chr${chromosome}_binary.txt"
     echo -e "${cell_type}\tchr${chromosome}" > "${dense_file}" 
     echo "${mark_name}_dense" >> "${dense_file}"

     sparse_file="${output_directory}/sparse/${cell_type}_chr${chromosome}_binary.txt"
     echo -e "${cell_type}\tchr${chromosome}" > "${sparse_file}" 
     echo "${mark_name}_sparse" >> "${sparse_file}"

     Rscript "$RSCRIPT_DIR/binarize.R" \
       "${input_directory}/bin_counts/chromosome${chromosome}.bed" \
       "${dense_file}" \
       "${sparse_file}"

     gunzip "${dense_file}" "${sparse_file}"
   done

   module purge
}
