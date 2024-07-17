#!/bin/bash

intersect_intersectBSWithONT() {
  ONT_bed_file=$1
  BS_bed_file=$2
  output_file_path=$3

logs "${DEBUG_MODE:0}" \
"Intersecting:
Bisulphite sequencing file: ${BS_bed_file}
with
Oxford Nanopore sequencing file: ${ONT_bed_file}."

  module purge
  module load BEDTools

  bedtools intersect \
    -wo \
    -a "${ONT_bed_file}" \
    -b "${BS_bed_file}" | \
    awk \
    'function convert_to_percent(reads, total_reads) {
      return (int(reads / total_reads * 10000) / 100)
    }
    {OFS="\t"} 
    {print $1,$2,$3,$4,$5,$7,$12,convert_to_percent($11,$12)}' > \
      "${output_file_path}"

    module purge
}
