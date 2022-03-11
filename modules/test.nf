#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process STARALIGN {

  tag "STAR align on $sample_id"
  publishDir params.outdir, mode: 'copy'

  input:
    tuple val(sample_id), file(reads)
    path index

  output:
    file 'sample_id*'

  script:
    """
    STAR  --genomeDir $index \
          --runMode alignReads \
          --readFilesIn ${reads[0]} ${reads[1]} \
          --outFileNamePrefix $sample_id \
          --runThreadN 16 \
          --outSAMtype BAM SortedByCoordinate \
          --outFilterMultimapNmax 100 \
          --winAnchorMultimapNmax 100
    """
}

