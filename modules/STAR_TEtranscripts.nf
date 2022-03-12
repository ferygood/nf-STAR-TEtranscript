#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process STARALIGN {

  tag "STAR align on $sample_id"
  publishDir "$params.outdir", mode: 'copy'

  input:
    tuple val(sample_id), file(reads)
    path index

  output:
    file '*.bam'

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


process TEcount {

  tag "TEcount quantified on $bam_id"
  publishDir "$params.quantdir", mode: 'copy'

  input:
    tuple val(bam_id), file(bam_file)
    file gtf
    file rmsk_ind

  output:
    file '*.cntTable'

  script:
    """
    TEcount -b $bam_file \
            --GTF $gtf \
            --TE $rmsk_ind \
            --sortByPos \
            --project $bam_id
    """
    
}
