#!/usr/bin/env nextflow

params.index = "/data/scratch/yaochung41/genomeIndex/hg19"
params.reads = "/data/scratch/yaochung41/kap1/data/hm/*_R{1,2}.fastq"
params.outdir = "$baseDir/results"

log.info """\
         STAR - Aligning Reads for TE analysis
         =====================================
         index          : ${params.index}
         reads          : ${params.reads}
         outdir         : ${params.outdir}
         """
         .stripIndent()

Channel
    .fromFilePairs( params.reads )
    .set { reads_ch }

Channel
    .fromPath( params.index )
    .set { index_ch }


process align {
    
    publishDir params.outdir, mode: 'copy'

    input:
        tuple val(sample_id), file(reads) from reads_ch
        path index from index_ch

    output:
        '*.bam'
    
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
