#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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


process align {
    
    publishDir params.outdir, mode: 'copy'

    input:
        path reads
        path index

    output:
        tuple val(sample_id), path("${sample_id}")
    
    script:
    """
    STAR  --genomeDir $index \
          --runMode alignReads \
          --readFilesIn ${reads[0]} ${reads[1]} \
          --outFileNamePrefix ${reads.baseName} \
          --runThreadN 16 \
          --outSAMtype BAM SortedByCoordinate \
          --outFilterMultimapNmax 100 \
          --winAnchorMultimapNmax 100
    """
}

Channel
    .fromPath( params.reads )
    .set { reads_ch }

Channel
    .fromPath( params.index )
    .set { index_ch }

workflow {
    align(reads_ch, index_ch)
}
