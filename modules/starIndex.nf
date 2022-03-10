#!/usr/bin/env nextflow

log.info """\
         STAR - building Index
         ============================
         GTF          : ${params.gtf}
         FASTA        : ${params.fasta}
         outDirPrefix : ${params.outDirPrefix}
         """
         .stripIndent()

nextflow.enable.dsl=2

// declare paramters
params.gtf = "/data/scratch/yaochung41/gtf_file/hg19.ensGene.gtf"
params.fasta = "/data/scratch/yaochung41/reference/hg19.fa"
params.outDirPrefix = "/data/scratch/yaochung41/genomeIndex/hg19"


process buildIndex {

    input:
        path gtf 
        path fasta 
        path outDirPrefix 
        
    output: 
    
    script:
    """
    STAR --runThreadN 20 \
         --runMode genomeGenerate \
         --genomeDir $outDirPrefix \
         --genomeFastaFiles $fasta \
         --sjdbGTFfile $gtf \
         --sjdbOverhang 100

    """
    
}

workflow {
    gtf_ch = Channel.fromPath( params.gtf, checkIfExists: true)
    fasta_ch = Channel.fromPath( params.fasta, checkIfExists: true)
    outDirPrefix_ch = Channel.fromPath( params.outDirPrefix )
    buildIndex(gtf_ch, fasta_ch, outDirPrefix_ch)
}
