#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Default pipeline parameters
 */
params.reads = '/data/scratch/yaochung41/kap1/data/hm/*_R{1,2}.fastq' 
params.indexDir = '/data/scratch/yaochung41/genomeIndex/hg19'
params.outdir = './results'
params.gtf = "/data/scratch/yaochung41/gtf_file/hg19.ensGene.gtf"
params.rmsk_ind = "/data/scratch/yaochung41/gtf_file/hg19_rmsk_TE.gtf.ind"
params.quantdir = "./quantResults"

log.info """\
STAR - TEtranscripts
=============================
reads            : ${params.reads}
indexDir         : ${params.indexDir}
outdir           : ${params.outdir}
gtf              : ${params.gtf}
rmsk_ind         : ${params.rmsk_ind}
quantdir         : ${params.quantdir}
"""

//import modules
include { STARALIGN; TEcount } from './modules/STAR_TEtranscripts.nf'

workflow buildStarIndex {
    STARINDEX(indexDir_ch, fasta_ch, gtf_ch)
}

workflow {
    indexDir_ch = channel.fromPath( params.indexDir )
    gtf_ch = channel.fromPath( params.gtf )
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
    rmsk_ind_ch = channel.fromPath( params.rmsk_ind )
    STARALIGN( read_pairs_ch, indexDir_ch.toList() )
    TEcount( STARALIGN.out, gtf_ch.toList(), rmsk_ind_ch.toList() )
}

