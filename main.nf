#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Default pipeline parameters
 */
params.fasta = '/data/scratch/yaochung41/reference/hg19.fa'
params.reads = '/data/scratch/yaochung41/kap1/data/hm/*_R{1,2}.fastq' 
params.indexDirPrefix = '/data/scratch/yaochung41/genomeIndex/hg19'
params.outdir = './results'
params.gtf = "/data/scratch/yaochung41/gtf_file/hg19.ensGene.gtf"
params.rmsk_ind = "/data/scratch/yaochung41/gtf_file/hg19_rmsk_TE.gtf.ind"
params.quantdir = "./quantResults"

log.info """\
STAR - TEtranscripts
=============================
fasta            : ${params.fasta}
reads            : ${params.reads}
indexDirPrefix   : ${params.indexDirPrefix}
outdir           : ${params.outdir}
gtf              : ${params.gtf}
rmsk_ind         : ${params.rmsk_ind}
quantdir         : ${params.quantdir}
"""

//import modules
include { STARINDEX; STARALIGN; TEcount } from './modules/STAR_TEtranscripts.nf'

workflow {
    indexDirPrefix_ch = channel.fromPath( params.indexDirPrefix )
    fasta_ch = channel.fromPath( params.fasta )
    gtf_ch = channel.fromPath( params.gtf )
    STARINDEX( indexDirPrefix_ch, fasta_ch, gtf_ch)
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
    index_ch = channel.fromPath( params.indexDirPrefix )
    STARALIGN( read_pairs_ch, index_ch.toList() )
    rmsk_ind_ch = channel.fromPath( params.rmsk_ind )
    TEcount( STARALIGN.out, gtf_ch.toList(), rmsk_ind_ch.toList() )
}

