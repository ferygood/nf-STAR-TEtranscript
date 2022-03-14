#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Default pipeline parameters
 */

params.reads = '/data/scratch/yaochung41/kap1/data/hm/*_R{1,2}.fastq' 
params.index = '/data/scratch/yaochung41/genomeIndex/hg19'
params.outdir = './results'
params.gtf = "/data/scratch/yaochung41/gtf_file/hg19.ensGene.gtf"
params.rmsk_ind = "/data/scratch/yaochung41/gtf_file/hg19_rmsk_TE.gtf.ind"
params.quantdir = "./quantResults"

log.info """\
STAR - TEtranscripts
=============================
fasta    : ${params.fasta}
index    : ${params.index}
reads    : ${params.reads}
outdir   : ${params.outdir}
gtf      : ${params.gtf}
rmsk_ind : ${params.rmsk_ind}
quantdir : ${params.quantdir}
"""

//import modules
include { STARALIGN; TEcount } from './modules/STAR_TEtranscripts.nf'

workflow {
    fasta_ch = channel.fromPath( params.fasta )
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
    index_ch = channel.fromPath( params.index )
    gtf_ch = channel.fromPath( params.gtf )
    STARINDEX( index_ch, fasta_ch, gtf_ch)
    STARALIGN( read_pairs_ch, index_ch.toList() )
    rmsk_ind_ch = channel.fromPath( params.rmsk_ind )
    TEcount( STARALIGN.out, gtf_ch.toList(), rmsk_ind_ch.toList() )
}

