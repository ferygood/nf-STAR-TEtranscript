#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Default pipeline parameters
 */

params.reads = '/home/yaochung41/github/nf-STAR-TEtranscript/data/*_R{1,2}.fastq' 
params.index = '/home/yaochung41/starIndex/hg19'
params.outdir = './results'

log.info """\
STAR - TEtranscripts
=============================
index    : ${params.index}
reads    : ${params.reads}
outdir   : ${params.outdir}
"""

//import modules
include { STARALIGN } from './modules/test.nf'

workflow {
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
    index_ch = channel.fromPath( params.index )
    STARALIGN( read_pairs_ch, index_ch )
}