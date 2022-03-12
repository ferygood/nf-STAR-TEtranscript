#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Default pipeline parameters
 */

params.reads = '/data/scratch/yaochung41/kap1/data/hm/*_R{1,2}.fastq' 
params.index = '/data/scratch/yaochung41/genomeIndex/hg19'
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
    STARALIGN( read_pairs_ch, index_ch.toList() )
}
