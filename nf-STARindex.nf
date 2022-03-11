#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\
         STAR - building Index
         ============================
         """

// parameter
gtf_ch = Channel.fromPath( '/data/scratch/yaochung41/gtf_file/hg19.ensGene.gtf', checkIfExists: true )
fasta_ch = Channel.fromPath( '/data/scratch/yaochung41/reference/hg19.fa', checkIfExists: true )

//process
process starIndex {

    tag "Building star index"

    input:
    file gtf from gtf_ch
    file fasta from fasta_ch

    script:
    """
    STAR --runThreadN 20 \
         --runMode genomeGenerate \
         --genomeDir /data/scratch/yaochung41/genomeIndex/hg19 \
         --genomeFastaFiles $fasta \
         --sjdbGTFfile $gtf \
         --sjdbOverhang 100

    """

}