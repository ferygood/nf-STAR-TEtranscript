#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process STARINDEX {

  file(params.indexDir).mkdir()
  publishDir "${params.indexDir}", mode: 'copy'

  input:
    path indexDir
    file fasta
    file gtf

  script:
    """
    STAR  --runThreadN 20 \
          --runMode genomeGenerate \
          --genomeDir $indexDir \
          --genomeFastaFiles $fasta \
          --sjdbGTFfile $gtf \
          --sjdbOverhang 100
    """

}


workflow {
    params.indexDir = '/data/scratch/yaochung41/genomeIndex/hg19'
    params.fasta = '/data/scratch/yaochung41/reference/hg19.fa'
    params.gtf = "/data/scratch/yaochung41/gtf_file/hg19.ensGene.gtf"
    
    indexDir_ch = channel.fromPath( params.indexDir )
    fasta_ch = channel.fromPath( params.fasta )
    gtf_ch = channel.fromPath( params.gtf )

    STARINDEX( indexDir_ch, fasta_ch, gtf_ch)
}