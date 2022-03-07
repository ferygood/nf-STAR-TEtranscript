#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\
         STAR - TEtranscript Pipeline
         ============================
         """

// Define parameters

process star_index {
    
    input:
      path gtf
      path fasta

    output:
      path 'index'
    
    script:
    """
    STAR --runThreadN 20 \
         --runMode genomeGenerate \
         --genomeDir /data/scratch/yaochung41/genomeIndex/hg19 \
         --genomeFastaFiles ${fasta} \
         --sjdbGTFfile ${gtf} \
         --sjdbOverhang 100

    """
}

process star_align {
    
    publishDir "results/quant", mode: 'symlink'

    input:
      tuple(val(sample_id), path(reads)) 
      each index

    output:
    tuple val(sample_id), path("${sample_id}")

    script:
    """
    STAR  --genomeDir ${index} \
          --runMode alignReads \
          --readFilesIn ${reads[0]} ${reads[1]} \
          --outFileNamePrefix ${sample_id} \
          --runThreadN 16 \
          --outSAMtype BAM SortedByCoordinate \
          --outFilterMultimapNmax 100 \
          --winAnchorMultimapNmax 100
    """

}


workflow {
    gtf_ch = Channel.fromPath( '/data/scratch/yaochung41/gtf_file/hg19.ensGene.gtf', checkIfExists: true )
    fasta_ch = Channel.fromPath( '/data/scratch/yaochung41/reference/hg19.fa', checkIfExists: true )
    reads_ch = Channel.fromFilePairs( '/data/scratch/yaochung41/kap1/data/hm/*_R{1,2}.fastq', checkIfExists: true )
    // index_ch = Channel.fromPath( '/data/scratch/yaochung41/genomeIndex/hg19' )
    index_ch = star_index(gtf_ch, fasta_ch)
    star_align(reads_ch, index_ch)
}

