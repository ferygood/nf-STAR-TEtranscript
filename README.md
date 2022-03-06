# nf-STAR-TEtranscript 

## Introduction  

This pipeline is built using [Nextflow](https://www.nextflow.io/). The purpose is to handle RNA-seq fastq file from aligning using [STAR](https://github.com/alexdobin/STAR) and estiamte their transposable elements expression using [TEtranscript](https://github.com/mhammell-laboratory/TEtranscripts). 

## Quick Start  
1. Install `Nextflow` (>= 21.10.6)
2. We recommend using this pipeline with `Singularity`.
3. Downlaod the pipeline and test it with built-in example dataset:  
```sh
nextflow run nf-core/star-tetranscript -profile test,YOURPROFILE --outdir <OUTDIR>
```
4. Start running your analysis:

## Documentation. 
The nf-core pipeline comes with documentation about the pipeline usage, parameters and output.
 
## Credits  
 
## Contributions and Supports  

## References  
- Nextflow enables reproducible computational workflows (DOI: 10.1038/nbt.3820)
- STAR: ultrafast universal RNA-seq aligner (DOI: 10.1093/bioinformatics/bts635)
- TEtranscripts: a package for including transposable elements in differential expression analysis of RNA-seq datasets (DOI: 10.1093/bioinformatics/btv422)
