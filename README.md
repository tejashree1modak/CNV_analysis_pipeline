## Exceptional copy number variation within and among populations of the eastern oyster (*Crassostrea virginica*) 
### Tejashree H. Modak, Rachel S. Schwartz

#### Reference genome and resequencing data:
- [Reference genome](https://www.ncbi.nlm.nih.gov/genome/398)
- Resequence data: 91 individuals of *C.virginica* across 16 locations along the eastern coast of the United States were sampled (details in the upcoming Genome paper).
#### CNV discovery
- **Software:** CNVs were called using 'germline SV calling' with default parameters in [DELLY2](https://github.com/dellytools/delly).
- **Requirements:** BAM files per sample and reference genome fasta. 
- All steps performed as described in the Germline SV calling section.
- Convert output file format from .bcf to .vcf using [bcftools](http://samtools.github.io/bcftools/bcftools.html#view). 

#### R scripts for data processing and analysis 
**Requirements:** Combined vcf output from DELLY for all samples. 

#### 1. Filtration

This step filters duplications identified by delly using the following criteria:

1. Duplications that are present in >90% of samples hence likely fixed in the population.   
2. Duplications that overlap >10% with a repeat region identified in the reference genome. 

Script: filter_dups.R

#### 2. Characterization

This step is used to characterize filtered duplications 

1. across locations
2. across genome

Script 1: post_filtration_analyses.R 

Script 2: cnv_analyses.R

#### 3. Annotation

This step annotates duplications using the reference genome annotation.

1. GO and KEGG annotation
2. Mapping duplications to different features in the genome

Script: annotation_dups.R
