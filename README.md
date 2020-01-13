## Exceptional copy number variation within and among populations of the eastern oyster (*Crassostrea virginica*) 
### Tejashree H. Modak, Rachel S. Schwartz

#### Reference genome and resequencing data:
- [Reference genome](https://www.ncbi.nlm.nih.gov/genome/398)
- Resequence data: 91 individuals of *C.virginica* across 16 locations along the eastern coast of the United States were sampled (details in the upcoming Genome paper).
#### CNV discovery
- **Software:** CNVs were called using default parameters in [DELLY2](https://github.com/dellytools/delly).
- **Requirements:** BAM files per sample and reference genome fasta. 
- All steps performed as described in the Germline SV calling section.
- Convert output file format from .bcf to .vcf using [bcftools](http://samtools.github.io/bcftools/bcftools.html#view). 

#### R scripts
##### R scripts for data processing and analysis 
**Requirements:** Combined vcf output from DELLY for all samples.  
