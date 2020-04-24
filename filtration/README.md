## Exceptional copy number variation within and among populations of the eastern oyster (*Crassostrea virginica*) 
### Tejashree H. Modak, Rachel S. Schwartz

#### Duplication filtration:
#### Filter 1: Duplications overlapping with repeat regions identified in the reference genome
##### Duplications that overlap >10% with a repeat region were filtered from the analyses
- [Repeats identified in C.virginica genome](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/022/765/GCA_002022765.4_C_virginica-3.0)
Repeats file from NCBI was preprocessed to make a bedfile of repeats as follows:
Modify chromosome names to match gff3 to get Cvir_genome_repeats_mod.bed

```shell
 awk -v OFS='\t' 'NR>3{
    print $5,$6,$7,$5"_"$7
 }' GCA_002022765.4_C_virginica-3.0_rm.out  | \
 sed -e 's/CM008241.1/NC_035780.1/g' -e 's/CM008242.1/NC_035781.1/g' \
     -e 's/CM008243.1/NC_035782.1/g' -e 's/CM008244.1/NC_035783.1/g' \
     -e 's/CM008245.1/NC_035784.1/g -e 's/CM008246.1/NC_035785.1/g'  \
     -e 's/CM008247.1/NC_035786.1/g' -e 's/CM008248.1/NC_035787.1/g' \
     -e 's/CM008249.1/NC_035788.1/g' -e 's/CM008250.1/NC_035789.1/g' \
     > Cvir_genome_repeats_mod.bed
```

##### BEDTools used to obtain overlap between duplications and repeats
Use bedtools merge to merge repeats as follows:
bedtools merge -i Cvir_genome_repeats_mod.bed -c 1,5 -o count,collapse > Cvir_genome_repeats_merged.bed
Use bedtools using the -wo flag to obtain overlaps between merged repeats and dups
bedtools intersect -a oysterduplicate_sort.bed -b Cvir_repeats_merged.bed -wo > dup_repeat_merged_overlap_mod.bed
This file was used to filter duplications that overlap >10% with a repeat region.
Use filter_repeats.R for this step. 
