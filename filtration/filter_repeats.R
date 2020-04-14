#This script is a part of filtration of duplicatons that overlap with repeat regions
# Duplications that overlap >10% with a repeat region were filtered from the analyses. 

# Dependencies installation and loading
for (i in c("tidyverse" , "here")) {
  if(!require(i, character.only=TRUE)) {
    install.packages(i, dependencies=TRUE)
    require(i, character.only=TRUE)
  }
}

#Repeatmasker output was obtained from the ftp server where C.virginica genome files are at
#ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/022/765/GCA_002022765.4_C_virginica-3.0
#The file is GCA_002022765.4_C_virginica-3.0_rm.out.gz
#Input file generated from NCBI as preprocessed to make a bedfile of repeats as follows:
#Input file 1: awk -v OFS='\t' '{print $6, $7}' GCA_002022765.4_C_virginica-3.0_rm.out > Cvir_genome_repeats.txt
#Input file 2: awk -v OFS='\t' 'NR>3{print $5,$6,$7,$5"_"$7}' GCA_002022765.4_C_virginica-3.0_rm.out >  Cvir_genome_repeats.bed

#Modify chromosome names to match gff3 and get Cvir_genome_repeats_mod.bed
#sed 's/CM008241.1/NC_035780.1/g' , sed 's/CM008242.1/NC_035781.1/g', sed 's/CM008243.1/NC_035782.1/g', 
#sed 's/CM008244.1/NC_035783.1/g', sed 's/CM008245.1/NC_035784.1/g, sed 's/CM008246.1/NC_035785.1/g', 
#sed 's/CM008247.1/NC_035786.1/g', sed 's/CM008248.1/NC_035787.1/g', sed 's/CM008249.1/NC_035788.1/g', 
#sed 's/CM008250.1/NC_035789.1/g',
  
#Get overlap 
#Use bedtools merge to merge repeats
##bedtools merge -i Cvir_genome_repeats_mod.bed -c 1,5 -o count,collapse > Cvir_genome_repeats_merged.bed
#Use bedtools using the -wo flag to obtain overlaps between merged repeats and dups
#bedtools intersect -a oysterduplicate_sort.bed -b Cvir_repeats_merged.bed -wo > dup_repeat_merged_overlap_mod.bed

#all vcf data for each individual for each duplication obtained from DELLY
oysterdup <- read.table(here("filtration/germline_nohead_dup.vcf"),stringsAsFactors = FALSE)
header <- strsplit("CHROM POS ID      REF     ALT     QUAL    FILTER  
                   INFO    FORMAT  CL_1    CL_2    CL_3    CL_4    CL_5    CL_6    CLP_1   CLP_2   
                   CLP_3   CLP_4   CLP_5   CLP_6   CS_1    CS_2    CS_3    CS_5    CS_6    CS_7    
                   DEBY_1  DEBY_2  DEBY_3  DEBY_4  DEBY_5  DEBY_6  HC_1    HC_3    HC_4    HC_5    
                   HC_6    HC_7    HCVA_1 HCVA_2 HCVA_3 HCVA_4 HCVA_5 HCVA_6 HG_HG0F2       
                   HG_HG2F1        HG_HG2M5
                   HI_1    HI_2    HI_3    HI_4    HI_5    HI_6    LM_1_pool       LM_3    LM_4    
                   LM_7    LM_8    LOLA_1  LOLA_2  LOLA_3  LOLA_4  LOLA_5  LOLA_6  NEH_1   NEH_2   
                   NEH_3   NEH_4   NEH_5   NEH_6   NG_NH0H4        NG_NH2F6        NG_NH2F8        
                   NG_NH2M1        OBOYS2_1        OBOYS2_2        OBOYS2_3        OBOYS2_4        
                   OBOYS2_5        OBOYS2_6        SL_1    SL_2    SL_3    SL_4    SL_5    SL_6
                   SM_10   SM_11   SM_12   SM_7    SM_8    SM_9    UMFS_1  UMFS_2  UMFS_3  UMFS_4  
                   UMFS_5  UMFS_6", "\\s+")[[1]]
colnames(oysterdup)<-header
oysterdup <-dplyr::filter(oysterdup,FILTER=="PASS")
oysterdup$end <- str_split(oysterdup$INFO, ';') %>%
  map_chr(5) %>% str_split('=') %>% map_chr(2) %>% as.integer()
oysterdup$length <- oysterdup$end - oysterdup$POS

# Read in bedtools Ouput of intersect between repeat regions in reference genome and duplications  
dup_repeat_overlap <- read.table(here("filtration/dup_repeat_merged_overlap_mod.bed"), 
                                 sep="\t" , stringsAsFactors = FALSE)
colnames(dup_repeat_overlap) <- c("CHROM", "POS","end","ID","R_POS","R_end","R_ID","l")
#Number of repeats mapped to each duplicate
dup_repeat_overlap %>% select("ID","l") %>% group_by(ID) %>% tally() %>% View()
#Total len of repeats included in dups
overlap_total_len <- dup_repeat_overlap %>% select("ID","l") %>% group_by(ID) %>% summarise(sum(l))
colnames(overlap_total_len) <- c("ID","total_len")
#Get % of dup len covered by repeats for each dup that overlaps with repeats
percent_overlap <- oysterdup %>% select(ID,length) %>% left_join(overlap_total_len,by = 'ID') %>% na.omit() 
percent_overlap$percent <- (percent_overlap$total_len/percent_overlap$length)*100
#dups with >10% repeat coverage
percent_overlap %>% filter(percent > 10) %>% nrow() #filter out 1778 dups
repeat_filter_dups <- percent_overlap %>% filter(percent > 10) %>% select("ID")
# These duplications were filtered from the final set. 