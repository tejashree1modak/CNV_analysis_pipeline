# This script uses filtered duplication data obtained from filtration for further characterization


# Dependencies installation and loading
for (i in c("UpSetR","tidyverse","here")) {
  if(!require(i, character.only=TRUE)) {
    install.packages(i, dependencies=TRUE)
    require(i, character.only=TRUE)
  }
}

#### Genome coverage ####
# Merged duplications to avoid counting overlapping duplication multiple times
# The output of filter_dups.R is used as input for bedtools
# bedtools merge -i cvir_filtered_dups.bed -c 1,2,3 -o count,collapse,collapse  > characterization/cvir_filtered_dups_merged.bed
dups_fil_merged <- read.table(here("characterization/cvir_filtered_dups_merged.bed"), 
                              sep="\t" , stringsAsFactors = FALSE)
colnames(dups_fil_merged) <- c("CHROM", "POS","end","count","POS_collapse","end_collapse")
# Length of merged dups
dups_fil_merged$len <- (dups_fil_merged$end - dups_fil_merged$POS) + 1
#Total number of bases of all duplications
sum(dups_fil_merged$len) #118732426
# % of genome covered by duplications
# Genome size of C.virginica genome can be found at https://www.ncbi.nlm.nih.gov/genome/398
(sum(dups_fil_merged$len)/684000000)*100 #17% 
# basic stats
median(dups_fil_merged$len)
mean(dups_fil_merged$len)
summary(dups_fil_merged$len)


#### Genotypes for Duplication  ####
#Extracting genotype will inform us of presence/absence of a duplication per sample
# Genotype 0/0 = homozygous for absence of duplication, genotype 0/1 or 1/0 = heterozygouse for presence of duplication
# genotype 1/1 = homozygous for presence of duplication
gett <- function(bedout_col){
  sep_out <- str_split( bedout_col, ':') %>% as.data.frame()
  g <- sep_out[1,] %>% unname() %>% t()
  g2<-g[,1]
  table(g2)
}
#run above function on all samples (columns)
gtypes <- select(oysterdup_fil,CL_1:UMFS_6) %>% map_dfr(gett)
gtypes2 <-as.data.frame(t(gtypes))
colnames(gtypes2) <- c("./.", "0/0",  "0/1",  "1/1")
gtypes2$pop <- map_chr(str_split(rownames(gtypes2),"_"),1) #add pop info (no sample num)
gtypes3 <- gather(gtypes2,genotype,number,-pop) #make long
#proportions
gtypes_p <- gtypes2 %>% mutate(sum=rowSums(select(gtypes2,("0/0":"1/1")))) %>%
  mutate(p0 = gtypes2$"0/0"/sum) %>% mutate(p01 = gtypes2$"0/1"/sum) %>%
  mutate(p1 = gtypes2$"1/1"/sum) %>% select(pop,p0,p01,p1)
gtypesp2 <- gather(gtypes_p,genotype,number,-pop)
gtypesp2$pop <- as.factor(gtypesp2$pop)
levels(gtypesp2$pop) <- c("OBOYS2","UMFS","NEH","DEBY","LOLA" ,
                          "HI","SM","HC","HCVA",  "CS", "CLP",
                          "SL","CL","LM",
                          "HG","NG")  #reorder pops
# Fig 2c from paper: Proportion of genotypes for duplications per population
ggplot(gtypesp2,aes(genotype,number,color=pop))+geom_boxplot() + labs(x="Genotype", y="Proportion") + theme_classic() +
  theme(axis.text.x  = element_text(size=14), axis.text.y  = element_text(size=14), 
        axis.title.x  = element_text(face = "bold", size=16), axis.title.y  = element_text(face = "bold", size=16),
        legend.title = element_blank())


#### Duplication lengths ####
oysterdup_fil <- read.table(here("characterization/oysterdup_fil"), 
              sep="\t" , stringsAsFactors = FALSE, header = TRUE)
# Fig 1a from paper: Frequency distribution of duplication lengths
ggplot(oysterdup_fil, aes(length))+geom_histogram(binwidth = 60,fill="steelblue")+ylim(c(0,100))+
  xlim(c(0,10000)) + labs(x="Length of duplications", y="Frequency") + theme_classic() +
  theme(axis.text.x  = element_text(size=18), axis.text.y  = element_text(size=18), axis.title.x  = element_text(face = "bold", size=20), axis.title.y  = element_text(face = "bold", size=20)) 
# Fig 1b from paper: Distribution of duplication lengths per population
pop_num_alts_present_fil <- read.table(here("characterization/pop_num_alts_present_fil"), 
              sep="\t" , stringsAsFactors = FALSE, header = TRUE)
ggplot(pop_num_alts_present_fil, aes(pop,length)) +geom_violin() + ylim(c(0,2500))
# Fig 1c from paper: Mean lengths of duplications compared across populations
meanl <- group_by(pop_num_alts_present_fil,pop) %>% summarize(mean_len = mean(length),sd = sd(length))
ggplot(meanl,aes(pop,mean_len))+geom_point()+
  geom_errorbar(aes(ymin=mean_len+sd,ymax=mean_len-sd))


#### Duplication comparison across locations ####
## Upset plot of duplications across locations POST FILTERATION ##
# get a de-duplicated list of locus id's
ids <- unique(pop_num_alts_present_fil$ID)
# for each id, get a binary indicator of whether it is in a pop and bind to one dataframe
pops <- unique(pop_num_alts_present_fil$pop)
binaries <- pops %>% 
  map_dfc(~ ifelse(ids %in% filter(pop_num_alts_present_fil, pop == .x)$ID, 1, 0) %>% 
            as.data.frame) # UpSetR doesn't like tibbles
# set column names
names(binaries) <- pops
# have a look at the data
head(binaries)  
# how many duplications are present in more than 3 locations
filter(binaries,rowSums(binaries)>3) %>% nrow() #6769 ie 6%
# Fig 2b from paper: UpSet plot of the intersected duplications across locations
upset(binaries, nsets = length(pops), main.bar.color = "SteelBlue", sets.bar.color = "DarkCyan", 
      sets.x.label = "Number duplicate loci", text.scale = c(rep(1.4, 5), 2), order.by = "freq")

# Dups per location POST FILTRATION #
pop_sum_fil <- as.data.frame(colSums(binaries))
pop_sum_fil <- data.frame(pop = names(binaries),total_dups=colSums(binaries))
#get proportion of duplications by:
# dividing total duplications per location by total number of duplications across locations
pop_sum_fil$prop <- pop_sum_fil$total_dups/length(oysterdup_fil$ID)  #number of filtered dups are 11339 
# Fig 2a from paper: Proportion of duplications per location
ggplot(pop_sum_fil, aes(x=pop,y=prop, color=pop)) + geom_bar(stat = "identity", fill="white") + 
  labs(x="Populations", y="Proportion of total duplications per population", title ="Post filteration") 
#+ scale_color_manual(values=values,labels=labels)

##Comparison of dups within samples for inbred pop ##
getg <- function(bedout_col){
  str_split( bedout_col, ':') %>% map_chr(1)
}

gtypes_only_fil <- map_dfr(select(oysterdup_fil,CL_1:UMFS_6),getg)
gtypes_only_fil$ID <- oysterdup_fil$ID
gtypes_long_fil <- gather(gtypes_only_fil,key=sample,value=gtype,-ID)
gtypes_long_fil$pop <- str_split(gtypes_long_fil$sample,'_') %>% map(1) %>% as.character()
gtypes_long_fil$pop <- as.vector(gtypes_long_fil$pop)
gtypes_long_fil$num_alts <- str_split(gtypes_long_fil$gtype,'/') %>% 
  map(as.integer) %>% 
  map_int(sum)
# For HG
gtypes_long_HG <- gtypes_long_fil %>% filter(num_alts > 0) %>% 
  filter(pop == 'HG')
ids_HG <- unique(gtypes_long_HG$ID)
samples <- unique(gtypes_long_HG$sample)
binaries_HG <- samples %>% 
  map_dfc(~ ifelse(ids_HG %in% filter(gtypes_long_HG, sample == .x)$ID, 1, 0) %>% 
            as.data.frame)
names(binaries_HG) <- samples
# For NG
gtypes_long_NG <- gtypes_long_fil %>% filter(num_alts > 0) %>% 
  filter(pop == 'NG')
ids_NG <- unique(gtypes_long_NG$ID)
samples_NG <- unique(gtypes_long_NG$sample)
binaries_NG <- samples_NG %>% 
  map_dfc(~ ifelse(ids_NG %in% filter(gtypes_long_NG, sample == .x)$ID, 1, 0) %>% 
            as.data.frame)
names(binaries_NG) <- samples_NG
# Fig 3b from paper: UpSet plot of the intersected duplications among inbred populations
upset(binaries_HG, nsets = length(samples), main.bar.color = "SteelBlue", sets.bar.color = "DarkCyan", 
      sets.x.label = "Number duplicate loci", text.scale = c(rep(1.4, 5), 2), order.by = "freq")
upset(binaries_NG, nsets = length(samples_NG), main.bar.color = "SteelBlue", sets.bar.color = "DarkCyan", 
      sets.x.label = "Number duplicate loci", text.scale = c(rep(1.4, 5), 2), order.by = "freq")

#### Frequency of Duplications per chromosome comparison across locations ####
# get chromosome locations
chrom_pos <- oysterdup_fil %>% select(CHROM, POS)
# get chromosome information including lengths (can be obtained from https://www.ncbi.nlm.nih.gov/genome/398)
chrom_len <- data.frame(CHROM=c("NC_035780.1","NC_035781.1","NC_035782.1","NC_035783.1","NC_035784.1","NC_035785.1",
                                "NC_035786.1", "NC_035787.1","NC_035788.1","NC_035789.1"), 
                        start=c(1,1,1,1,1,1,1,1,1,1), 
                        end=c(65668440,61752955,77061148,59691872,98698416,51258098,57830854,75944018,104168038,32650045))
chrom_len$len <- chrom_len$end - chrom_len$start
#get frequency of duplications per location
gtypes_pos_fil <- map_dfr(select(oysterdup_fil,CL_1:UMFS_6),getg)
gtypes_pos_fil$POS <- oysterdup_fil$POS
gtypes_pos_long_fil <- gather(gtypes_pos_fil,key=sample,value=gtype,-POS)
gtypes_pos_long_fil$pop <- str_split(gtypes_pos_long_fil$sample,'_') %>% map(1) %>% as.character()
gtypes_pos_long_fil$pop <- as.vector(gtypes_pos_long_fil$pop)
gtypes_pos_long_fil$num_alts <- str_split(gtypes_pos_long_fil$gtype,'/') %>% 
  map(as.integer) %>% 
  map_int(sum)
#adding dups in all individuals of same pop to give pop count
pop_num_pos_alts_fil <- gtypes_pos_long_fil %>% filter(!is.na(num_alts)) %>%
  group_by(pop,POS) %>%
  summarize(num_alts = sum(num_alts))
pop_num_pos_alts_present_fil <- filter(pop_num_pos_alts_fil,num_alts >0)
#join to get chromosome information 
pop_num_pos_alts_present_chrom_fil <- left_join(pop_num_pos_alts_present_fil, chrom_pos, by = "POS")
pop_alts_per_chrom_fil <- pop_num_pos_alts_present_chrom_fil %>% group_by(pop,CHROM) %>% 
  summarize(num_alts = sum(num_alts)) 
pop_alts_per_chrom_len_fil <- left_join(pop_alts_per_chrom_fil, chrom_len, by = "CHROM")
pop_alts_per_chrom_len_fil$pop <- factor (as.character(pop_alts_per_chrom_len_fil$pop), 
                                          levels=c("HI","SM","CS","HC","HCVA","CLP","CL","SL","LM","UMFS","NEH","HG","NG","DEBY","LOLA","OBOYS2"))

# Fig 4 from paper: Frequency of duplications per chromosome across locations normalized by chromosome length
ggplot(pop_alts_per_chrom_len_fil, aes(x=CHROM,y=(num_alts/len), color=pop)) + geom_bar(stat = "identity", fill="white") + 
  labs(x="Chromosome Number", y="Frequency of CNVs",title = "Post filteration")
#+ scale_color_manual(values=values,labels=labels)

