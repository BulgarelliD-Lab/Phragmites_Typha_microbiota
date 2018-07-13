#############################################################
#
# Ref to the ARTICLE
# 
# Laura Pietrangelo and Davide Bulgarelli
# 
# laa.piee@gmail.com
# d.bulgarelli@dundee.ac.uk
# 
# Revison July 2018
# 
# script to reproduce calculations and figures presented in the manuscript
# 
# Disclaimer: the manuscript is currently submitted, the script might be subjected to changes derived from the revision process
#
#############################################################


#############################################################
# Libraries and functions required
#############################################################

#These initial commands are required to clean-up the memory and start a new session
rm(list=ls())
dev.off()


#set working directory
#From Laura CPU
setwd("C:/Users/LaUrA/Desktop/JH06_L/")
#From Davide CPU
#setwd("C:/Users/DB42008/Box Sync/Davide_lab_manuscripts/Laura_TyphaPhragmites_2018/R_revision/")

#1st time Phyloseq installation
#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
#biocLite("DESeq2")
#biocLite("VennDiagram")
#biocLite("PMCMR")

#load the required packages (all included in Phyloseq, but is necessary to invoke them)
library("phyloseq")
library("DESeq2")
library("ggplot2")
library("vegan")
library ("ape")
library("VennDiagram")
library("PMCMR")
library("ade4")


#############################################################

#############################################################
#import the count matrix and the desing file

#OTU table this file has been generated using QIIME 1.9.0. In the OTU ids, OTU abundance information has been removed
dat_info <- read.delim("JH06_otu_table_nc2.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
colnames(dat_info)

#inspect the file using the command dim
dim(dat_info)

#extract the total number of reads clustered at OTU 97% identiy (the number underneath the Hv identifier represents the total number of reads clustered for that sample) 
OTU_97_reads <- sort(colSums(dat_info[, 1:40]))
OTU_97_reads

#total reads
OTU_97_reads_sum <- sum(colSums(dat_info[, 1:40]))
OTU_97_reads_sum

#design file
design <- read.delim("Map_JH06_detailed_phyloseq_L.txt", sep = "\t", header=TRUE, row.names=1)
design

#remove chloroplast e mitochondria OTUs from the original dataset
#to achieve this task we use the command grepl that creates two new datasets: one with chloroplast only and the other one with mitochondria only OTUs
Chloroplast <- dat_info[grepl("Chloroplast", dat_info$ConsensusLineage), ]
dim(Chloroplast)

mitochondria <- dat_info[grepl("mitochondria", dat_info$ConsensusLineage), ]
dim(mitochondria)


#set a difference between the row names of the the three dataset: this information will be used to filter Plant derived OTUs from the phyloseq object
noPlants <- setdiff(rownames(dat_info), c(rownames(Chloroplast), rownames(mitochondria)))
#and check the results
length(rownames(dat_info))
length(noPlants)

#save this piece of information for qiime purposes. We need to create a a new OTU table in QIIME to generate the taxa tables
write(noPlants, "JH06_noPlant_OTUs_id.txt")

#and generate a new OTU table depleted with chloroplast and mitochondria OTUs (for comparison purposes)
dat_info_noPlants <- dat_info[noPlants, ]

dat_info_noPlants

#create a new count matrix without OTUs assigned to Choloplast and Mitochondria
dat_count <- dat_info[, rownames(design)]
dat_count_noplants <- dat_info[noPlants, rownames(design)]
dim(dat_count_noplants)

#and a new taxa table
dat_tax_noPlants <- as.data.frame(dat_info[rownames(dat_count_noplants), 41])
rownames(dat_tax_noPlants) <- rownames(dat_count_noplants)
#now we need to save the above file and in excel generate a tax table where column represents a taxonomic rank
write.table(dat_tax_noPlants, file="JH06_dat_tax_noPlants.txt", sep="\t")


#check the effect of mitochondria/chloroplast depletion on the new OTU table

#with plant sequences
dim(dat_count)

#w/o plant sequences
dim(dat_count_noplants)

#total number of reads w/o plant sequences
OTU_97_reads_noPlants <- colSums(dat_count_noplants)
OTU_97_reads_noPlants

#now sort per sample
sort(OTU_97_reads_noPlants)

#total number of reads
OTU_97_reads_noPlants_sum <- sum(OTU_97_reads_noPlants)
OTU_97_reads_noPlants_sum 

#now determine the proportion of non-plant reads in the original dataset (this might be a recurrent question: recall Tim Mauchline concerns, which is the proportion of contaminant sequences (e.g., plant-derived) in your datases?)
useful_reads <- (OTU_97_reads_noPlants_sum/sum(colSums(dat_count)))*100
useful_reads

#determine total amount of reads of sample included in Sample set 1 and 2 which are included in the ms analysis
set1_samples <- rownames(design)[which(design$selection == 1)]
set2_samples <- rownames(design)[which(design$replicate ==1)]
set1_2_samples <- unique(c(set1_samples,set2_samples))
#total reads
set1_2_samples_reads <- sum(colSums(dat_count[, set1_2_samples]))
set1_2_samples_reads
#total reads w/o plants 
set1_2_samples_reads_noPlants <- sum(colSums(dat_count_noplants[, set1_2_samples]))
set1_2_samples_reads_noPlants
#% of retained reads
retained_reads <- (set1_2_samples_reads_noPlants/set1_2_samples_reads)*100
retained_reads
#############################################################

#############################################################
#generate a phyloseq object

#a) The OTU Table counts
JH06_OTU <- otu_table(dat_count_noplants, taxa_are_rows=TRUE)
#verify the dimension of the files
dim(dat_count_noplants)
dim(otu_table(JH06_OTU))

#b) The taxonomy information
#Note this is a new file generated in excell from the output of the command of lines 93-97
#it is a tab-delimited file with 8 columns, the header names are: OTU id; "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",
#and the empty cells are filled with the term 'Unassigned'
JH06_taxa_ordered <- read.delim ("JH06_dat_tax_noPlants_ordered.txt", sep = "\t", row.names=1, header=T, blank.lines.skip = FALSE)
JH06_taxa <- tax_table(as.matrix(JH06_taxa_ordered))
dim(JH06_taxa)

#c) The mapping file 
JH06_map <- sample_data(design)

#d) The phylogenetic tree: the OTU table has been generated using a closed reference approach agains the greengenes 13_05 database use the corresponding phylogenetic tree 
JH06_tree <- read_tree_greengenes("97_otus.tree.gz")

#check whether the tree is rooted
is.rooted(JH06_tree)

#merge the files and create the phyloseq object
JH06_data_phyloseq <- merge_phyloseq(JH06_OTU, JH06_taxa, JH06_map,  JH06_tree)

#inspect the generated data
JH06_data_phyloseq
sum(colSums(otu_table(JH06_data_phyloseq)))
dim(dat_count_noplants)
sum(colSums(dat_count_noplants))

#############################################################

##########
#Identify contaminating OTUs
##########

#first step
#remove samples that have less than 10,000 reads
JH06_data_phyloseq_2 <- prune_samples(sample_sums(JH06_data_phyloseq)>=10000, JH06_data_phyloseq)
JH06_data_phyloseq_2

#inspect the dataset using betadiversity
#Transform the count in relative abundance cpm
JH06_data_phyloseq_prop <- transform_sample_counts(JH06_data_phyloseq_2,  function(x) 1e+03 * x/sum(x))

#PCoA Bray distance
JH06_data_phyloseq_prop_bray <- ordinate(JH06_data_phyloseq_prop, "PCoA", "bray")
plot_ordination(JH06_data_phyloseq_prop, JH06_data_phyloseq_prop_bray , color = "Microhabitat")

#assign shapes to Species and color to Microhabitat
p=plot_ordination(JH06_data_phyloseq_prop, JH06_data_phyloseq_prop_bray , shape ="Species", color = "Microhabitat")
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_colour_manual(values = c("red","black", "darkblue", "purple"))
p + ggtitle("PCoA 16S data, Bray distance")

#PCoA weighted unifrac distance
JH06_data_phyloseq_prop_wunifrac <- ordinate(JH06_data_phyloseq_prop, "PCoA", "unifrac", weighted = TRUE)
plot_ordination(JH06_data_phyloseq_prop, JH06_data_phyloseq_prop_wunifrac , color = "Microhabitat")

#assign shapes to Species and color to Microhabitat
p=plot_ordination(JH06_data_phyloseq_prop, JH06_data_phyloseq_prop_wunifrac , shape ="Species", color = "Microhabitat")
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_colour_manual(values = c("red","black", "darkblue", "purple"))
p + ggtitle("PCoA 16S data, weighted Unifrac distance")

#there is a clear separation between samples and controls
#assess the significance of these differences
#WU distance
WU <- phyloseq::distance(JH06_data_phyloseq_prop, "unifrac", weighted= TRUE)
#inspect the dissimilarity matrix: this is a file where differences in termo of OTUs abundances and presence is converted in a coefficient ranging from 0 (identical samples) to 1 (none of the features=OTUs conserved)
WU
#statistical test
adonis(WU ~ Description, data= design[row.names(sample_data(JH06_data_phyloseq_prop)), ], permutations = 5000)

#test 1: remove OTUs that have at least 0.01% in the negative control
#create a phyloseq object with the two sets of samples
#Controls
JH06_data_phyloseq_controls <- subset_samples(JH06_data_phyloseq_prop, Description =="control")
JH06_data_phyloseq_controls

#Filter OTUs from controls: let's consider OTUs with a relative abundance above 0.01% as contaminat, If you change the number 1 below it might affect the impact on the dataset
JH06_data_phyloseq_controls_01 = prune_taxa(rowMeans(otu_table(JH06_data_phyloseq_controls)) > 0.1, JH06_data_phyloseq_controls)
JH06_data_phyloseq_controls_01
#identify the contaminant OTUs
contaminat_OTUs <- row.names(otu_table(JH06_data_phyloseq_controls_01))
#inspect the taxonomy information
dat_tax_noPlants[contaminat_OTUs, ]

#amount of reads assigned to "contaminat OTUs"
sample_sums(JH06_data_phyloseq_controls_01)

#not bad we got more than 99% of the reads assigned to controls

#filter the OTUs from the initial dataset of the samples
#not the code if you use '==' in the option of the subset_samples you will keep all the samples flagged by that attribute (see line 312 above). 
#Conversely, if you use the "!=" option you will keep everything which is not flagged by that attribute
JH06_data_phyloseq_3_samples <- subset_samples(JH06_data_phyloseq_2, Description != "control")
JH06_data_phyloseq_3_samples_no_contaminant <- prune_taxa(setdiff(rownames(otu_table(JH06_data_phyloseq_3_samples)), contaminat_OTUs), JH06_data_phyloseq_3_samples)   
#check if the two datasets differ of 56 OTUs
JH06_data_phyloseq_3_samples
JH06_data_phyloseq_3_samples_no_contaminant

#effect of contamination on number of reads
sum(sample_sums(JH06_data_phyloseq_3_samples_no_contaminant))/sum(sample_sums(JH06_data_phyloseq_3_samples))*100

#use betadiversity to inspect the samples only
#transform in relative abundance
JH06_data_phyloseq_3_samples_no_contaminant_prop <- transform_sample_counts(JH06_data_phyloseq_3_samples_no_contaminant,  function(x) 1e+03 * x/sum(x))

#PCoA Bray distance
JH06_samples_no_contaminant_bray <- ordinate(JH06_data_phyloseq_3_samples_no_contaminant_prop, "PCoA", "bray")
plot_ordination(JH06_data_phyloseq_3_samples_no_contaminant_prop, JH06_samples_no_contaminant_bray, color = "Microhabitat")

#assign shapes to Species and color to Microhabitat
p=plot_ordination(JH06_data_phyloseq_3_samples_no_contaminant_prop, JH06_samples_no_contaminant_bray, shape ="Microhabitat", color = "Species")
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_colour_manual(values = c("red", "darkblue"))
p + ggtitle("PCoA 16S data, Bray distance")

#PCoA weighted unifrac distance
JH06_samples_no_contaminant_wunifrac <- ordinate(JH06_data_phyloseq_3_samples_no_contaminant_prop, "PCoA", "unifrac", weighted = TRUE)
plot_ordination(JH06_data_phyloseq_3_samples_no_contaminant_prop, JH06_samples_no_contaminant_wunifrac, color = "Microhabitat")

#assign shapes to Species and color to Microhabitat
p=plot_ordination(JH06_data_phyloseq_3_samples_no_contaminant_prop, JH06_samples_no_contaminant_wunifrac, shape ="Microhabitat", color = "Species")
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_colour_manual(values = c("red", "darkblue"))
p + ggtitle("PCoA 16S data, weighted Unifrac distance")

#assess the significance of these differences
BC <- phyloseq::distance(JH06_data_phyloseq_prop, "bray")
WU <- phyloseq::distance(JH06_data_phyloseq_prop, "unifrac", weighted= TRUE)
#statistical test
adonis(BC ~ Microhabitat * Species, data= design[row.names(sample_data(JH06_data_phyloseq_prop)), ], permutations = 5000)
adonis(WU ~ Microhabitat * Species, data= design[row.names(sample_data(JH06_data_phyloseq_prop)), ], permutations = 5000)



######################## #######################
########################## samples_set1 ##########################

#OK we should remove the technical replicates: we could use them to define the minimum number of reads to be retained for an OTU to be considered reproducible
#firstly let's see the results with the replicates 1 flagged in the map file as selection 
JH06_data_phyloseq_3_samples_no_contaminant_T1 <- subset_samples(JH06_data_phyloseq_3_samples_no_contaminant, selection == 1)

#how many reads have we retained?
sort(sample_sums(JH06_data_phyloseq_3_samples_no_contaminant_T1))

#Remove OTUs not seen more than 25 times in at least 20% of the samples
JH06_data_phyloseq_4 = filter_taxa(JH06_data_phyloseq_3_samples_no_contaminant_T1, function(x) sum(x > 25) > (0.2*length(x)), TRUE)
JH06_data_phyloseq_4 

#how many reads are left?
sample_sums(JH06_data_phyloseq_4)
sort(sample_sums(JH06_data_phyloseq_4))

#let's rarefy the dataset at 66,000 reads
JH06_data_phyloseq_66K_set1 <- rarefy_even_depth(JH06_data_phyloseq_4 , sample.size = 66000)

#extract the data
JH06_data_phyloseq_66K_set1_table <- as.data.frame(otu_table(JH06_data_phyloseq_66K_set1))
dim(JH06_data_phyloseq_66K_set1_table)

#save the file for the reproducibility of the code
write.table(JH06_data_phyloseq_66K_set1_table, file="JH06_data_phyloseq_66K_set1_table.txt", sep="\t")

#import the rarefied OTU counts (note file name counts2.txt this file has been generated in excel and includes the #OTU ID as header of the first column)
dat_count_JH06_66K_set1 <- read.delim("JH06_data_phyloseq_66K_set1_table.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
dim(dat_count_JH06_66K_set1)

#Note: the above table will be used to generate Supplementary Figure 1-set 1 in the revised manuscript
#generate a PICRUSt table containing all OTUs
#prune OTU table below the threshold from the taxonomy information
dat_tax_noPlants_66K <- as.data.frame(dat_tax_noPlants[rownames(dat_count_JH06_66K_set1), ])
colnames(dat_tax_noPlants_66K) <- c("ConsensusLineage") 
colnames(dat_tax_noPlants_66K)
#merge count and taxonomy files
dat_info_rare_66K <- cbind(dat_count_JH06_66K_set1, dat_tax_noPlants_66K)
dim(dat_info_rare_66K)
colnames(dat_info_rare_66K)

#generate a design file set 1 to generate taxa chart including set 1 samples only
design_set1 <- design[colnames(dat_count_JH06_66K_set1), ]
design_set1_qiime <- as.data.frame(design_set1[, 1:4])
rownames(design_set1_qiime)
colnames(design_set1_qiime)
dim(design_set1_qiime)

#save the files for the reproducibility of the code
write.table(dat_info_rare_66K, file="JH06_info_rare_66K_set1.txt", sep="\t")
#Note: for taxonomy calculation in QIIME the file JH06_info_rare_66K_set1_2.txt, with the header of the first column '#OTU ID' was generated in excel
write.table(design_set1_qiime, file="JH06_map_set1_qiime.txt", sep="\t")
#Note: for taxonomy calculation in QIIME the file JH06_map_set1_qiime_2.txt, with the header of the first column '#SampleID' was generated in excel

#Import average value phylum level
dat_info_taxa_Phylum_set_1 <- read.delim("Description_otu_table_L2_set1.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
colnames(dat_info_taxa_Phylum_set_1)
dim(dat_info_taxa_Phylum_set_1)

# transform the data in % for visualisation
dat_norm_Phylum_set_1 <- dat_info_taxa_Phylum_set_1 * 100 
dat_norm_Phylum_set_1[1:5, ]
colSums(dat_norm_Phylum_set_1)

# determine the average % of reads for each phylum
Phylum_mean_sorted_set_1 <- dat_norm_Phylum_set_1[(order(-rowSums(dat_norm_Phylum_set_1))), ] 

#Identifythe phyla whose average abundance in above 1% and their aggregated relative abundance
Phylum_mean_topRank_set_1 <- Phylum_mean_sorted_set_1[rownames(Phylum_mean_sorted_set_1)[which(rowMeans(Phylum_mean_sorted_set_1) > 1)], ]
dim(Phylum_mean_topRank_set_1)
colSums(Phylum_mean_topRank_set_1)
rowMeans(Phylum_mean_topRank_set_1)

#transform in matrix for plotting purposes
Phylum_mean_topRank_set_1 <- as.matrix(Phylum_mean_topRank_set_1) 
colnames(Phylum_mean_topRank_set_1) 

#re-arrange the order for plotting purposes
col.order <- c("waterP","waterhizoP","biofilmP","waterT","waterhizoT","biofilmT")
Phylum_mean_topRank_set_1 <- Phylum_mean_topRank_set_1[,col.order]

# Stacked Bar Plot with Colors and Legend
barplot(Phylum_mean_topRank_set_1, main="Phylum Distribution",
        xlab="Samples", ylab = "Reads %", ylim = c(0,100), col=c("darkblue","red", "green", "cyan", "orange", "brown", "black", "magenta", "grey"), beside=FALSE,   legend = rownames(Phylum_mean_topRank_set_1))

#Due to size limits the legend covers part of the graph, save the graph as .eps file and in illustrator uou can adjust this (and other)  graphical issues
#barplot_no legend: you can save this as image and then reconstruct the leged using the previous figure. Note: the order of samples can be inferred from the command on line 256
barplot(Phylum_mean_topRank_set_1, main="Phylum Distribution",
        xlab="Samples", ylab = "Reads %", ylim = c(0,100), col=c("darkblue","red", "green", "cyan", "orange", "brown", "black", "magenta", "grey"), beside=FALSE)

#convert the file into a phyloseq-compatible OTU table
JH06_OTU_5 <- otu_table(dat_count_JH06_66K_set1, taxa_are_rows=TRUE)

#create a new phyloseq object
#Note: the number of samples/OTUs in the otu table written drive the composition of the phyloseq object
#For this reason we could use older mapping, taxonomy and phylogenetic tree: they will be "pruned" of the OTUs/samples not present in the OTU table
JH06_data_phyloseq_5 <- merge_phyloseq(JH06_OTU_5, JH06_taxa, JH06_map,  JH06_tree)

#inspect the created object
JH06_data_phyloseq_5

################################################################################################
#OK we will use the JH06_data_phyloseq_5 for the analysis and figure generation in the manuscript
###############################################################################################

#first thing to check is alpha diversity
#Index calculations
JH06_alpha_rare <-  estimate_richness(JH06_data_phyloseq_5, measures = c("Observed", "Chao1", "Shannon"))

#generate a new data frame combining alphadiversity calculation and sample information for data visualisation
#Sample information

#Microhabitat
design_microhabitat <- as.data.frame(sample_data(JH06_data_phyloseq_5)[, 1])
rownames(design_microhabitat) <- rownames(sample_data(JH06_data_phyloseq_5))
colnames(design_microhabitat) <- c("Microhabitat")

#Species
design_Species <- as.data.frame(sample_data(JH06_data_phyloseq_5)[, 2])
rownames(design_Species) <- rownames(sample_data(JH06_data_phyloseq_5))
colnames(design_Species) <- c("Species")

#combine the file
design_MS <- cbind(design_microhabitat, design_Species)

##############
#Observed OTUs
#############

JH06_alpha_rare_rare_Observed <- as.data.frame(JH06_alpha_rare[ ,1])
rownames(JH06_alpha_rare_rare_Observed) <- rownames(JH06_alpha_rare)
colnames(JH06_alpha_rare_rare_Observed) <- c("Observed")

#Combine the dataset sample description and Observed OTUs
JH06_alpha_rare_rare_Observed_MS <- cbind(design_MS, JH06_alpha_rare_rare_Observed)
JH06_alpha_rare_rare_Observed_MS

#Order the levels according to a defined order
JH06_alpha_rare_rare_Observed_MS$Microhabitat <- ordered(JH06_alpha_rare_rare_Observed_MS$Microhabitat, levels=c("water","waterhizo","biofilm")) 
JH06_alpha_rare_rare_Observed_MS$Species <- ordered(JH06_alpha_rare_rare_Observed_MS$Species, levels=c("P","T")) 

#################################################
#Figure Observed OTUs
with(JH06_alpha_rare_rare_Observed_MS, boxplot(Observed ~ Microhabitat * Species, xlab = "Samples", ylab = "Observed OTUs",   main = "OTU richness", ylim = c(800, 1800)))
#################################################


#Observed OTUs set1 stripchart
stripchart(Observed ~ Microhabitat * Species,
           data=JH06_alpha_rare_rare_Observed_MS,
           main="OTUs richness",
           xlab="Samples",
           ylab="Observed OTUs",
           col="blue", 
           group.names=c("water.P","waterhizo.P","biofilm.P","water.T","waterhizo.T","biofilm.T"),
           vertical=TRUE,
           pch=16)


#statistical analysis conducted on biofilm and waterrhizo only
#remove water samples from the analysis
nowater_set1 <- rownames(JH06_alpha_rare_rare_Observed_MS)[which(JH06_alpha_rare_rare_Observed_MS$Microhabitat!= "water")]
JH06_alpha_rare_rare_Observed_MS_nowater_set1 <- JH06_alpha_rare_rare_Observed_MS[nowater_set1, ]

#test the normality of the observed data
shapiro.test(JH06_alpha_rare_rare_Observed_MS_nowater_set1$Observed)


# data not normally distributed (p-value < 0.05); use a non parametric test to investigate the microhabitat effect
wilcox.test(Observed ~ Microhabitat, data = JH06_alpha_rare_rare_Observed_MS_nowater_set1)

#since the pValue from Shapiro test is close to the alpha level, use also the two-way anova to asses significant differences between samples but simply as check
Observed_stats_set1 <- aov(JH06_alpha_rare_rare_Observed_MS_nowater_set1$Observed ~ Microhabitat* Species, data = JH06_alpha_rare_rare_Observed_MS_nowater_set1)
summary(Observed_stats_set1)


##############
#Chao 1
#############

JH06_alpha_rare_rare_Chao1 <- as.data.frame(JH06_alpha_rare[ ,2])
rownames(JH06_alpha_rare_rare_Chao1) <- rownames(JH06_alpha_rare)
colnames(JH06_alpha_rare_rare_Chao1) <- c("Chao1")

#Combine the dataset sample description and Chao1 OTUs
JH06_alpha_rare_rare_Chao1_MS <- cbind(design_MS, JH06_alpha_rare_rare_Chao1)
JH06_alpha_rare_rare_Chao1_MS

#Order the levels according to a defined order
JH06_alpha_rare_rare_Chao1_MS$Microhabitat <- ordered(JH06_alpha_rare_rare_Chao1_MS$Microhabitat, levels=c("water","waterhizo","biofilm")) 
JH06_alpha_rare_rare_Chao1_MS$Species <- ordered(JH06_alpha_rare_rare_Chao1_MS$Species, levels=c("P","T")) 


#############################################################
#Figure Chao1 OTUs
with(JH06_alpha_rare_rare_Chao1_MS, boxplot(Chao1 ~ Microhabitat * Species, xlab = "Samples", ylab = "Chao1 OTUs",   main = "OTU richness", ylim = c(900, 1900)))
#############################################################


#Chao1 set1 stripchart
stripchart(Chao1 ~ Microhabitat * Species,
           data=JH06_alpha_rare_rare_Chao1_MS,
           main="OTUs richness",
           xlab="Samples",
           ylab="Chao1 OTUs",
           col="darkgreen",
           group.names=c("water.P","waterhizo.P","biofilm.P","water.T","waterhizo.T","biofilm.T"),
           vertical=TRUE,
           pch=16)

#statistical analysis conducted on biofilm and waterrhizo only
#remove water samples from the analysis
nowater_Chao1_set1 <- rownames(JH06_alpha_rare_rare_Chao1_MS)[which(JH06_alpha_rare_rare_Chao1_MS$Microhabitat!= "water")]
JH06_alpha_rare_rare_Chao1_MS_nowater_set1 <- JH06_alpha_rare_rare_Chao1_MS[nowater_Chao1_set1, ]

#test the normality of the Chao1 data
shapiro.test(JH06_alpha_rare_rare_Chao1_MS_nowater_set1$Chao1)

# data not normally distributed (p-value < 0.05); use a non parametric test to investigate the microhabitat effect
wilcox.test(Chao1 ~ Microhabitat, data = JH06_alpha_rare_rare_Chao1_MS_nowater_set1)


#since the pValue from Shapiro test is close to the alpha level, use also the two-way anova to asses significant differences between samples but simply as check
Chao1_stats_set1 <- aov(JH06_alpha_rare_rare_Chao1_MS_nowater_set1$Chao1 ~ Microhabitat* Species, data = JH06_alpha_rare_rare_Chao1_MS_nowater_set1)
summary(Chao1_stats_set1)


##############
#Shannon
#############

JH06_alpha_rare_rare_Shannon <- as.data.frame(JH06_alpha_rare[ ,4])
rownames(JH06_alpha_rare_rare_Shannon) <- rownames(JH06_alpha_rare)
colnames(JH06_alpha_rare_rare_Shannon) <- c("Shannon")

#Combine the dataset sample description and Shannon OTUs
JH06_alpha_rare_rare_Shannon_MS <- cbind(design_MS, JH06_alpha_rare_rare_Shannon)
JH06_alpha_rare_rare_Shannon_MS

#Order the levels according to a defined order
JH06_alpha_rare_rare_Shannon_MS$Microhabitat <- ordered(JH06_alpha_rare_rare_Shannon_MS$Microhabitat, levels=c("water","waterhizo","biofilm")) 
JH06_alpha_rare_rare_Shannon_MS$Species <- ordered(JH06_alpha_rare_rare_Shannon_MS$Species, levels=c("P","T")) 

##########################################################################
#Figure Shannon OTUs
with(JH06_alpha_rare_rare_Shannon_MS, boxplot(Shannon ~ Microhabitat * Species, xlab = "Samples", ylab = "Shannon OTUs",   main = "OTU evenness", ylim = c(3, 7)))
##########################################################################

#Evenness set 1 stripchart
stripchart(Shannon ~ Microhabitat * Species,
           data=JH06_alpha_rare_rare_Shannon_MS,
           main="OTUs evenness",
           xlab="Samples",
           ylab="Shannon OTUs",
           col="purple",
           group.names=c("water.P","waterhizo.P","biofilm.P","water.T","waterhizo.T","biofilm.T"),
           vertical=TRUE,
           pch=16)

#statistical analysis conducted on biofilm and waterrhizo only
#remove water samples from the analysis
nowater_Shannon_set1 <- rownames(JH06_alpha_rare_rare_Shannon_MS)[which(JH06_alpha_rare_rare_Shannon_MS$Microhabitat!= "water")]
JH06_alpha_rare_rare_Shannon_MS_nowater_set1 <- JH06_alpha_rare_rare_Shannon_MS[nowater_Shannon_set1, ]

#test the normality of the Shannon data
shapiro.test(JH06_alpha_rare_rare_Shannon_MS_nowater_set1$Shannon)

#use a non parametric test to investigate the microhabitat effect
wilcox.test(Shannon ~ Microhabitat, data = JH06_alpha_rare_rare_Shannon_MS_nowater_set1)


########################
#betadiversity calculation
#######################

#PCoA bray distance
JH06_5_bray <- ordinate(JH06_data_phyloseq_5, "PCoA", "bray")
plot_ordination(JH06_data_phyloseq_5, JH06_5_bray, color = "Microhabitat")

#assign shapes to Microhabitat and color to Species
p=plot_ordination(JH06_data_phyloseq_5, JH06_5_bray , shape ="Microhabitat", color = "Species")
p = p + geom_point(size = 8, alpha = 0.75)
p = p + scale_colour_manual(values = c("red", "blue"))
p + ggtitle("PCoA 16S data, Bray distance")


#assign shapes to Species and color to Microhabitat
p=plot_ordination(JH06_data_phyloseq_5, JH06_5_bray , shape ="Species", color = "Microhabitat")
p = p + geom_point(size = 8, alpha = 0.75)
p = p + scale_colour_manual(values = c("red", "blue", "purple"))
p + ggtitle("PCoA 16S data, Bray distance")


#PCoA weighted unifrac distance
JH06_5_unifrac <- ordinate(JH06_data_phyloseq_5, "PCoA", "unifrac", weighted = TRUE)
plot_ordination(JH06_data_phyloseq_5, JH06_5_unifrac, color = "Microhabitat")

#assign shapes to Microhabitat and color to Species
p=plot_ordination(JH06_data_phyloseq_5, JH06_5_unifrac, shape ="Microhabitat", color = "Species")
p = p + geom_point(size = 8, alpha = 0.75)
p = p + scale_colour_manual(values = c("red", "blue"))
p + ggtitle("PCoA 16S data, weighted Unifrac distance")

#assign shapes to Species and color to Microhabitat
p=plot_ordination(JH06_data_phyloseq_5, JH06_5_unifrac, shape= "Species", color = "Microhabitat")
p = p + geom_point(size = 8, alpha = 0.75)
p = p + scale_colour_manual(values = c("red", "blue", "purple"))
p + ggtitle("PCoA 16S data, weighted Unifrac distance")

#assess the significance of these differences
BC <- phyloseq::distance(JH06_data_phyloseq_5, "bray")
WU <- phyloseq::distance(JH06_data_phyloseq_5, "unifrac", weighted= TRUE)
#statistical test
adonis(BC ~ Microhabitat * Species, data= design[row.names(sample_data(JH06_data_phyloseq_5)), ], permutations = 5000)
adonis(WU ~ Microhabitat * Species, data= design[row.names(sample_data(JH06_data_phyloseq_5)), ], permutations = 5000)

#Statistical analysis using Deseq2
#extract count data 
JH06_OTU_counts_integer <- otu_table(JH06_data_phyloseq_5)
countData = as.data.frame(JH06_OTU_counts_integer)
colnames(JH06_OTU_counts_integer)

#the design file containing sample information
colData = design[colnames(JH06_OTU_counts_integer), ]
rownames(colData)

#construct a DESeq dataset combining count data and sample information
JH06_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Description)
JH06_cds

#execute the differential count analysis with the function DESeq 
JH06_cds_test <- DESeq(JH06_cds, fitType="local", betaPrior = FALSE) 

#define the OTUs significantly enriched in the biofilm samples
Biofilm_NB_Phragmites <- results(JH06_cds_test, contrast = c("Description",  "biofilmP", "waterhizoP")) 
Biofilm_NB_Typha <- results(JH06_cds_test, contrast = c("Description",  "biofilmT", "waterhizoT")) 

#inspect a result file
Biofilm_NB_Phragmites  
mcols(Biofilm_NB_Phragmites  , use.names=TRUE)
Biofilm_NB_Typha  
mcols(Biofilm_NB_Typha  , use.names=TRUE)

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
Biofilm_NB_Phragmites_FDR_005 <- Biofilm_NB_Phragmites[(rownames(Biofilm_NB_Phragmites)[which(Biofilm_NB_Phragmites$padj <0.05)]), ]
Biofilm_NB_Typha_FDR_005 <- Biofilm_NB_Typha[(rownames(Biofilm_NB_Typha)[which(Biofilm_NB_Typha$padj <0.05)]), ]

#enriched in biofilm (positive fold change)
Biofilm_NB_Phragmites_enriched <-  Biofilm_NB_Phragmites[(rownames(Biofilm_NB_Phragmites)[which(Biofilm_NB_Phragmites$log2FoldChange > 0)]), ]
Biofilm_NB_Typha_enriched <-  Biofilm_NB_Typha[(rownames(Biofilm_NB_Typha)[which(Biofilm_NB_Typha$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
Biofilm_NB_Phragmites_enriched_FDR005 <- intersect(rownames(Biofilm_NB_Phragmites_FDR_005), rownames(Biofilm_NB_Phragmites_enriched))
Biofilm_NB_Typha_enriched_FDR005 <- intersect(rownames(Biofilm_NB_Typha_FDR_005), rownames(Biofilm_NB_Typha_enriched))

#Define the number of OTUs significantly enriched in the biofilm of each plant species
length(Biofilm_NB_Phragmites_enriched_FDR005)
length(Biofilm_NB_Typha_enriched_FDR005)

#visualise the relationships among enriched OTUs
#define the areas of the diagram
length(Biofilm_NB_Phragmites_enriched_FDR005)
length(Biofilm_NB_Typha_enriched_FDR005)
length(intersect(Biofilm_NB_Phragmites_enriched_FDR005, Biofilm_NB_Typha_enriched_FDR005))       

#draw the diagram
#venn diagram
dev.off()
draw.pairwise.venn(119, 92, 36, category = c("P. australis", "T. laifolia"), lty = rep("blank", 2), fill = c("red2", "royalblue3"), scaled = TRUE)

#identify the OTUs conserved across the two plants
intersect_set1 <- dat_tax_noPlants[intersect(Biofilm_NB_Phragmites_enriched_FDR005, Biofilm_NB_Typha_enriched_FDR005), ]
intersect_set1 

#save the file for the reproducibility of the code
write.table(intersect_set1, file="intersect_set1.txt", sep="\t")

#create dataframes containing both statistical and taxonomic information of the relevant OTUs
Biofilm_NB_Phragmites_enriched_taxa <- JH06_taxa_ordered[Biofilm_NB_Phragmites_enriched_FDR005, ]
Biofilm_NB_Phragmites_enriched_FDR005_taxa <- cbind(as.data.frame(Biofilm_NB_Phragmites[Biofilm_NB_Phragmites_enriched_FDR005, ]), Biofilm_NB_Phragmites_enriched_taxa)
Biofilm_NB_Typha_enriched_taxa <- JH06_taxa_ordered[Biofilm_NB_Typha_enriched_FDR005, ]
Biofilm_NB_Typha_enriched_FDR005_taxa <- cbind(as.data.frame(Biofilm_NB_Typha[Biofilm_NB_Typha_enriched_FDR005, ]), Biofilm_NB_Typha_enriched_taxa)

#save to be open in excel
write.table(Biofilm_NB_Phragmites_enriched_FDR005_taxa, file="Biofilm_NB_Phragmites_enriched_FDR005_taxa.txt ", sep="\t")
write.table(Biofilm_NB_Typha_enriched_FDR005_taxa, file="Biofilm_NB_Typha_enriched_FDR005_taxa.txt ", sep="\t")


#visualise the 32 OTUs enriched in both biofilm
#create a phyloseq object with the 36 OTUs
JH06_data_phyloseq_3_samples_biofilm <- prune_taxa(intersect(Biofilm_NB_Phragmites_enriched_FDR005, Biofilm_NB_Typha_enriched_FDR005), JH06_data_phyloseq_3_samples_no_contaminant)   
JH06_data_phyloseq_3_samples_biofilm_2 <- subset_samples(JH06_data_phyloseq_3_samples_biofilm, Microhabitat == "biofilm")
JH06_data_phyloseq_3_samples_biofilm_2




###############################################################################
##################   samples_set2 ############################

#design file
design <- read.delim("Map_JH06_detailed_phyloseq_L.txt", sep = "\t", header=TRUE, row.names=1)
design

#OK we should remove the technical replicates: we could use them to define the minimum number of reads to be retained for an OTU to be considered reproducible
#firstly let's see the results with the replicates 2 flagged in the map file as replicate 
JH06_data_phyloseq_3_samples_no_contaminant_T2 <- subset_samples(JH06_data_phyloseq_3_samples_no_contaminant, replicate ==1)

#how many reads have we retained?
sort(sample_sums(JH06_data_phyloseq_3_samples_no_contaminant_T2))

#Remove OTUs not seen more than 25 times in at least 20% of the samples
JH06_data_phyloseq_4_1 = filter_taxa(JH06_data_phyloseq_3_samples_no_contaminant_T2, function(x) sum(x > 25) > (0.2*length(x)), TRUE)
JH06_data_phyloseq_4_1 

#how many reads are left?
sample_sums(JH06_data_phyloseq_4_1)
sort(sample_sums(JH06_data_phyloseq_4_1))

#let's rarefy the dataset at 66,000 reads
JH06_data_phyloseq_66K_set2 <- rarefy_even_depth(JH06_data_phyloseq_4_1 , sample.size = 66000)

#extract the data
JH06_data_phyloseq_66K_set2_table <- as.data.frame(otu_table(JH06_data_phyloseq_66K_set2))
dim(JH06_data_phyloseq_66K_set2_table)

#save the file for the reproducibility of the code
write.table(JH06_data_phyloseq_66K_set2_table, file="JH06_data_phyloseq_66K_set2_table.txt", sep="\t")

#import the rarefied OTU counts (note file name counts2.txt this file has been generated in excel and includes the #OTU ID as header of the first column)
dat_count_JH06_66K_set2 <- read.delim("JH06_data_phyloseq_66K_set2_table.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
dim(dat_count_JH06_66K_set2)

#Note: the above table will be used to generate Supplementary Figure 1-set 2 in the revised manuscript
#generate a PICRUSt table containing all OTUs
#prune OTU table below the threshold from the taxonomy information
dat_tax_noPlants_66K_2 <- as.data.frame(dat_tax_noPlants[rownames(dat_count_JH06_66K_set2), ])
colnames(dat_tax_noPlants_66K_2) <- c("ConsensusLineage") 
colnames(dat_tax_noPlants_66K_2)
#merge count and taxonomy files
dat_info_rare_66K_2 <- cbind(dat_count_JH06_66K_set2, dat_tax_noPlants_66K_2)
dim(dat_info_rare_66K_2)
colnames(dat_info_rare_66K_2)

#generate a design file set 1 to generate taxa chart including set 2 samples only
design_set2 <- design[colnames(dat_count_JH06_66K_set2), ]
design_set2_qiime <- as.data.frame(design_set1[, 1:4])
rownames(design_set2_qiime)
colnames(design_set2_qiime)
dim(design_set2_qiime)

#save the files for the reproducibility of the code
write.table(dat_info_rare_66K, file="JH06_info_rare_66K_set2.txt", sep="\t")
#Note: for taxonomy calculation in QIIME the file JH06_info_rare_66K_set2_2.txt, with the header of the first column '#OTU ID' was generated in excel
write.table(design_set2_qiime, file="JH06_map_set2_qiime.txt", sep="\t")
#Note: for taxonomy calculation in QIIME the file JH06_map_set2_qiime_2.txt, with the header of the first column '#SampleID' was generated in excel

#Import average value phylum level
dat_info_taxa_Phylum_set_2 <- read.delim("Description_otu_table_L2_set2.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
colnames(dat_info_taxa_Phylum_set_2)
dim(dat_info_taxa_Phylum_set_2)

# transform the data in % for visualisation
dat_norm_Phylum_set_2 <- dat_info_taxa_Phylum_set_2 * 100 
dat_norm_Phylum_set_2[1:5, ]
colSums(dat_norm_Phylum_set_2)

# determine the average % of reads for each phylum
Phylum_mean_sorted_set_2 <- dat_norm_Phylum_set_2[(order(-rowSums(dat_norm_Phylum_set_2))), ] 

#Identifythe phyla whose average abundance in above 1% and their aggregated relative abundance
Phylum_mean_topRank_set_2 <- Phylum_mean_sorted_set_2[rownames(Phylum_mean_sorted_set_2)[which(rowMeans(Phylum_mean_sorted_set_2) > 1)], ]
dim(Phylum_mean_topRank_set_2)
colSums(Phylum_mean_topRank_set_2)
rowMeans(Phylum_mean_topRank_set_2)

#transform in matrix for plotting purposes
Phylum_mean_topRank_set_2 <- as.matrix(Phylum_mean_topRank_set_2) 
colnames(Phylum_mean_topRank_set_2) 

#re-arrange the order for plotting purposes
col.order <- c("waterP","waterhizoP","biofilmP","waterT","waterhizoT","biofilmT")
Phylum_mean_topRank_set_2 <- Phylum_mean_topRank_set_2[,col.order]

# Stacked Bar Plot with Colors and Legend
barplot(Phylum_mean_topRank_set_2, main="Phylum Distribution",
        xlab="Samples", ylab = "Reads %", ylim = c(0,100), col=c("darkblue","red", "green", "cyan", "orange", "brown", "black", "magenta", "grey"), beside=FALSE,   legend = rownames(Phylum_mean_topRank_set_2))

#Due to size limits the legend covers part of the graph, save the graph as .eps file and in illustrator uou can adjust this (and other)  graphical issues
#barplot_no legend: you can save this as image and then reconstruct the leged using the previous figure. Note: the order of samples can be inferred from the command on line 256
barplot(Phylum_mean_topRank_set_2, main="Phylum Distribution",
        xlab="Samples", ylab = "Reads %", ylim = c(0,100), col=c("darkblue","red", "green", "cyan", "orange", "brown", "black", "magenta", "grey"), beside=FALSE)


#convert the file into a phyloseq-compatible OTU table
JH06_OTU_6 <- otu_table(dat_count_JH06_66K_set2, taxa_are_rows=TRUE)

#create a new phyloseq object
#Note: the number of samples/OTUs in the otu table written drive the composition of the phyloseq object
#For this reason we could use older mapping, taxonomy and phylogenetic tree: they will be "pruned" of the OTUs/samples not present in the OTU table
JH06_data_phyloseq_6 <- merge_phyloseq(JH06_OTU_6, JH06_taxa, JH06_map,  JH06_tree)

#inspect the created object
JH06_data_phyloseq_6

################################################################################################
#OK we will use the JH06_data_phyloseq_6 for the analysis and figure generation in the manuscript
###############################################################################################

#first thing to check is alpha diversity
#Index calculations
JH06_alpha_rare1 <-  estimate_richness(JH06_data_phyloseq_6, measures = c("Observed", "Chao1", "Shannon"))

#generate a new data frame combining alphadiversity calculation and sample information for data visualisation
#Sample information

#Microhabitat
design_microhabitat1 <- as.data.frame(sample_data(JH06_data_phyloseq_6)[, 1])
rownames(design_microhabitat1) <- rownames(sample_data(JH06_data_phyloseq_6))
colnames(design_microhabitat1) <- c("Microhabitat")

#Species
design_Species1 <- as.data.frame(sample_data(JH06_data_phyloseq_6)[, 2])
rownames(design_Species1) <- rownames(sample_data(JH06_data_phyloseq_6))
colnames(design_Species1) <- c("Species")

#combine the file
design_MS1 <- cbind(design_microhabitat1, design_Species1)

##############
#Observed OTUs
#############

JH06_alpha_rare_rare_Observed1 <- as.data.frame(JH06_alpha_rare1[ ,1])
rownames(JH06_alpha_rare_rare_Observed1) <- rownames(JH06_alpha_rare1)
colnames(JH06_alpha_rare_rare_Observed1) <- c("Observed")

#Combine the dataset sample description and Observed OTUs
JH06_alpha_rare_rare_Observed_MS1 <- cbind(design_MS1, JH06_alpha_rare_rare_Observed1)
JH06_alpha_rare_rare_Observed_MS1

#Order the levels according to a defined order
JH06_alpha_rare_rare_Observed_MS1$Microhabitat <- ordered(JH06_alpha_rare_rare_Observed_MS1$Microhabitat, levels=c("water","waterhizo","biofilm")) 
JH06_alpha_rare_rare_Observed_MS1$Species <- ordered(JH06_alpha_rare_rare_Observed_MS1$Species, levels=c("P","T")) 

########################################################################
#Figure Observed OTUs
with(JH06_alpha_rare_rare_Observed_MS1, boxplot(Observed ~ Microhabitat * Species, xlab = "Samples", ylab = "Observed OTUs",   main = "OTU richness", ylim = c(800, 1800)))
#########################################################################


#Observed OTUs set2 stripchart
stripchart(Observed ~ Microhabitat * Species,
           data=JH06_alpha_rare_rare_Observed_MS1,
           main="OTUs richness",
           xlab="Samples",
           ylab="Observed OTUs",
           col="blue", 
           group.names=c("water.P","waterhizo.P","biofilm.P","water.T","waterhizo.T","biofilm.T"),
           vertical=TRUE,
           pch=16)



#statistical analysis conducted on biofilm and waterrhizo only
#remove water samples from the analysis
nowater_set2 <- rownames(JH06_alpha_rare_rare_Observed_MS1)[which(JH06_alpha_rare_rare_Observed_MS1$Microhabitat!= "water")]
JH06_alpha_rare_rare_Observed_MS1_nowater_set2 <- JH06_alpha_rare_rare_Observed_MS1[nowater_set2, ]

#test the normality of the observed data
shapiro.test(JH06_alpha_rare_rare_Observed_MS1_nowater_set2$Observed)

# data not normally distributed (p-value < 0.05); use a non parametric test to investigate the microhabitat effect
wilcox.test(Observed ~ Microhabitat, data = JH06_alpha_rare_rare_Observed_MS1_nowater_set2)



##############
#Chao 1
#############

JH06_alpha_rare_rare_Chao1_1 <- as.data.frame(JH06_alpha_rare1[ ,2])
rownames(JH06_alpha_rare_rare_Chao1_1) <- rownames(JH06_alpha_rare1)
colnames(JH06_alpha_rare_rare_Chao1_1) <- c("Chao1")

#Combine the dataset sample description and Chao1 OTUs
JH06_alpha_rare_rare_Chao1_MS1 <- cbind(design_MS1, JH06_alpha_rare_rare_Chao1_1)
JH06_alpha_rare_rare_Chao1_MS1

#Order the levels according to a defined order
JH06_alpha_rare_rare_Chao1_MS1$Microhabitat <- ordered(JH06_alpha_rare_rare_Chao1_MS1$Microhabitat, levels=c("water","waterhizo","biofilm")) 
JH06_alpha_rare_rare_Chao1_MS1$Species <- ordered(JH06_alpha_rare_rare_Chao1_MS1$Species, levels=c("P","T")) 


###################################################
#Figure Chao1 OTUs
with(JH06_alpha_rare_rare_Chao1_MS1, boxplot(Chao1 ~ Microhabitat * Species, xlab = "Samples", ylab = "Chao1 OTUs",   main = "OTU richness", ylim = c(900, 1900)))
###################################################



#Chao1 set2 stripchart
stripchart(Chao1 ~ Microhabitat * Species,
           data=JH06_alpha_rare_rare_Chao1_MS1,
           main="OTUs richness",
           xlab="Samples",
           ylab="Chao1 OTUs",
           col="darkgreen",
           group.names=c("water.P","waterhizo.P","biofilm.P","water.T","waterhizo.T","biofilm.T"),
           vertical=TRUE,
           pch=16)


#statistical analysis conducted on biofilm and waterrhizo only
#remove water samples from the analysis
nowater_set2 <- rownames(JH06_alpha_rare_rare_Chao1_MS1)[which(JH06_alpha_rare_rare_Chao1_MS1$Microhabitat!= "water")]
JH06_alpha_rare_rare_Chao1_MS1_nowater_set2 <- JH06_alpha_rare_rare_Chao1_MS1[nowater_set2, ]

#test the normality of the observed data
shapiro.test(JH06_alpha_rare_rare_Chao1_MS1_nowater_set2$Chao1)

# data not normally distributed (p-value < 0.05); use a non parametric test to investigate the microhabitat effect
wilcox.test(Chao1 ~ Microhabitat, data = JH06_alpha_rare_rare_Chao1_MS1_nowater_set2)


##############
#Shannon
#############

JH06_alpha_rare_rare_Shannon1 <- as.data.frame(JH06_alpha_rare1[ ,4])
rownames(JH06_alpha_rare_rare_Shannon1) <- rownames(JH06_alpha_rare1)
colnames(JH06_alpha_rare_rare_Shannon1) <- c("Shannon")

#Combine the dataset sample description and Shannon OTUs
JH06_alpha_rare_rare_Shannon_MS1 <- cbind(design_MS1, JH06_alpha_rare_rare_Shannon1)
JH06_alpha_rare_rare_Shannon_MS1

#Order the levels according to a defined order
JH06_alpha_rare_rare_Shannon_MS1$Microhabitat <- ordered(JH06_alpha_rare_rare_Shannon_MS1$Microhabitat, levels=c("water","waterhizo","biofilm")) 
JH06_alpha_rare_rare_Shannon_MS1$Species <- ordered(JH06_alpha_rare_rare_Shannon_MS1$Species, levels=c("P","T")) 

###############################################################
#Figure Shannon OTUs
with(JH06_alpha_rare_rare_Shannon_MS1, boxplot(Shannon ~ Microhabitat * Species, xlab = "Samples", ylab = "Shannon OTUs",   main = "OTU evenness", ylim = c(3, 7)))
###############################################################


#Evenness set 2 stripchart
stripchart(Shannon ~ Microhabitat * Species,
           data=JH06_alpha_rare_rare_Shannon_MS1,
           main="OTUs evenness",
           xlab="Samples",
           ylab="Shannon OTUs",
           col="purple",
           group.names=c("water.P","waterhizo.P","biofilm.P","water.T","waterhizo.T","biofilm.T"),
           vertical=TRUE,
           pch=16)


#statistical analysis conducted on biofilm and waterrhizo only
#remove water samples from the analysis
nowater_set2 <- rownames(JH06_alpha_rare_rare_Shannon_MS1)[which(JH06_alpha_rare_rare_Shannon_MS1$Microhabitat!= "water")]
JH06_alpha_rare_rare_Shannon_MS1_nowater_set2 <- JH06_alpha_rare_rare_Shannon_MS1[nowater_set2, ]

#test the normality of the observed data
shapiro.test(JH06_alpha_rare_rare_Shannon_MS1_nowater_set2$Shannon)

#use a non parametric test to investigate the microhabitat effect
wilcox.test(Shannon ~ Microhabitat, data = JH06_alpha_rare_rare_Shannon_MS1_nowater_set2)


########################
#betadiversity calculation
#######################

#PCoA weighted bray distance
JH06_6_bray <- ordinate(JH06_data_phyloseq_6, "PCoA", "bray")
plot_ordination(JH06_data_phyloseq_6, JH06_6_bray, color = "Microhabitat")

#assign shapes to Microhabitat and color to Species
p=plot_ordination(JH06_data_phyloseq_6, JH06_6_bray , shape ="Microhabitat", color = "Species")
p = p + geom_point(size = 8, alpha = 0.75)
p = p + scale_colour_manual(values = c("red", "blue"))
p + ggtitle("PCoA 16S data, Bray distance")

#assign shapes to Species and color to Microhabitat
p=plot_ordination(JH06_data_phyloseq_6, JH06_6_bray , shape ="Species", color = "Microhabitat")
p = p + geom_point(size = 8, alpha = 0.75)
p = p + scale_colour_manual(values = c("red", "blue", "purple"))
p + ggtitle("PCoA 16S data, Bray distance")

#PCoA weighted unifrac distance
JH06_6_unifrac <- ordinate(JH06_data_phyloseq_6, "PCoA", "unifrac", weighted = TRUE)
plot_ordination(JH06_data_phyloseq_6, JH06_6_unifrac, color = "Microhabitat")

#assign shapes to Microhabitat and color to Species
p=plot_ordination(JH06_data_phyloseq_6, JH06_6_unifrac, shape ="Microhabitat", color = "Species")
p = p + geom_point(size = 8, alpha = 0.75)
p = p + scale_colour_manual(values = c("red", "blue"))
p + ggtitle("PCoA 16S data, weighted Unifrac distance")

#assign shapes to Species and color to Microhabitat
p=plot_ordination(JH06_data_phyloseq_6, JH06_6_unifrac, shape= "Species", color = "Microhabitat")
p = p + geom_point(size = 8, alpha = 0.75)
p = p + scale_colour_manual(values = c("red", "blue", "purple"))
p + ggtitle("PCoA 16S data, weighted Unifrac distance")

#assess the significance of these differences
BC1 <- phyloseq::distance(JH06_data_phyloseq_6, "bray")
WU1 <- phyloseq::distance(JH06_data_phyloseq_6, "unifrac", weighted= TRUE)
#statistical test
adonis(BC1 ~ Microhabitat * Species, data= design[row.names(sample_data(JH06_data_phyloseq_6)), ], permutations = 5000)
adonis(WU1 ~ Microhabitat * Species, data= design[row.names(sample_data(JH06_data_phyloseq_6)), ], permutations = 5000)

#Statistical analysis using Deseq2
#extract count data 
JH06_OTU_counts_integer1 <- otu_table(JH06_data_phyloseq_6)
countData1 = as.data.frame(JH06_OTU_counts_integer1)
colnames(JH06_OTU_counts_integer1)

#the design file containing sample information
colData1 = design[colnames(JH06_OTU_counts_integer1), ]
rownames(colData1)

#construct a DESeq dataset combining count data and sample information
JH06_cds1 <- DESeqDataSetFromMatrix(countData =countData1, colData=colData1 , design= ~Description)
JH06_cds1

#execute the differential count analysis with the function DESeq 
JH06_cds_test1 <- DESeq(JH06_cds1, fitType="local", betaPrior = FALSE) 

#define the OTUs significantly enriched in the biofilm samples
Biofilm_NB_Phragmites1 <- results(JH06_cds_test1, contrast = c("Description",  "biofilmP", "waterhizoP")) 
Biofilm_NB_Typha1 <- results(JH06_cds_test1, contrast = c("Description",  "biofilmT", "waterhizoT")) 

#inspect a result file
Biofilm_NB_Phragmites1  
mcols(Biofilm_NB_Phragmites1  , use.names=TRUE)
Biofilm_NB_Typha1  
mcols(Biofilm_NB_Typha1  , use.names=TRUE)

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
Biofilm_NB_Phragmites_FDR_005_1 <- Biofilm_NB_Phragmites1[(rownames(Biofilm_NB_Phragmites1)[which(Biofilm_NB_Phragmites1$padj <0.05)]), ]
Biofilm_NB_Typha_FDR_005_1 <- Biofilm_NB_Typha1[(rownames(Biofilm_NB_Typha1)[which(Biofilm_NB_Typha1$padj <0.05)]), ]

#enriched in biofilm (positive fold change)
Biofilm_NB_Phragmites_enriched1 <-  Biofilm_NB_Phragmites1[(rownames(Biofilm_NB_Phragmites1)[which(Biofilm_NB_Phragmites1$log2FoldChange > 0)]), ]
Biofilm_NB_Typha_enriched1 <-  Biofilm_NB_Typha1[(rownames(Biofilm_NB_Typha1)[which(Biofilm_NB_Typha1$log2FoldChange > 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
Biofilm_NB_Phragmites_enriched_FDR005_1 <- intersect(rownames(Biofilm_NB_Phragmites_FDR_005_1), rownames(Biofilm_NB_Phragmites_enriched1))
Biofilm_NB_Typha_enriched_FDR005_1 <- intersect(rownames(Biofilm_NB_Typha_FDR_005_1), rownames(Biofilm_NB_Typha_enriched1))

#Define the number of OTUs significantly enriched in the biofilm of each plant species
length(Biofilm_NB_Phragmites_enriched_FDR005_1)
length(Biofilm_NB_Typha_enriched_FDR005_1)

#visualise the relationships among enriched OTUs
#define the areas of the diagram
length(Biofilm_NB_Phragmites_enriched_FDR005_1)
length(Biofilm_NB_Typha_enriched_FDR005_1)
length(intersect(Biofilm_NB_Phragmites_enriched_FDR005_1, Biofilm_NB_Typha_enriched_FDR005_1))       

#draw the diagram
#venn diagram
dev.off()
draw.pairwise.venn(120, 131, 39, category = c("P. australis", "T. latifolia"), lty = rep("blank", 2), fill = c("red2", "royalblue3"), scaled = TRUE)

#identify the OTUs conserved across the two plants
intersect_set2 <- dat_tax_noPlants[intersect(Biofilm_NB_Phragmites_enriched_FDR005_1, Biofilm_NB_Typha_enriched_FDR005_1), ]

#save the file for the reproducibility of the code
write.table(intersect_set2, file="intersect_set2.txt", sep="\t")

#create dataframes containing both statistical and taxonomic information of the relevant OTUs
Biofilm_NB_Phragmites_enriched_taxa1 <- JH06_taxa_ordered[Biofilm_NB_Phragmites_enriched_FDR005_1, ]
Biofilm_NB_Phragmites_enriched_FDR005_taxa1 <- cbind(as.data.frame(Biofilm_NB_Phragmites1[Biofilm_NB_Phragmites_enriched_FDR005_1, ]), Biofilm_NB_Phragmites_enriched_taxa1)
Biofilm_NB_Typha_enriched_taxa1 <- JH06_taxa_ordered[Biofilm_NB_Typha_enriched_FDR005_1, ]
Biofilm_NB_Typha_enriched_FDR005_taxa1 <- cbind(as.data.frame(Biofilm_NB_Typha1[Biofilm_NB_Typha_enriched_FDR005_1, ]), Biofilm_NB_Typha_enriched_taxa1)

#save to be open in excel
write.table(Biofilm_NB_Phragmites_enriched_FDR005_taxa1, file="Biofilm_NB_Phragmites_enriched_FDR005_taxa1.txt ", sep="\t")
write.table(Biofilm_NB_Typha_enriched_FDR005_taxa1, file="Biofilm_NB_Typha_enriched_FDR005_taxa1.txt ", sep="\t")


#To identify OTUs significantly enriched in both replicates sets for each plant
Biofilm_NB_Phragmites_enriched_set1_set2 <- intersect(rownames(Biofilm_NB_Phragmites_enriched_FDR005_taxa), rownames(Biofilm_NB_Phragmites_enriched_FDR005_taxa1))
Biofilm_NB_Typha_enriched_set1_set2 <- intersect(rownames(Biofilm_NB_Typha_enriched_FDR005_taxa), rownames(Biofilm_NB_Typha_enriched_FDR005_taxa1))


#identify the OTUs conserved across the two plants
enriched_Phragmites_intersect_set1_set2 <- dat_tax_noPlants[(Biofilm_NB_Phragmites_enriched_set1_set2), ]
enriched_Typha_intersect_set1_set2 <- dat_tax_noPlants[(Biofilm_NB_Typha_enriched_set1_set2), ]

#save to be open in excel
write.table(enriched_Phragmites_intersect_set1_set2, file="enriched_Phragmites_intersect_set1_set2.txt ", sep="\t")
write.table(enriched_Typha_intersect_set1_set2, file="enriched_Typha_intersect_set1_set2.txt ", sep="\t")


#identify the OTUs conserved across the two plants
intersect_enriched_set1_set2 <- dat_tax_noPlants[intersect(Biofilm_NB_Phragmites_enriched_set1_set2, Biofilm_NB_Typha_enriched_set1_set2), ]

#save to be open in excel
write.table(intersect_enriched_set1_set2, file="intersect_enriched_set1_set2.txt ", sep="\t")


#visualise the relationships among enriched OTUs
#define the areas of the diagram
length(Biofilm_NB_Phragmites_enriched_set1_set2)
length(Biofilm_NB_Typha_enriched_set1_set2)
length(intersect(Biofilm_NB_Phragmites_enriched_set1_set2, Biofilm_NB_Typha_enriched_set1_set2))       


#draw the diagram
#venn diagram
dev.off()
draw.pairwise.venn(78, 66, 23, category = c("P. australis", "T. latifolia"), lty = rep("blank", 2), fill = c("red2", "royalblue3"), scaled = TRUE)

