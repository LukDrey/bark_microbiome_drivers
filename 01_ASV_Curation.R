#################################################################
##                          Section 1                          ##
##                       Package Loading                       ##
#################################################################

library(here); packageVersion("here")
# 1.0.1
library(decontam); packageVersion("decontam")
#1.16.0
library(phyloseq); packageVersion("phyloseq")
#1.40.0
library(lulu); packageVersion("lulu")
#0.1.0
library(Biostrings); packageVersion("Biostrings")
#2.64.0
library(tidyverse); packageVersion("tidyverse")
#1.3.2

##----------------------------------------------------------------
##                        Custom Functions                       -
##----------------------------------------------------------------

# Function to create a FASTA file from the curated algae asvs. 
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"seq_name_algae"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"sequence_algae"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

##################################################################
##                          Section 2                           ##
##                Decontam Run and LULU curation                ##
##################################################################

##---------
##  Algae  
##---------

###
# Decontam
###

####
# Data Loading
####

# Load ASV table for algae (available as supplementary data).
algae_asv <- utils::read.csv(here::here("Data", 'asv_table_algae.txt'), header = T, sep = '\t')

# Load FASTA file to get ASV_IDs.
algae_seqs_fasta <- Biostrings::readDNAStringSet(here::here("Data", 'ASVs_algae.fa'))

# Make a dataframe of the sequences and their ASV ID. 
seq_name_algae <- base::names(algae_seqs_fasta)
sequence_algae <- base::paste(algae_seqs_fasta)
algae_rep_seqs <- base::data.frame(seq_name_algae, sequence_algae)

# Join the ASV table and the representative sequences.
algae_asv_IDs <- dplyr::left_join(algae_rep_seqs, algae_asv, by = 'sequence_algae')

# Naming the Samples. 
algae_asv_IDs <- algae_asv_IDs %>% 
  dplyr::rename_with(~paste0("Sample_", .), -c(1:2))

# Remove the dot from PCR negative controls (PCN) and multiplex controls (MPC). 
base::colnames(algae_asv_IDs) <- sub("[.]", "_", x = base::colnames(algae_asv_IDs))

# set ASV ID as the rownames.
base::rownames(algae_asv_IDs) <- algae_asv_IDs$seq_name_algae
algae_asv_IDs$seq_name_algae <- NULL

# Load the DNA concentration necessary for decontam (available as supplementary data). 
algae_conc <- utils::read.csv(here("Data", 'algae_decontam_conc_big.csv'), header = T, sep =',')
base::rownames(algae_conc) <- paste0("Sample_", algae_conc$sample_ID)
base::rownames(algae_conc) <- sub("[-]", "_", x = base::rownames(algae_conc))

# Make the parts of the phyloseq object.
ASV_mat_decontam <- base::data.matrix(algae_asv_IDs)
ASV_algae_decontam <- phyloseq::otu_table(ASV_mat_decontam, taxa_are_rows = TRUE)
sampledata_algae_decontam <- phyloseq::sample_data(algae_conc)

# Combine with phyloseq
ps_algae_decontam <- phyloseq::phyloseq(ASV_algae_decontam, sampledata_algae_decontam)
ps_algae_decontam

# Check the library sizes of the Samples and Controls.
df_algae <- base::as.data.frame(sample_data(ps_algae_decontam)) # Put sample_data into a ggplot-friendly data.frame
df_algae$LibrarySize <- phyloseq::sample_sums(ps_algae_decontam)
df_algae <- df_algae[base::order(df_algae$LibrarySize),]
df_algae$Index <- base::seq(base::nrow(df_algae))
ggplot2::ggplot(data=df_algae, ggplot2::aes(x=Index, y=LibrarySize, color=Sample_or_Control)) +
  ggplot2::geom_point()

# Decontam with combined frequency and prevalence approach. 
phyloseq::sample_data(ps_algae_decontam)$is.neg <- phyloseq::sample_data(ps_algae_decontam)$Sample_or_Control == "Control Sample"
contam_algae_combi <- decontam::isContaminant(ps_algae_decontam,
                                    method = 'combined',
                                    neg = 'is.neg',
                                    conc = 'quant_reading')
base::table(contam_algae_combi$contaminant)
base::which(contam_algae_combi$contaminant)

# Trim the identified contaminants from the phyloseq object.
# In this case there are 15 identified contaminants. 
ps_algae_noncontam <- phyloseq::prune_taxa(!contam_algae_combi$contaminant,
                                 ps_algae_decontam)
ps_algae_noncontam

###
#LULU curation
###

# Remove taxa without reads.
ps_algae_noncontam_pruned <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_algae_noncontam) > 0,
                                        ps_algae_noncontam)

# Transform the ASV table from Phyloseq to a data frame. 
ASV_table_algae <- base::as.data.frame(phyloseq::otu_table(ps_algae_noncontam_pruned))

# Load the matchlist created on the server. 
algae_matchlist <- utils::read.table(here::here("Data", 'match_list_algae.txt'),
                              header = F,
                              as.is = T,
                              stringsAsFactors = F)

# Run the LULU algorithm. 
ASV_table_algae_cur <- lulu::lulu(ASV_table_algae, algae_matchlist)

ASV_table_algae_cur$curated_count
ASV_table_algae_cur$discarded_count

##---------
##  Fungi  
##---------

###
# Decontam
###

# Load ASV table for fungi (available as supplementary data).
fungi_asv <- utils::read.csv(here("Data", 'asv_table_fungi.txt'), header = T, sep = '\t')

# Load FASTA file to get ASV_IDs.
fungi_seqs_fasta <- Biostrings::readDNAStringSet(here::here("Data", 'ASVs_fungi.fa'))

# Make a dataframe of the sequences and their ASV ID. 
seq_name_fungi <- base::names(fungi_seqs_fasta)
sequence_fungi <- base::paste(fungi_seqs_fasta)
fungi_rep_seqs <- base::data.frame(seq_name_fungi, sequence_fungi)

# Join the ASV table and the representative sequences.
fungi_asv_IDs <- dplyr::left_join(fungi_rep_seqs, fungi_asv, by = 'sequence_fungi')

# Naming the Samples. 
fungi_asv_IDs <- fungi_asv_IDs %>% 
  dplyr::rename_with(~base::paste0("Sample_", .), -c(1:2))

# set ASV ID as the rownames.
base::rownames(fungi_asv_IDs) <- fungi_asv_IDs$seq_name_fungi
fungi_asv_IDs$seq_name_fungi <- NULL

# Load the DNA concentration necessary for decontam (available as supplementary data). 
fungi_conc <- utils::read.csv(here("Data", 'fungi_decontam_conc_big.csv'), header = T, sep =',')
base::rownames(fungi_conc) <- base::paste0("Sample_", fungi_conc$sample_ID)
base::rownames(fungi_conc) <- base::sub("[-]", "_", x = base::rownames(fungi_conc))

# Make the parts of the phyloseq object.
ASV_mat_decontam_fun <- base::data.matrix(fungi_asv_IDs)
ASV_fungi_decontam <- phyloseq::otu_table(ASV_mat_decontam_fun, taxa_are_rows = TRUE)
sampledata_fungi_decontam <- phyloseq::sample_data(fungi_conc)

# Combine with phyloseq
ps_fungi_decontam <- phyloseq::phyloseq(ASV_fungi_decontam, sampledata_fungi_decontam)
ps_fungi_decontam

# Check the library sizes of the Samples and Controls.
df_fungi <- base::as.data.frame(phyloseq::sample_data(ps_fungi_decontam)) # Put sample_data into a ggplot-friendly data.frame
df_fungi$LibrarySize <- phyloseq::sample_sums(ps_fungi_decontam)
df_fungi <- df_fungi[base::order(df_fungi$LibrarySize),]
df_fungi$Index <- base::seq(base::nrow(df_fungi))
ggplot2::ggplot(data=df_fungi, ggplot2::aes(x=Index, y=LibrarySize, color=Sample_or_Control)) +
  ggplot2::geom_point()

# Decontam with combined frequency and prevalence approach. 
phyloseq::sample_data(ps_fungi_decontam)$is.neg <- phyloseq::sample_data(ps_fungi_decontam)$Sample_or_Control == "Control Sample"
contam_fungi_combi <- decontam::isContaminant(ps_fungi_decontam,
                                    method = 'combined',
                                    neg = 'is.neg',
                                    conc = 'quant_reading')
base::table(contam_fungi_combi$contaminant)
base::which(contam_fungi_combi$contaminant)

# Trim the identified contaminants from the phyloseq object 
ps_fungi_noncontam <- phyloseq::prune_taxa(!contam_fungi_combi$contaminant,
                                 ps_fungi_decontam)
ps_fungi_noncontam

###
#LULU curation
###

# Remove taxa without reads.
ps_fungi_noncontam_pruned <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_fungi_noncontam) > 0,
                                        ps_fungi_noncontam)
ASV_table_fungi <- base::as.data.frame(phyloseq::otu_table(ps_fungi_noncontam_pruned))
fungi_matchlist <- utils::read.table(here::here("Data", 'match_list_fungi.txt'),
                              header = F,
                              as.is = T,
                              stringsAsFactors = F)

# Run the LULU algorithm. 

ASV_table_fungi_cur <- lulu::lulu(ASV_table_fungi, fungi_matchlist)

ASV_table_fungi_cur$curated_count
ASV_table_fungi_cur$discarded_count

##---------
## Bacteria  
##---------

###
# Decontam
###

# Load ASV table for bacteria (available as supplementary data).
bacteria_asv <- utils::read.csv(here::here("Data", 'asv_table_bacteria.txt'), header = T, sep = '\t')

# Load FASTA file to get ASV_IDs.
bacteria_seqs_fasta <- Biostrings::readDNAStringSet(here::here("Data", 'ASVs_bacteria.fa'))

# Make a dataframe of the sequences and their ASV ID. 
seq_name_bacteria <- base::names(bacteria_seqs_fasta)
sequence_bacteria <- base::paste(bacteria_seqs_fasta)
bacteria_rep_seqs <- base::data.frame(seq_name_bacteria, sequence_bacteria)

# Join the ASV table and the representative sequences.
bacteria_asv_IDs <- dplyr::left_join(bacteria_rep_seqs, bacteria_asv, by = 'sequence_bacteria')

# Naming the Samples. 
bacteria_asv_IDs <- bacteria_asv_IDs %>% 
  dplyr::rename_with(~base::paste0("Sample_", .), -c(1:2))

# set ASV ID as the rownames.
base::rownames(bacteria_asv_IDs) <- bacteria_asv_IDs$seq_name_bacteria
bacteria_asv_IDs$seq_name_bacteria <- NULL

# Load the DNA concentration necessary for decontam (available as supplementary data). 
bacteria_conc <- utils::read.csv(here::here("Data", 'bacteria_decontam_conc_big.csv'), header = T, sep =',')
base::rownames(bacteria_conc) <- base::paste0("Sample_", bacteria_conc$sample_ID)
base::rownames(bacteria_conc) <- base::sub("[-]", "_", x = base::rownames(bacteria_conc))

# Make the parts of the phyloseq object.
ASV_mat_decontam_fun <- base::data.matrix(bacteria_asv_IDs)
ASV_bacteria_decontam <- phyloseq::otu_table(ASV_mat_decontam_fun, taxa_are_rows = TRUE)
sampledata_bacteria_decontam <- phyloseq::sample_data(bacteria_conc)

# Combine with phyloseq
ps_bacteria_decontam <- phyloseq::phyloseq(ASV_bacteria_decontam, sampledata_bacteria_decontam)
ps_bacteria_decontam

# Check the library sizes of the Samples and Controls.
df_bacteria <- base::as.data.frame(phyloseq::sample_data(ps_bacteria_decontam)) # Put sample_data into a ggplot-friendly data.frame
df_bacteria$LibrarySize <- phyloseq::sample_sums(ps_bacteria_decontam)
df_bacteria <- df_bacteria[base::order(df_bacteria$LibrarySize),]
df_bacteria$Index <- base::seq(base::nrow(df_bacteria))
ggplot2::ggplot(data=df_bacteria, ggplot2::aes(x=Index, y=LibrarySize, color=Sample_or_Control)) +
  ggplot2::geom_point()

# Decontam with combined frequency and prevalence approach. 
phyloseq::sample_data(ps_bacteria_decontam)$is.neg <- phyloseq::sample_data(ps_bacteria_decontam)$Sample_or_Control == "Control Sample"
contam_bacteria_combi <- decontam::isContaminant(ps_bacteria_decontam,
                                       method = 'combined',
                                       neg = 'is.neg',
                                       conc = 'quant_reading')
base::table(contam_bacteria_combi$contaminant)
base::which(contam_bacteria_combi$contaminant)

# Trim the identified contaminants from the phyloseq object. 
ps_bacteria_noncontam <- phyloseq::prune_taxa(!contam_bacteria_combi$contaminant,
                                    ps_bacteria_decontam)
ps_bacteria_noncontam

###
#LULU curation
###

# Remove taxa without reads.
ps_bacteria_noncontam_pruned <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_bacteria_noncontam) > 0,
                                           ps_bacteria_noncontam)
ASV_table_bacteria <- base::as.data.frame(phyloseq::otu_table(ps_bacteria_noncontam_pruned))
bacteria_matchlist <- utils::read.table(here::here("Data", 'match_list_bacteria.txt'),
                                 header = F,
                                 as.is = T,
                                 stringsAsFactors = F)

# Run the LULU algorithm. 

ASV_table_bacteria_cur <- lulu::lulu(ASV_table_bacteria, bacteria_matchlist)

ASV_table_bacteria_cur$curated_count
ASV_table_bacteria_cur$discarded_count

##################################################################
##                          Section 3                           ##
##              Save Data needed for Data Cleaning              ##
##################################################################

saveRDS(ASV_table_algae_cur, "ASV_table_algae_cur.rds")

saveRDS(ASV_table_bacteria_cur, "ASV_table_bacteria_cur.rds")

saveRDS(ASV_table_fungi_cur, "ASV_table_fungi_cur.rds")

saveRDS(algae_asv_IDs, "algae_asv_IDs.rds")

