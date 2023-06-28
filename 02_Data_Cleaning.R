#################################################################
##                          Section 1                          ##
##                       Package Loading                       ##
#################################################################

library(here); packageVersion("here")
#1.0.1
library(dplyr); packageVersion("dplyr")
#1.0.9
library(phyloseq); packageVersion("phyloseq")
#1.40.0
library(stringr); packageVersion("stringr")
#1.4.0
library(GGally); packageVersion("GGally")
#2.1.2
library(heplots); packageVersion("heplots")
#1.4.2
##################################################################
##                          Section 2                           ##
##         Data Loading, Data Cleaning and Initialization       ##
##################################################################

##----------------------------------------------------------------
##                           Metadata                            -
##----------------------------------------------------------------

# Load the metadata table.
own_metadata <- utils::read.csv(here::here("Data","sample_data_2.csv"), header = T, sep = ",")

# Get the EP_Plot_ID from the Sample names
own_metadata <- own_metadata %>% tidyr::separate(sample, c(NA, 'Plot_ID', NA) , sep = '_', remove = F)

# Insert EW after the first character of the string in a given cell. 
own_metadata$Plot_ID <- base::gsub("^(.{1})(.*)$",
                             "\\1EW\\2",
                             own_metadata$Plot_ID)

# Load the Climate data. 
climate <- utils::read.csv(here::here("TestData","plots_climate_weekly.csv"), header = T, sep = ",")

# Subset to the sampling weeks and 2 weeks before for each Explo. Alb = weeks 16-18,
# Hai = weeks 17-19, Sch = weeks 18-20. 
sampling_weeks_clim <- climate %>% 
  dplyr::filter(base::grepl("AEW", plotID) & datetime %in% c("2021W16", "2021W17", "2021W18") |
                base::grepl("HEW", plotID) & datetime %in% c("2021W17", "2021W18", "2021W19") |
                base::grepl("SEW", plotID) & datetime %in% c("2021W18", "2021W19", "2021W20"))

# Remove the datetime column since it is not needed anymore.
sampling_weeks_clim$datetime <- NULL

# Calculate the average of the climate variable for each plot. 
sampling_weeks_clim_avg <- sampling_weeks_clim %>% 
  dplyr::group_by(plotID) %>%  
  dplyr::summarise_each(dplyr::funs(base::round(base::mean(., na.rm = T), 2))) %>% 
  dplyr::rename(Plot_ID = plotID)

# One of the plots (HEW33) had a broken thermometer so we replace the temperature 
# value with that of the nearest station (HEW26).
sampling_weeks_clim_avg$Ta_200[sampling_weeks_clim_avg$Plot_ID == "HEW33"] <-
  sampling_weeks_clim_avg$Ta_200[sampling_weeks_clim_avg$Plot_ID == "HEW26"]

# Load the stand structure dataset. 
stand_structure_full <-  utils::read.csv(here::here("TestData","22766_3_data.csv"), header = T, sep = ";")

# Subset to the variables we are interested in and give them more meaningful names. 
stand_structure_variables <- stand_structure_full %>% 
  dplyr::select(EP_Plotid, ssm_N, ssm_BA, sp_BA_1D, d_m, d_SD, d_gini) %>% 
  dplyr::rename(Plot_ID = EP_Plotid, 
                stand_density_abundance = ssm_N,
                stand_density_basal_area = ssm_BA,
                stand_evenness_basal_area = sp_BA_1D,
                DBH_avg = d_m)

# Load the dataset of the effective number of forest layers = measure of vertical heterogeneity. 
vertical_heterogeneity <- utils::read.csv(here::here("TestData","27826_3_data.csv"), header = T, sep = ";")

# We are only interested in the newest data. 
vertical_heterogeneity_new <- vertical_heterogeneity %>% 
  dplyr::select(EP, enl_2019) %>% 
  dplyr::rename(Plot_ID = EP)

# Load the dataset describing the proportion of forest in a 2 km radius around the plot to see 
# if the plot is situated within a forest or is a smaller patch of trees. 
surroundings <- utils::read.csv(here::here("TestData","15929_2_data.csv"), header = T, sep = ";")

# Subset to only keep info on surrounding forest areas. 
surroundings_forest <- surroundings %>% 
  dplyr::select(Plot, RA_forest) %>% 
  dplyr::rename(Plot_ID = Plot)

# Load the Canopy Openness dataset. 
canopy_openness <- utils::read.csv(here::here("TestData","canopy_openness_2019.csv"), header = T, sep = ",") %>% 
  dplyr::rename(Plot_ID = EP)

# Load in the inventory data of single trees.
single_trees <- utils::read.csv(here::here('TestData', '21426_3_data.csv'), header = T, sep = ';')

# Calculate plot wise number of present tree species. 
plot_struc <- single_trees %>% 
  dplyr::group_by(EP_Plotid) %>% 
  dplyr::count(species) %>% 
  dplyr::rename(number = n, Plot_ID = EP_Plotid)

# Find the most abundant tree species on each plot.
plot_dom <- plot_struc %>% 
  dplyr::group_by(Plot_ID) %>%
  dplyr::mutate(num_dom = max(number)) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(number == num_dom) %>%
  dplyr::select(-number) %>%
  dplyr::rename(dom_species = species)

# Get the total number of trees on each plot. 
plot_trees <- plot_struc %>% 
  dplyr::group_by(Plot_ID) %>%
  dplyr::summarise(num_tot = sum(number))

# Combine both into one table. 
plot_comp <- dplyr::left_join(plot_dom, plot_trees, by = 'Plot_ID')

# Calculate most abundant tree proportion. 
plot_comp <- plot_comp %>%
  dplyr::mutate(dom_tot_ratio = base::round(num_dom/num_tot , 2))

# Remove intermediate columns we are not interested in. 
plot_comp_final <- plot_comp %>% 
  dplyr::select("Plot_ID", "dom_tot_ratio")

# Load the coordinates for the plots. 
geo <- utils::read.csv(here::here("basic_plot_info.csv"), sep = ";") %>% 
  dplyr::filter(Landuse == "Forest") %>% 
  dplyr::filter(EP_Plot_ID != "na") %>% 
  dplyr::rename(Plot = EP_Plot_ID) %>% 
  dplyr::select(Plot, Latitude, Longitude)

# Combine all of our metadata into one big metadata table. 
metadata_full <- plyr::join_all(base::list(own_metadata,
                                     sampling_weeks_clim_avg,
                                     stand_structure_variables,
                                     vertical_heterogeneity_new, 
                                     surroundings_forest, 
                                     canopy_openness, 
                                     plot_comp_final),
                                by = "Plot_ID",
                                type = "inner")

# Set the sample IDs as the rownames.
base::rownames(metadata_full) <- metadata_full$sample
metadata_full$sample <- NULL

# Load in the Algae ASV IDs. 
algae_asv_IDs <- base::readRDS("algae_asv_IDs.rds")

# Get the sample names from the algal ASV table. 
sample_names <- base::colnames(algae_asv_IDs) %>% 
  stringr::str_subset(pattern = ("sequence|MPC|PCN|BL"), negate = T)

# Some Plots had been clear cut and need to be removed. 
# An additional one has very low read counts for algae and fungi and also needs to be removed.
metadata_full <- metadata_full[!(metadata_full$Plot_ID == "HEW3" |
                                   metadata_full$Plot_ID == "HEW13" |
                                   metadata_full$Plot_ID == "HEW12"),]

# Few plots have a dominant tree species that is not beech, pine or fir. 
# This takes away too much statistical power, thus we remove them. 

metadata_full_tree_filtered <- metadata_full %>% 
  dplyr::filter(dominant_tree == "Pinus_sylvestris" | 
                  dominant_tree == "Picea_abies" | 
                  dominant_tree == "Fagus_sylvatica")

# We need to append the coordinates as well for the GDMs later on. 

metadata_full_tree_filtered <- metadata_full_tree_filtered %>% 
  tibble::rownames_to_column(var = "Sample_ID") %>% 
  dplyr::inner_join(., geo %>% rename(Plot_ID = Plot))

# Set the sample names as the rownames 

rownames(metadata_full_tree_filtered) <- metadata_full_tree_filtered$Sample_ID

metadata_full_tree_filtered$Sample_ID <- NULL

##----------------------------------------------------------------
##                          ASV tables                           -
##----------------------------------------------------------------

##---------
##  Algae  
##---------
ASV_table_algae_cur <- base::readRDS("ASV_table_algae_cur.rds")

# Keep only samples that do represent real tree swabs. Cut Controls.
asv_algae <- ASV_table_algae_cur$curated_table %>%
  dplyr::select(all_of(sample_names))

##---------
##  Bacteria  
##---------
ASV_table_bacteria_cur <- base::readRDS("ASV_table_bacteria_cur.rds")

# Keep only samples that do represent real tree swabs. Cut Controls.
asv_bacteria <- ASV_table_bacteria_cur$curated_table %>% 
  dplyr::select(all_of(sample_names))

##---------
##  Fungi  
##---------
ASV_table_fungi_cur <- base::readRDS("ASV_table_fungi_cur.rds")

# Keep only samples that do represent real tree swabs. Cut Controls. 
asv_fungi <- ASV_table_fungi_cur$curated_table %>%
  dplyr::select(all_of(sample_names))

##---------------------------------------------------------------
##                        Taxonomy Tables                       -
##---------------------------------------------------------------

##---------
##  Algae  
##---------

tax_algae <- base::readRDS(here::here("Data", 'tax_table_algae.rds'))

tax_algae <- base::as.data.frame(tax_algae)

tax_algae$taxaID <- NULL

# Filter out any uncultured or environmental sequences. 
tax_algae_no_uncultured <- tax_algae %>% 
  dplyr::filter(!base::grepl('uncultured', species) ) %>% 
  dplyr::filter(!base::grepl('environmental', species) )

# Filter anything under 95 % percent identity.
tax_algae_good_blast <- tax_algae_no_uncultured %>% 
  dplyr::filter(percent_ident > 95.00)

# Keep only the assignment with the highest percent identity for each ASV. 
tax_algae_highest_ident <- tax_algae_good_blast %>% 
  dplyr::group_by(ASV_ID) %>% 
  dplyr::top_n(1, percent_ident)

# Keep everything with a unique taxaID. 
tax_algae_unique_taxaID <- tax_algae_highest_ident %>% 
  dplyr::distinct(taxaID, .keep_all = T)

# Keep only algae (all Chlorphyta, Klebsormidiophyceae, Xanthophyceae, Chrysophyceae).
tax_clean_algae_dups <- tax_algae_unique_taxaID %>% 
  dplyr::filter(phylum == "Chlorophyta" | class %in% c("Klebsormidiophyceae", "Xanthophyceae", "Chrysophyceae"))

# The data still contains duplicates since some ASVs can be assigned with the same level of confidence.
# Remove duplicated ASV_IDs to only keep the top match in NCBI. 
tax_clean_algae <- tax_clean_algae_dups %>% 
  dplyr::distinct(ASV_ID, .keep_all = T) %>% 
  base::data.frame()

# Set the ASV ID as the rownames.
base::row.names(tax_clean_algae) <- tax_clean_algae$ASV_ID

#Remove all columns that are unnecessary now. 
tax_clean_algae[, c("ASV_ID","accession", "percent_ident", "length", "mismatches", 
                    "gapopen", "qsstart", "qend", "subject_start", 
                    "subject_end", "evalue", "bitscore", "taxaID")] <- base::list(NULL)

# Add a kingdom column. 
tax_clean_algae <- base::cbind(Kingdom = "Viridiplantae", tax_clean_algae)

# remove Superkingdom. 
tax_clean_algae$superkingdom <- NULL

# NCBI writes the headers in lowercase. We change that to make it comparable with the two other groups. 
tax_clean_algae <- tax_clean_algae %>% 
  dplyr::rename_with(stringr::str_to_title)

# Load algael sequences. 
algae_seqs_fasta <- Biostrings::readDNAStringSet(here::here("Data", 'curated_ASVs_algae.fa'))

# Make a dataframe of the sequences and their ASV ID. 
seq_name_algae <- base::names(algae_seqs_fasta)
sequence_algae <- base::paste(algae_seqs_fasta)
algae_rep_seqs <- base::data.frame(seq_name_algae, sequence_algae)

##---------
##  Bacteria  
##---------

# Load the bacterial taxonomy table. 
# (Available as supplementary data)
tax_bacteria <- base::readRDS(here::here("Data", 'tax_table_bacteria.rds'))
tax_bacteria <- base::as.data.frame(tax_bacteria) %>%
  tibble::rownames_to_column('sequence')
tax_bacteria <- tax_bacteria %>%
  dplyr::rename(sequence_bacteria = sequence)

# Load bacterial sequences. 
bacteria_seqs_fasta <- Biostrings::readDNAStringSet(here::here("Data", 'ASVs_bacteria.fa'))

# Make a dataframe of the sequences and their ASV ID. 
seq_name_bacteria <- base::names(bacteria_seqs_fasta)
sequence_bacteria <- base::paste(bacteria_seqs_fasta)
bacteria_rep_seqs <- base::data.frame(seq_name_bacteria, sequence_bacteria)

# Join the taxonomy table and the representative sequences
tax_clean_bacteria <- dplyr::left_join(tax_bacteria, bacteria_rep_seqs, by = 'sequence_bacteria')

# Set rownames.
base::row.names(tax_clean_bacteria) <- tax_clean_bacteria$seq_name_bacteria

# Remove any Chloroplast and mitochondrial Sequences
bacteria_tax_fin_raw <- dplyr::filter(tax_clean_bacteria, Order != 'Chloroplast')
bacteria_tax_fin_raw <- dplyr::filter(tax_clean_bacteria, Family != 'Mitochondria')

bacteria_tax_fin_raw$seq_name_bacteria <- NULL
bacteria_tax_fin_raw$sequence_bacteria <- NULL

##---------
##  Fungi  
##---------

# Load the fungal taxonomy table.
# (Available as supplementary data)
tax_fungi <- base::readRDS(here::here("Data", 'tax_table_fungi_new.rds'))
tax_fungi <- base::as.data.frame(tax_fungi) %>%
  tibble::rownames_to_column('sequence')
tax_fungi <- tax_fungi %>%
  dplyr::rename(sequence_fungi = sequence)

# Load the fungal reads.
fungi_seqs_fasta <- Biostrings::readDNAStringSet(here::here("Data", 'ASVs_fungi.fa'))

# Make a dataframe of the sequences and their ASV ID. 
seq_name_fungi <- base::names(fungi_seqs_fasta)
sequence_fungi <- base::paste(fungi_seqs_fasta)
fungi_rep_seqs <- base::data.frame(seq_name_fungi, sequence_fungi)

# Join the taxonomy table and the representative sequences
tax_clean_fungi <- dplyr::left_join(tax_fungi, fungi_rep_seqs, by = 'sequence_fungi')

# Split the taxonomy into different columns of taxonomic levels.

fungi_tax_fin <- tidyr::separate(tax_clean_fungi, Kingdom, c(NA, 'Kingdom') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Phylum, c(NA, 'Phylum') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Class, c(NA, 'Class') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Order, c(NA, 'Order') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Family, c(NA, 'Family') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Genus, c(NA, 'Genus') , sep = '__')
fungi_tax_fin <- tidyr::separate(fungi_tax_fin, Species, c(NA, 'Species') , sep = '__')

# Keep only Fungi 
fungi_tax_fin <- fungi_tax_fin %>% 
  dplyr::filter(Kingdom == "Fungi")

# Rename the ASV_ID column. 
fungi_tax_fin <- dplyr::rename(fungi_tax_fin, ASV_ID = seq_name_fungi)

# Set rownames.
base::row.names(fungi_tax_fin) <- fungi_tax_fin$ASV_ID

fungi_tax_fin$sequence_fungi <- NULL
fungi_tax_fin$ASV_ID <- NULL

##---------------------------------------------------------------
##               Create the Phyloseq Objects                    -
##---------------------------------------------------------------

##---------
##  Algae  
##---------

# Transform dataframe to matrix.
asvmat_algae <- base::as.matrix(asv_algae) 

# Transform dataframe to matrix.
taxmat_algae <- base::as.matrix(tax_clean_algae) 

# Create ASV table for phyloseq.
ASV_ALG <- phyloseq::otu_table(asvmat_algae, taxa_are_rows = T)

# Create taxonomy table for phyloseq.
TAX_ALG <- phyloseq::tax_table(taxmat_algae) 

# Metadata for phyloseq. 
sampledata <- phyloseq::sample_data(metadata_full_tree_filtered) 

# Combine in phyloseq object. 
phy_algae <- phyloseq::phyloseq(ASV_ALG, TAX_ALG, sampledata) 
phy_algae

# Remove taxa that have zero reads remaining. 
phy_algae <- phyloseq::prune_taxa(phyloseq::taxa_sums(phy_algae) != 0, phy_algae)

##---------
##  Bacteria  
##---------

# Transform dataframe to matrix.
asvmat_bacteria <- base::as.matrix(asv_bacteria) 

# Transform dataframe to matrix.
taxmat_bacteria <- base::as.matrix(bacteria_tax_fin_raw) 

# Create OTU table for phyloseq.
ASV_BAC <- phyloseq::otu_table(asvmat_bacteria, taxa_are_rows = T) 

# Create taxonomy table for phyloseq.
TAX_BAC <- phyloseq::tax_table(taxmat_bacteria) 

# Metadata for phyloseq.
sampledata <- phyloseq::sample_data(metadata_full_tree_filtered)  

# Combine in phyloseq object. 
phy_bacteria <- phyloseq::phyloseq(ASV_BAC, TAX_BAC, sampledata) 
phy_bacteria

# Remove taxa that have zero reads remaining. 
phy_bacteria <- phyloseq::prune_taxa(phyloseq::taxa_sums(phy_bacteria) != 0, phy_bacteria)

##---------
##  Fungi  
##---------

# Transform dataframe to matrix.
asvmat_fungi <- base::as.matrix(asv_fungi)

# Transform dataframe to matrix.
taxmat_fungi <- base::as.matrix(fungi_tax_fin) 

# Create ASV table for phyloseq.
ASV_FUN <- phyloseq::otu_table(asvmat_fungi, taxa_are_rows = T) 

# Create taxonomy table for phyloseq.
TAX_FUN <- phyloseq::tax_table(taxmat_fungi) 

# Metadata for phyloseq. 
sampledata <- phyloseq::sample_data(metadata_full_tree_filtered) 

# Combine in phyloseq object. 
phy_fungi <- phyloseq::phyloseq(ASV_FUN, TAX_FUN, sampledata) 
phy_fungi

# Remove taxa that have zero reads remaining. 
phy_fungi <- phyloseq::prune_taxa(phyloseq::taxa_sums(phy_fungi) != 0, phy_fungi)

# Subset the phyloseq object to split the dataset into the bark and soil samples. 
##---------
##  Algae  
##---------

phy_algae_bark <- phyloseq::subset_samples(phy_algae, substrate == "bark") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

phy_algae_soil <- phyloseq::subset_samples(phy_algae, substrate == "soil") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

##---------
##  Bacteria  
##---------

phy_bacteria_bark <- phyloseq::subset_samples(phy_bacteria, substrate == "bark") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

phy_bacteria_soil <- phyloseq::subset_samples(phy_bacteria, substrate == "soil") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)


##---------
##  Fungi  
##---------

phy_fungi_bark <- phyloseq::subset_samples(phy_fungi, substrate == "bark")  %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)

phy_fungi_soil <- phyloseq::subset_samples(phy_fungi, substrate == "soil") %>% 
  phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0,.)


# Check rarefaction curves for sample completeness. 
#vegan::rarecurve(t(phyloseq::otu_table(phy_algae_bark)), label = F)

#vegan::rarecurve(t(phyloseq::otu_table(phy_algae_soil)), label = F)

#vegan::rarecurve(t(phyloseq::otu_table(phy_bacteria_bark)))

#vegan::rarecurve(t(phyloseq::otu_table(phy_bacteria_soil)))

#vegan::rarecurve(t(phyloseq::otu_table(phy_fungi_bark)))

#vegan::rarecurve(t(phyloseq::otu_table(phy_fungi_soil)))

##---------------------------------------------------------------
##     Subset hypothesized variables based on Correlations      -
##---------------------------------------------------------------

# Pairs plot 
# Select only the variables.
# Remove soil samples. 
metadata_pairs <- metadata_full_tree_filtered %>% 
  dplyr::filter(substrate == "bark") %>% 
  dplyr::select(-one_of("Plot_ID", "substrate")) 

#GGally::ggpairs(metadata_pairs)

# Have a second look at the pairs plot only including variables that we want to include.
# So we remove others that were correlated etc. 
metadata_pairs_2 <- metadata_pairs %>%  
  dplyr::select(-one_of("tree_type", "stand_density_abundance", "precipitation_radolan", "enl_2019",
                        "stand_evenness_basal_area", "d_SD", "PAR"))

#GGally::ggpairs(metadata_pairs_2)

# Correlation between categorical and continous variables is not straightforward.
# We are interested in the individual influence of climate and stand variables, but they are 
# highly related to the tree type. Excactly how much we want to find out so we can justify excluding 
# tree type from the analysis. 

manova <- stats::manova(cbind(stand_density_basal_area,
                              DBH_avg, d_gini, 
                              canopy_openness_2019, dom_tot_ratio,
                              rH_200, Ta_200) ~ dominant_tree, 
                       data = metadata_pairs_2)

summary(manova)


# Plot the correlation of tree species and tree dependent variables. 
metadata_pairs_2 %>% dplyr::select("dominant_tree", "rH_200", "Ta_200", "stand_density_basal_area",
                                   "DBH_avg", "d_gini", "canopy_openness_2019", "dom_tot_ratio")  
  GGally::ggpairs()

#################################################################
##                          Section 4                          ##
##                         Data Saving                         ##
#################################################################

base::saveRDS(metadata_bark, "metadata_bark.rds")

base::saveRDS(phy_algae_bark, "phy_algae_bark.rds")

base::saveRDS(phy_bacteria_bark, "phy_bacteria_bark.rds")

base::saveRDS(phy_fungi_bark, "phy_fungi_bark.rds")

#################################################################
##                          Section 5                          ##
##                   Miscallenous Data Export                  ##
#################################################################

###
# Species Lists
###
algae_species_list <- metagMisc::phyloseq_to_df(phy_algae_bark, addtax = T, sorting = "taxonomy") %>% 
  dplyr::select("OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") %>% 
  dplyr::left_join(., dplyr::rename(algae_rep_seqs, OTU = seq_name_algae)) %>% 
  dplyr::rename(Sequence = sequence_algae, ASV_ID = OTU)

fungi_species_list <- metagMisc::phyloseq_to_df(phy_fungi_bark, addtax = T, sorting = "taxonomy") %>% 
  dplyr::select("OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") %>% 
  dplyr::left_join(., dplyr::rename(fungi_rep_seqs, OTU = seq_name_fungi)) %>% 
  dplyr::rename(Sequence = sequence_fungi, ASV_ID = OTU)

bacteria_species_list <- metagMisc::phyloseq_to_df(phy_bacteria_bark, addtax = T, sorting = "taxonomy") %>% 
  dplyr::select("OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus") %>% 
  dplyr::left_join(., dplyr::rename(bacteria_rep_seqs, OTU = seq_name_bacteria)) %>% 
  dplyr::rename(Sequence = sequence_bacteria, ASV_ID = OTU)

write.csv(algae_species_list, "algae_species_list.csv", row.names = F)
write.csv(fungi_species_list, "fungi_species_list.csv", row.names = F)
write.csv(bacteria_species_list, "bacteria_species_list.csv", row.names = F)

###
# ASV tables
###



###
# Fasta Files for Genbank
###

# Filter the original sequences to contain only sequences from the curation and after removing non-algae.
filtered_algae_reps_bark <- algae_rep_seqs %>% 
  dplyr::filter(seq_name_algae %in% taxa_names(phy_algae_bark)) %>% 
  dplyr::rename(name = seq_name_algae,
                seq = sequence_algae)

filtered_algae_reps_soil <- algae_rep_seqs %>% 
  dplyr::filter(seq_name_algae %in% taxa_names(phy_algae_soil)) %>% 
  dplyr::rename(name = seq_name_algae,
                seq = sequence_algae)

# Filter the original sequences to contain only sequences from the curation and after removing non-fungi.
filtered_fungi_reps_bark <- fungi_rep_seqs %>% 
  dplyr::filter(seq_name_fungi %in% taxa_names(phy_fungi_bark)) %>% 
  dplyr::rename(name = seq_name_fungi,
                seq = sequence_fungi)

filtered_fungi_reps_soil <- fungi_rep_seqs %>% 
  dplyr::filter(seq_name_fungi %in% taxa_names(phy_fungi_soil)) %>% 
  dplyr::rename(name = seq_name_fungi,
                seq = sequence_fungi)

# Filter the original sequences to contain only sequences from the curation and after removing non-bacteria.
filtered_bacteria_reps_bark <- bacteria_rep_seqs %>% 
  dplyr::filter(seq_name_bacteria %in% taxa_names(phy_bacteria_bark)) %>% 
  dplyr::rename(name = seq_name_bacteria,
                seq = sequence_bacteria)

filtered_bacteria_reps_soil <- bacteria_rep_seqs %>% 
  dplyr::filter(seq_name_bacteria %in% taxa_names(phy_bacteria_soil)) %>% 
  dplyr::rename(name = seq_name_bacteria,
                seq = sequence_bacteria)

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}


writeFasta(filtered_algae_reps_bark, "ASVs_filtered_algae_bark.fa")

writeFasta(filtered_fungi_reps_bark, "ASVs_filtered_fungi_bark.fa")

writeFasta(filtered_bacteria_reps_bark, "ASVs_filtered_bacteria_bark.fa")

writeFasta(filtered_algae_reps_soil, "ASVs_filtered_algae_soil.fa")

writeFasta(filtered_fungi_reps_soil, "ASVs_filtered_fungi_soil.fa")

writeFasta(filtered_bacteria_reps_soil, "ASVs_filtered_bacteria_soil.fa")
