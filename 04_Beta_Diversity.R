#################################################################
##                          Section 1                          ##
##                       Package Loading                       ##
#################################################################

library(here); packageVersion("here")
#1.0.1
library(vegan); packageVersion("vegan")
#2.6.2
library(tidyverse); packageVersion("tidyverse")
#1.3.2
library(hillR); packageVersion("hillR")
#0.5.1
library(gdm); packageVersion("gdm")
#1.5.0.7
library(betapart); packageVersion("betapart")
#1.5.6
library(phyloseq); packageVersion("phyloseq")
#1.40.0
library(effects); packageVersion("effects")
#4.2.1

##----------------------------------------------------------------
##                        Custom Functions                       -
##----------------------------------------------------------------

# Calculates the Sorensen diversity partitions of two communities. 
# See Alberdi & Gilbert 2019 MolEcolResources and Chao, Chiu & Hsieh 2012 in Ecology.
# Sorensen type overlap. Here for community similarity, the GDM needs dissimilarity so we need to 
# take the one-complements 1-CqN. 
by_hand_cqn <- function(beta, q, N){(((1/beta)^(q-1)) - ((1/N)^(q-1))) / (1-((1/N)^(q-1)))}

# When q = 1 the formula is not defined. Therefore we need to set the Horn Index as an approximation
# of q approaching 1. See Chao, Chiu & Hsieh 2012 in Ecology.
cqn_q1 <- function(beta, N){1-log(beta)/log(N)}

#################################################################
##                          Section 2                          ##
##                        Data Loading                         ##
#################################################################

phy_algae_bark <- base::readRDS("phy_algae_bark.rds")

phy_bacteria_bark <- base::readRDS("phy_bacteria_bark.rds")

phy_fungi_bark <- base::readRDS("phy_fungi_bark.rds")

metadata_bark <- base::readRDS("metadata_bark.rds")

metadata_full_tree_filtered_gdm <- base::readRDS("metadata_full_tree_filtered_gdm.rds")

#################################################################
##                          Section 3                          ##
##    Generalized Dissimilarity Models of Beta Diversity       ##
#################################################################

# Here we fit Generalized Dissimilarity Models (GDM) of the Hill Numbers for beta diversity. 
# It allows us to find the drivers of beta diversity along environmental gradients. 
# For a nice guide see Mokany et. al 2022 in Global Ecology and Biogeography


# First we need to calculate the beta diversity Hill numbers for our plot pairs. 
# q = 0 is sorensen index of dissimilarity.
# q = 1 is Horn index of dissimilarity. 
# q = 2 is Morisita-Horn index of dissimilarity.

# Algae
asv_algae_hill_bark <- base::data.frame(phyloseq::otu_table(phy_algae_bark)) 

alg_q0_bark_beta <- hillR::hill_taxa_parti_pairwise(comm = t(asv_algae_hill_bark),
                                                    q = 0, pairs = "full")
alg_q1_bark_beta <- hillR::hill_taxa_parti_pairwise(comm = t(asv_algae_hill_bark),
                                                    q = 1, pairs = "full")
alg_q2_bark_beta <- hillR::hill_taxa_parti_pairwise(comm = t(asv_algae_hill_bark),
                                                    q = 2, pairs = "full")
# Bacteria
asv_bacteria_hill_bark <- base::data.frame(phyloseq::otu_table(phy_bacteria_bark)) 

bac_q0_bark_beta <- hillR::hill_taxa_parti_pairwise(comm = t(asv_bacteria_hill_bark),
                                                    q = 0, pairs = "full")
bac_q1_bark_beta <- hillR::hill_taxa_parti_pairwise(comm = t(asv_bacteria_hill_bark),
                                                    q = 1, pairs = "full")
bac_q2_bark_beta <- hillR::hill_taxa_parti_pairwise(comm = t(asv_bacteria_hill_bark),
                                                    q = 2, pairs = "full")

# Fungi
asv_fungi_hill_bark <- base::data.frame(phyloseq::otu_table(phy_fungi_bark)) 

fun_q0_bark_beta <- hillR::hill_taxa_parti_pairwise(comm = t(asv_fungi_hill_bark),
                                                    q = 0, pairs = "full")
fun_q1_bark_beta <- hillR::hill_taxa_parti_pairwise(comm = t(asv_fungi_hill_bark),
                                                    q = 1, pairs = "full")
fun_q2_bark_beta <- hillR::hill_taxa_parti_pairwise(comm = t(asv_fungi_hill_bark),
                                                    q = 2, pairs = "full")

# Then we need to transform them into community dissimilarities so we can plug them into the GDM. 
# Calculate Sorensen Type Overlap as community dissimilarities. 

# Algae
alg_q0_bark_beta_complement_cqn <- base::data.frame(site1 = bac_q0_bark_beta$site1,
                                              site2 = bac_q0_bark_beta$site2,
                                              alg_q0_cqn = (1 - by_hand_cqn(alg_q0_bark_beta$TD_beta, q = 0, N = 2)))

alg_q1_bark_beta_complement_cqn <- base::data.frame(site1 = bac_q0_bark_beta$site1,
                                              site2 = bac_q0_bark_beta$site2,
                                              alg_q1_cqn = 1 - cqn_q1(beta = alg_q1_bark_beta$TD_beta, N = 2))

alg_q2_bark_beta_complement_cqn <- base::data.frame(site1 = bac_q0_bark_beta$site1,
                                              site2 = bac_q0_bark_beta$site2,
                                              alg_q2_cqn = (1 - by_hand_cqn(alg_q2_bark_beta$TD_beta, q = 2, N = 2)))

# Bacteria
bac_q0_bark_beta_complement_cqn <- base::data.frame(site1 = bac_q0_bark_beta$site1,
                                              site2 = bac_q0_bark_beta$site2,
                                              bac_q0_cqn = (1 - by_hand_cqn(bac_q0_bark_beta$TD_beta, q = 0, N = 2)))

bac_q1_bark_beta_complement_cqn <- base::data.frame(site1 = bac_q0_bark_beta$site1,
                                              site2 = bac_q0_bark_beta$site2,
                                              bac_q1_cqn = 1 - cqn_q1(beta = bac_q1_bark_beta$TD_beta, N = 2))

bac_q2_bark_beta_complement_cqn <- base::data.frame(site1 = bac_q0_bark_beta$site1,
                                              site2 = bac_q0_bark_beta$site2,
                                              bac_q2_cqn = (1 - by_hand_cqn(bac_q2_bark_beta$TD_beta, q = 2, N = 2)))

# Fungi
fun_q0_bark_beta_complement_cqn <- base::data.frame(site1 = bac_q0_bark_beta$site1,
                                              site2 = bac_q0_bark_beta$site2,
                                              fun_q0_cqn = (1 - by_hand_cqn(fun_q0_bark_beta$TD_beta, q = 0, N = 2)))

fun_q1_bark_beta_complement_cqn <- base::data.frame(site1 = bac_q0_bark_beta$site1,
                                              site2 = bac_q0_bark_beta$site2,
                                              fun_q1_cqn = 1 - cqn_q1(beta = fun_q1_bark_beta$TD_beta, N = 2))

fun_q2_bark_beta_complement_cqn <- base::data.frame(site1 = bac_q0_bark_beta$site1,
                                              site2 = bac_q0_bark_beta$site2,
                                              fun_q2_cqn = (1 - by_hand_cqn(fun_q2_bark_beta$TD_beta, q = 2, N = 2)))


# Get the predictors in a format we need for the GDM.
metadata_bark_beta <- dplyr::rename(metadata_bark, 
                                    site1 = Sample_ID,
                                    yCoord = Latitude,
                                    xCoord = Longitude) %>% 
  dplyr::select(., c("site1", "rH_200", "Ta_200","stand_density_basal_area", 
                     "DBH_avg", "d_gini", "RA_forest", "canopy_openness_2019",
                     "dom_tot_ratio", "yCoord", "xCoord")) %>% 
  dplyr::mutate(site2 = site1)

# Get the geographical data in a format we need for the GDM.
metadata_bark_geo <- metadata_bark  %>% 
  dplyr::select(., c("Sample_ID", "Longitude", "Latitude")) %>% 
  dplyr::rename(site1 = Sample_ID) %>% 
  dplyr::mutate(site2 = site1) %>% 
  dplyr::rename(yCoord = Latitude,
                xCoord = Longitude)

# Combine everything in one long table per value of q to go into the GDM. 

# q = 0 
data_gdm_q0 <- alg_q0_bark_beta_complement_cqn %>% 
  dplyr::inner_join(., bac_q0_bark_beta_complement_cqn, by = c("site1", "site2")) %>% 
  dplyr::inner_join(., fun_q0_bark_beta_complement_cqn, by = c("site1", "site2")) %>% 
  dplyr::left_join(., metadata_bark_beta %>% 
                     dplyr::select("site1", "rH_200", "Ta_200","stand_density_basal_area", 
                                   "DBH_avg", "d_gini", "RA_forest", "canopy_openness_2019",
                                   "dom_tot_ratio", "yCoord", "xCoord"), by = "site1") %>% 
  dplyr::rename(s1.rH_200 = rH_200,
                s1.Ta_200 = Ta_200,
                s1.stand_density_basal_area = stand_density_basal_area,
                s1.DBH_avg = DBH_avg,
                s1.d_gini = d_gini,
                s1.RA_forest = RA_forest,
                s1.canopy_openness_2019 = canopy_openness_2019,
                s1.dom_tot_ratio = dom_tot_ratio,
                s1.yCoord = yCoord,
                s1.xCoord = xCoord)  %>% 
  dplyr::left_join(., metadata_bark_beta %>% 
                     dplyr::select("site2", "rH_200", "Ta_200","stand_density_basal_area", 
                                   "DBH_avg", "d_gini", "RA_forest", "canopy_openness_2019",
                                   "dom_tot_ratio", "yCoord", "xCoord"), by = "site2") %>% 
  dplyr::rename(s2.rH_200 = rH_200,
                s2.Ta_200 =Ta_200,
                s2.stand_density_basal_area = stand_density_basal_area,
                s2.DBH_avg = DBH_avg,
                s2.d_gini = d_gini,
                s2.RA_forest = RA_forest,
                s2.canopy_openness_2019 = canopy_openness_2019,
                s2.dom_tot_ratio = dom_tot_ratio,
                s2.yCoord = yCoord,
                s2.xCoord = xCoord) 

# q = 1
data_gdm_q1 <- alg_q1_bark_beta_complement_cqn %>% 
  dplyr::inner_join(., bac_q1_bark_beta_complement_cqn, by = c("site1", "site2")) %>% 
  dplyr::inner_join(., fun_q1_bark_beta_complement_cqn, by = c("site1", "site2")) %>% 
  dplyr::left_join(., metadata_bark_beta %>% 
                     dplyr::select("site1", "rH_200", "Ta_200","stand_density_basal_area", 
                                   "DBH_avg", "d_gini", "RA_forest", "canopy_openness_2019",
                                   "dom_tot_ratio", "yCoord", "xCoord"), by = "site1") %>% 
  dplyr::rename(s1.rH_200 = rH_200,
                s1.Ta_200 = Ta_200,
                s1.stand_density_basal_area = stand_density_basal_area,
                s1.DBH_avg = DBH_avg,
                s1.d_gini = d_gini,
                s1.RA_forest = RA_forest,
                s1.canopy_openness_2019 = canopy_openness_2019,
                s1.dom_tot_ratio = dom_tot_ratio,
                s1.yCoord = yCoord,
                s1.xCoord = xCoord)  %>% 
  dplyr::left_join(., metadata_bark_beta %>% 
                     dplyr::select("site2", "rH_200", "Ta_200","stand_density_basal_area", 
                                   "DBH_avg", "d_gini", "RA_forest", "canopy_openness_2019",
                                   "dom_tot_ratio", "yCoord", "xCoord"), by = "site2") %>% 
  dplyr::rename(s2.rH_200 = rH_200,
                s2.Ta_200 =Ta_200,
                s2.stand_density_basal_area = stand_density_basal_area,
                s2.DBH_avg = DBH_avg,
                s2.d_gini = d_gini,
                s2.RA_forest = RA_forest,
                s2.canopy_openness_2019 = canopy_openness_2019,
                s2.dom_tot_ratio = dom_tot_ratio,
                s2.yCoord = yCoord,
                s2.xCoord = xCoord) 

# q = 2
data_gdm_q2 <- alg_q2_bark_beta_complement_cqn %>% 
  dplyr::inner_join(., bac_q2_bark_beta_complement_cqn, by = c("site1", "site2")) %>% 
  dplyr::inner_join(., fun_q2_bark_beta_complement_cqn, by = c("site1", "site2")) %>% 
  dplyr::left_join(., metadata_bark_beta %>% 
                     dplyr::select("site1", "rH_200", "Ta_200","stand_density_basal_area", 
                                   "DBH_avg", "d_gini", "RA_forest", "canopy_openness_2019",
                                   "dom_tot_ratio", "yCoord", "xCoord"), by = "site1") %>% 
  dplyr::rename(s1.rH_200 = rH_200,
                s1.Ta_200 = Ta_200,
                s1.stand_density_basal_area = stand_density_basal_area,
                s1.DBH_avg = DBH_avg,
                s1.d_gini = d_gini,
                s1.RA_forest = RA_forest,
                s1.canopy_openness_2019 = canopy_openness_2019,
                s1.dom_tot_ratio = dom_tot_ratio,
                s1.yCoord = yCoord,
                s1.xCoord = xCoord)  %>% 
  dplyr::left_join(., metadata_bark_beta %>% 
                     dplyr::select("site2", "rH_200", "Ta_200","stand_density_basal_area", 
                                   "DBH_avg", "d_gini", "RA_forest", "canopy_openness_2019",
                                   "dom_tot_ratio", "yCoord", "xCoord"), by = "site2") %>% 
  dplyr::rename(s2.rH_200 = rH_200,
                s2.Ta_200 =Ta_200,
                s2.stand_density_basal_area = stand_density_basal_area,
                s2.DBH_avg = DBH_avg,
                s2.d_gini = d_gini,
                s2.RA_forest = RA_forest,
                s2.canopy_openness_2019 = canopy_openness_2019,
                s2.dom_tot_ratio = dom_tot_ratio,
                s2.yCoord = yCoord,
                s2.xCoord = xCoord) 

# For the GDM we only need each combination of site1 & site2 once so we need to remove
# the rest from our input table. 

# For that we code a small workaround. 
dist <- base::data.frame(resp_test = base::seq(1, base::ncol(asv_bacteria_hill_bark)))
base::rownames(dist) <- base::factor(base::colnames(asv_bacteria_hill_bark))
h <- base::as.matrix(vegan::vegdist(dist))
h[!base::lower.tri(h)] <- 1000
h <- reshape2::melt(h, value.name = "testpredictor")
h <- base::data.frame(h)
h <- h[!h$testpredictor == 1000,] 
helper <- base::data.frame(helper = base::paste0(h$Var1, h$Var2))

# Keep only the needed pairs in the table. 

data_gdm_q0 <- data_gdm_q0 %>% 
  dplyr::mutate(helper_col = paste0(site1, site2)) %>% 
  dplyr::filter(helper_col %in% helper$helper) %>% 
  dplyr::select(!"helper_col")

data_gdm_q1 <- data_gdm_q1 %>% 
  dplyr::mutate(helper_col = paste0(site1, site2)) %>% 
  dplyr::filter(helper_col %in% helper$helper) %>% 
  dplyr::select(!"helper_col")

data_gdm_q2 <- data_gdm_q2 %>% 
  dplyr::mutate(helper_col = paste0(site1, site2)) %>% 
  dplyr::filter(helper_col %in% helper$helper) %>% 
  dplyr::select(!"helper_col")

# Get the tables per organismal group.
# Algae
data_gdm_q0_alg <- data_gdm_q0 %>% 
  dplyr::rename(distance = alg_q0_cqn,
                s2.bac_q0_cqn = bac_q0_cqn,
                s2.fun_q0_cqn = fun_q0_cqn) %>% 
  tibble::add_column(., s1.bac_q0_cqn = 0,
                     s1.fun_q0_cqn = 0, 
                     weights = 1) %>% 
  dplyr::select("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                "s1.bac_q0_cqn", "s1.fun_q0_cqn",
                "s1.rH_200", "s1.Ta_200", "s1.stand_density_basal_area",
                "s1.DBH_avg", "s1.d_gini", "s1.RA_forest",
                "s1.canopy_openness_2019", "s1.dom_tot_ratio", 
                "s2.bac_q0_cqn", "s2.fun_q0_cqn",
                "s2.rH_200", "s2.Ta_200", "s2.stand_density_basal_area",
                "s2.DBH_avg", "s2.d_gini", "s2.RA_forest",
                "s2.canopy_openness_2019", "s2.dom_tot_ratio") %>% 
  gdm::formatsitepair(bioData = ., bioFormat = 4, predData = .)


data_gdm_q1_alg <- data_gdm_q1 %>% 
  dplyr::rename(distance = alg_q1_cqn,
                s2.bac_q1_cqn = bac_q1_cqn,
                s2.fun_q1_cqn = fun_q1_cqn) %>% 
  tibble::add_column(., s1.bac_q1_cqn = 0,
                     s1.fun_q1_cqn = 0, 
                     weights = 1) %>% 
  dplyr::select("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                "s1.bac_q1_cqn", "s1.fun_q1_cqn",
                "s1.rH_200", "s1.Ta_200", "s1.stand_density_basal_area",
                "s1.DBH_avg", "s1.d_gini", "s1.RA_forest",
                "s1.canopy_openness_2019", "s1.dom_tot_ratio", 
                "s2.bac_q1_cqn", "s2.fun_q1_cqn",
                "s2.rH_200", "s2.Ta_200", "s2.stand_density_basal_area",
                "s2.DBH_avg", "s2.d_gini", "s2.RA_forest",
                "s2.canopy_openness_2019", "s2.dom_tot_ratio") %>% 
  gdm::formatsitepair(bioData = ., bioFormat = 4, predData = .)

data_gdm_q2_alg <- data_gdm_q2 %>% 
  dplyr::rename(distance = alg_q2_cqn,
                s2.bac_q2_cqn = bac_q2_cqn,
                s2.fun_q2_cqn = fun_q2_cqn) %>% 
  tibble::add_column(., s1.bac_q2_cqn = 0,
                     s1.fun_q2_cqn = 0, 
                     weights = 1) %>% 
  dplyr::select("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                "s1.bac_q2_cqn", "s1.fun_q2_cqn",
                "s1.rH_200", "s1.Ta_200", "s1.stand_density_basal_area",
                "s1.DBH_avg", "s1.d_gini", "s1.RA_forest",
                "s1.canopy_openness_2019", "s1.dom_tot_ratio", 
                "s2.bac_q2_cqn", "s2.fun_q2_cqn",
                "s2.rH_200", "s2.Ta_200", "s2.stand_density_basal_area",
                "s2.DBH_avg", "s2.d_gini", "s2.RA_forest",
                "s2.canopy_openness_2019", "s2.dom_tot_ratio") %>% 
  gdm::formatsitepair(bioData = ., bioFormat = 4, predData = .)

# Bacteria
data_gdm_q0_bac <- data_gdm_q0 %>% 
  dplyr::rename(distance = bac_q0_cqn,
                s2.alg_q0_cqn = alg_q0_cqn,
                s2.fun_q0_cqn = fun_q0_cqn) %>% 
  tibble::add_column(., s1.alg_q0_cqn = 0,
                     s1.fun_q0_cqn = 0, 
                     weights = 1) %>% 
  dplyr::select("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                "s1.alg_q0_cqn", "s1.fun_q0_cqn",
                "s1.rH_200", "s1.Ta_200", "s1.stand_density_basal_area",
                "s1.DBH_avg", "s1.d_gini", "s1.RA_forest",
                "s1.canopy_openness_2019", "s1.dom_tot_ratio", 
                "s2.alg_q0_cqn", "s2.fun_q0_cqn",
                "s2.rH_200", "s2.Ta_200", "s2.stand_density_basal_area",
                "s2.DBH_avg", "s2.d_gini", "s2.RA_forest",
                "s2.canopy_openness_2019", "s2.dom_tot_ratio") %>% 
  gdm::formatsitepair(bioData = ., bioFormat = 4, predData = .)


data_gdm_q1_bac <- data_gdm_q1 %>% 
  dplyr::rename(distance = bac_q1_cqn,
                s2.alg_q1_cqn = alg_q1_cqn,
                s2.fun_q1_cqn = fun_q1_cqn) %>% 
  tibble::add_column(., s1.alg_q1_cqn = 0,
                     s1.fun_q1_cqn = 0, 
                     weights = 1) %>% 
  dplyr::select("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                "s1.alg_q1_cqn", "s1.fun_q1_cqn",
                "s1.rH_200", "s1.Ta_200", "s1.stand_density_basal_area",
                "s1.DBH_avg", "s1.d_gini", "s1.RA_forest",
                "s1.canopy_openness_2019", "s1.dom_tot_ratio", 
                "s2.alg_q1_cqn", "s2.fun_q1_cqn",
                "s2.rH_200", "s2.Ta_200", "s2.stand_density_basal_area",
                "s2.DBH_avg", "s2.d_gini", "s2.RA_forest",
                "s2.canopy_openness_2019", "s2.dom_tot_ratio")%>% 
  gdm::formatsitepair(bioData = ., bioFormat = 4, predData = .)

data_gdm_q2_bac <- data_gdm_q2 %>% 
  dplyr::rename(distance = bac_q2_cqn,
                s2.alg_q2_cqn = alg_q2_cqn,
                s2.fun_q2_cqn = fun_q2_cqn) %>% 
  tibble::add_column(., s1.alg_q2_cqn = 0,
                     s1.fun_q2_cqn = 0, 
                     weights = 1) %>% 
  dplyr::select("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                "s1.alg_q2_cqn", "s1.fun_q2_cqn",
                "s1.rH_200", "s1.Ta_200", "s1.stand_density_basal_area",
                "s1.DBH_avg", "s1.d_gini", "s1.RA_forest",
                "s1.canopy_openness_2019", "s1.dom_tot_ratio", 
                "s2.alg_q2_cqn", "s2.fun_q2_cqn",
                "s2.rH_200", "s2.Ta_200", "s2.stand_density_basal_area",
                "s2.DBH_avg", "s2.d_gini", "s2.RA_forest",
                "s2.canopy_openness_2019", "s2.dom_tot_ratio") %>% 
  gdm::formatsitepair(bioData = ., bioFormat = 4, predData = .)


# Fungi
data_gdm_q0_fun <- data_gdm_q0 %>% 
  dplyr::rename(distance = fun_q0_cqn,
                s2.alg_q0_cqn = alg_q0_cqn,
                s2.bac_q0_cqn = bac_q0_cqn) %>% 
  tibble::add_column(., s1.alg_q0_cqn = 0,
                     s1.bac_q0_cqn = 0, 
                     weights = 1) %>% 
  dplyr::select("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                "s1.alg_q0_cqn", "s1.bac_q0_cqn",
                "s1.rH_200", "s1.Ta_200", "s1.stand_density_basal_area",
                "s1.DBH_avg", "s1.d_gini", "s1.RA_forest",
                "s1.canopy_openness_2019", "s1.dom_tot_ratio", 
                "s2.alg_q0_cqn", "s2.bac_q0_cqn",
                "s2.rH_200", "s2.Ta_200", "s2.stand_density_basal_area",
                "s2.DBH_avg", "s2.d_gini", "s2.RA_forest",
                "s2.canopy_openness_2019", "s2.dom_tot_ratio")%>% 
  gdm::formatsitepair(bioData = ., bioFormat = 4, predData = .)


data_gdm_q1_fun <- data_gdm_q1 %>% 
  dplyr::rename(distance = fun_q1_cqn,
                s2.alg_q1_cqn = alg_q1_cqn,
                s2.bac_q1_cqn = bac_q1_cqn) %>% 
  tibble::add_column(., s1.alg_q1_cqn = 0,
                     s1.bac_q1_cqn = 0, 
                     weights = 1) %>% 
  dplyr::select("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                "s1.alg_q1_cqn", "s1.bac_q1_cqn",
                "s1.rH_200", "s1.Ta_200", "s1.stand_density_basal_area",
                "s1.DBH_avg", "s1.d_gini", "s1.RA_forest",
                "s1.canopy_openness_2019", "s1.dom_tot_ratio", 
                "s2.alg_q1_cqn", "s2.bac_q1_cqn",
                "s2.rH_200", "s2.Ta_200", "s2.stand_density_basal_area",
                "s2.DBH_avg", "s2.d_gini", "s2.RA_forest",
                "s2.canopy_openness_2019", "s2.dom_tot_ratio") %>% 
  gdm::formatsitepair(bioData = ., bioFormat = 4, predData = .)

data_gdm_q2_fun <- data_gdm_q2 %>% 
  dplyr::rename(distance = fun_q2_cqn,
                s2.alg_q2_cqn = alg_q2_cqn,
                s2.bac_q2_cqn = bac_q2_cqn) %>% 
  tibble::add_column(., s1.alg_q2_cqn = 0,
                     s1.bac_q2_cqn = 0, 
                     weights = 1) %>% 
  dplyr::select("distance", "weights", "s1.xCoord", "s1.yCoord", "s2.xCoord", "s2.yCoord", 
                "s1.alg_q2_cqn", "s1.bac_q2_cqn",
                "s1.rH_200", "s1.Ta_200", "s1.stand_density_basal_area",
                "s1.DBH_avg", "s1.d_gini", "s1.RA_forest",
                "s1.canopy_openness_2019", "s1.dom_tot_ratio", 
                "s2.alg_q2_cqn", "s2.bac_q2_cqn",
                "s2.rH_200", "s2.Ta_200", "s2.stand_density_basal_area",
                "s2.DBH_avg", "s2.d_gini", "s2.RA_forest",
                "s2.canopy_openness_2019", "s2.dom_tot_ratio") %>% 
  gdm::formatsitepair(bioData = ., bioFormat = 4, predData = .)

##---------------------------------------------------------------
##                              GDMs                            -
##---------------------------------------------------------------

# q = 0
gdm_alg_q0 <- gdm::gdm(data_gdm_q0_alg, geo = T)
base::summary(gdm_alg_q0)
#base::plot(gdm_alg_q0)
gdm_alg_q0_perm <- gdm::gdm.varImp(data_gdm_q0_alg, geo = T, cores = 1, nPerm = 100)
gdm_alg_q0_perm

gdm_bac_q0 <- gdm::gdm(data_gdm_q0_bac, geo = T)
base::summary(gdm_bac_q0)
#base::plot(gdm_bac_q0)
gdm_bac_q0_perm <- gdm::gdm.varImp(data_gdm_q0_bac, geo = F, cores = 1, nPerm = 100)
gdm_bac_q0_perm

gdm_fun_q0 <- gdm::gdm(data_gdm_q0_fun, geo = T)
base::summary(gdm_fun_q0)
#base::plot(gdm_fun_q0)
gdm_fun_q0_perm <- gdm::gdm.varImp(data_gdm_q0_fun, geo = T, cores = 1, nPerm = 100)
gdm_fun_q0_perm

# q = 1
gdm_alg_q1 <- gdm::gdm(data_gdm_q1_alg, geo = T)
base::summary(gdm_alg_q1)
#base::plot(gdm_alg_q1)
gdm_alg_q1_perm <- gdm::gdm.varImp(data_gdm_q1_alg, geo = T, cores = 1, nPerm = 100)
gdm_alg_q1_perm

gdm_bac_q1 <- gdm::gdm(data_gdm_q1_bac, geo = T)
base::summary(gdm_bac_q1)
#base::plot(gdm_bac_q1)
gdm_bac_q1_perm <- gdm::gdm.varImp(data_gdm_q1_bac, geo = T, cores = 1, nPerm = 100)
gdm_bac_q1_perm

gdm_fun_q1 <- gdm::gdm(data_gdm_q1_fun, geo = T)
base::summary(gdm_fun_q1)
#base::plot(gdm_fun_q1)
gdm_fun_q1_perm <- gdm::gdm.varImp(data_gdm_q1_fun, geo = T, cores = 1, nPerm = 100)
gdm_fun_q1_perm

# q = 2
gdm_alg_q2 <- gdm::gdm(data_gdm_q2_alg, geo = T)
base::summary(gdm_alg_q2)
#base::plot(gdm_alg_q2)
gdm_alg_q2_perm <- gdm::gdm.varImp(data_gdm_q2_alg, geo = T, cores = 1, nPerm = 100)
gdm_alg_q2_perm

gdm_bac_q2 <- gdm::gdm(data_gdm_q2_bac, geo = T)
base::summary(gdm_bac_q2)
#base::plot(gdm_bac_q2)
gdm_bac_q2_perm <- gdm::gdm.varImp(data_gdm_q2_bac, geo = T, cores = 1, nPerm = 100)
gdm_bac_q2_perm

gdm_fun_q2 <- gdm::gdm(data_gdm_q2_fun, geo = T)
base::summary(gdm_fun_q2)
#base::plot(gdm_fun_q2)
gdm_fun_q2_perm <- gdm::gdm.varImp(data_gdm_q2_fun, geo = T, cores = 1, nPerm = 100)
gdm_fun_q2_perm

# Correct the p-values using the Benjamini-Hochberg procedure. 

p_gdm <- as.data.frame(gdm_alg_q0_perm$`Predictor p-values`)
p_gdm <- rbind(p_gdm, gdm_alg_q1_perm$`Predictor p-values`)
p_gdm <- rbind(p_gdm, gdm_alg_q2_perm$`Predictor p-values`)

p_gdm <- rbind(p_gdm, gdm_fun_q0_perm$`Predictor p-values`)
p_gdm <- rbind(p_gdm, gdm_fun_q1_perm$`Predictor p-values`)
p_gdm <- rbind(p_gdm, gdm_fun_q2_perm$`Predictor p-values`)

p_gdm <- rbind(p_gdm, gdm_bac_q0_perm$`Predictor p-values`)
p_gdm <- rbind(p_gdm, gdm_bac_q1_perm$`Predictor p-values`)
p_gdm <- rbind(p_gdm, gdm_bac_q2_perm$`Predictor p-values`)

p_adj_gdm <- stats::p.adjust(p_gdm$`All predictors`, method = "fdr")

p_vals_gdm <- data.frame(rownames(p_gdm) ,round(p_gdm$`All predictors`, 3), round(p_adj_gdm,3))

##---------------------------------------------------------------
##                    Variance partitioning                     -
##---------------------------------------------------------------

varSet <- vector("list", 2)
names(varSet) <- c("biotic", "abiotic")
varSet$abiotic <- c("rH_200", "Ta_200", "stand_density_basal_area",
                    "DBH_avg", "d_gini", "RA_forest", "canopy_openness_2019",
                    "dom_tot_ratio")

####
# Algae
####
#q0
varSet_alg_q0 <- vector("list", 2)
names(varSet_alg_q0) <- c("biotic", "abiotic")
varSet_alg_q0$biotic <- c("fun_q0_cqn", "bac_q0_cqn")
varSet_alg_q0$abiotic <- varSet$abiotic

alg_var_q0 <- gdm::gdm.partition.deviance(data_gdm_q0_alg, varSet_alg_q0, partSpace = T)

intersections <- alg_var_q0$VARIABLE_SET[8:15]

alg_var_q0 <- alg_var_q0 %>% 
  dplyr::filter(VARIABLE_SET %in% intersections) 

#q1
varSet_alg_q1 <- vector("list", 2)
names(varSet_alg_q1) <- c("biotic", "abiotic")
varSet_alg_q1$biotic <- c("fun_q1_cqn", "bac_q1_cqn")
varSet_alg_q1$abiotic <- varSet$abiotic

alg_var_q1 <- gdm::gdm.partition.deviance(data_gdm_q1_alg, varSet_alg_q1, partSpace = T) %>% 
  dplyr::filter(VARIABLE_SET %in% intersections)

#q2
varSet_alg_q2 <- vector("list", 2)
names(varSet_alg_q2) <- c("biotic", "abiotic")
varSet_alg_q2$biotic <- c("fun_q2_cqn", "bac_q2_cqn")
varSet_alg_q2$abiotic <- varSet$abiotic

alg_var_q2 <- gdm::gdm.partition.deviance(data_gdm_q2_alg, varSet_alg_q2, partSpace = T)%>% 
  dplyr::filter(VARIABLE_SET %in% intersections)

alg_variance <- rbind(alg_var_q0, alg_var_q1, alg_var_q2) %>% 
  cbind(q_lev = c(rep("q0", 8), rep("q1", 8), rep("q2", 8))) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "abiotic alone", "abiotic (a)")) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "biotic alone", "biotic (b)")) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "geo alone", "geographic (g)"))  %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "UNEXPLAINED", "unexplained"))  %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "biotic intersect abiotic, exclude geo", "b+a")) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "biotic intersect geo, exclude abiotic", "b+g")) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "geo intersect abiotic, exclude biotic", "a+g")) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "abiotic intersect biotic intersect geo", "a+b+g")) %>% 
  mutate(VARIABLE_SET = factor(VARIABLE_SET, levels = c("unexplained","a+b+g", "a+g", "b+g", "b+a",
                                                        "geographic (g)", "abiotic (a)", "biotic (b)")))  %>% 
  mutate(q_lev = factor(q_lev, levels = c("q2", "q1", "q0")))

# Fungi
####
#q0
varSet_fun_q0 <- vector("list", 2)
names(varSet_fun_q0) <- c("biotic", "abiotic")
varSet_fun_q0$biotic <- c("alg_q0_cqn", "bac_q0_cqn")
varSet_fun_q0$abiotic <- varSet$abiotic

fun_var_q0 <- gdm::gdm.partition.deviance(data_gdm_q0_fun, varSet_fun_q0, partSpace = T) %>% 
  dplyr::filter(VARIABLE_SET %in% intersections)

#q1
varSet_fun_q1 <- vector("list", 2)
names(varSet_fun_q1) <- c("biotic", "abiotic")
varSet_fun_q1$biotic <- c("alg_q1_cqn", "bac_q1_cqn")
varSet_fun_q1$abiotic <- varSet$abiotic

fun_var_q1 <- gdm::gdm.partition.deviance(data_gdm_q1_fun, varSet_fun_q1, partSpace = T) %>% 
  dplyr::filter(VARIABLE_SET %in% intersections)

#q2
varSet_fun_q2 <- vector("list", 2)
names(varSet_fun_q2) <- c("biotic", "abiotic")
varSet_fun_q2$biotic <- c("alg_q2_cqn", "bac_q2_cqn")
varSet_fun_q2$abiotic <- varSet$abiotic

fun_var_q2 <- gdm::gdm.partition.deviance(data_gdm_q2_fun, varSet_fun_q2, partSpace = T) %>% 
  dplyr::filter(VARIABLE_SET %in% intersections)

fun_variance <- rbind(fun_var_q0, fun_var_q1, fun_var_q2) %>% 
  cbind(q_lev = c(rep("q0", 8), rep("q1", 8), rep("q2", 8))) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "abiotic alone", "abiotic (a)")) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "biotic alone", "biotic (b)")) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "geo alone", "geographic (g)"))  %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "UNEXPLAINED", "unexplained"))  %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "biotic intersect abiotic, exclude geo", "b+a")) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "biotic intersect geo, exclude abiotic", "b+g")) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "geo intersect abiotic, exclude biotic", "a+g")) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "abiotic intersect biotic intersect geo", "a+b+g")) %>% 
  mutate(VARIABLE_SET = factor(VARIABLE_SET, levels = c("unexplained","a+b+g", "a+g", "b+g", "b+a",
                                                        "geographic (g)", "abiotic (a)", "biotic (b)")))  %>% 
  mutate(q_lev = factor(q_lev, levels = c("q2", "q1", "q0")))

####
# Bacteria
####
#q0
varSet_bac_q0 <- vector("list", 2)
names(varSet_bac_q0) <- c("biotic", "abiotic")
varSet_bac_q0$biotic <- c("fun_q0_cqn", "alg_q0_cqn")
varSet_bac_q0$abiotic <- varSet$abiotic

bac_var_q0 <- gdm::gdm.partition.deviance(data_gdm_q0_bac, varSet_bac_q0, partSpace = T) %>% 
  dplyr::filter(VARIABLE_SET %in% intersections)

#q1
varSet_bac_q1 <- vector("list", 2)
names(varSet_bac_q1) <- c("biotic", "abiotic")
varSet_bac_q1$biotic <- c("fun_q1_cqn", "alg_q1_cqn")
varSet_bac_q1$abiotic <- varSet$abiotic

bac_var_q1 <- gdm::gdm.partition.deviance(data_gdm_q1_bac, varSet_bac_q1, partSpace = T)  %>% 
  dplyr::filter(VARIABLE_SET %in% intersections)

#q2
varSet_bac_q2 <- vector("list", 2)
names(varSet_bac_q2) <- c("biotic", "abiotic")
varSet_bac_q2$biotic <- c("fun_q2_cqn", "alg_q2_cqn")
varSet_bac_q2$abiotic <- varSet$abiotic

bac_var_q2 <- gdm::gdm.partition.deviance(data_gdm_q2_bac, varSet_bac_q2, partSpace = T) %>% 
  dplyr::filter(VARIABLE_SET %in% intersections)

bac_variance <- rbind(bac_var_q0, bac_var_q1, bac_var_q2) %>% 
  cbind(q_lev = c(rep("q0", 8), rep("q1", 8), rep("q2", 8))) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "abiotic alone", "abiotic (a)")) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "biotic alone", "biotic (b)")) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "geo alone", "geographic (g)"))  %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "UNEXPLAINED", "unexplained"))  %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "biotic intersect abiotic, exclude geo", "b+a")) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "biotic intersect geo, exclude abiotic", "b+g")) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "geo intersect abiotic, exclude biotic", "a+g")) %>% 
  mutate(VARIABLE_SET = replace(VARIABLE_SET, VARIABLE_SET == "abiotic intersect biotic intersect geo", "a+b+g")) %>% 
  mutate(VARIABLE_SET = factor(VARIABLE_SET, levels = c("unexplained","a+b+g", "a+g", "b+g", "b+a",
                                                        "geographic (g)", "abiotic (a)", "biotic (b)")))  %>% 
  mutate(q_lev = factor(q_lev, levels = c("q2", "q1", "q0")))


##---------------------------------------------------------------
##   Linear Models of beta diversity to get effect direction    -
##---------------------------------------------------------------

# First we need to calculate the absolute differences in our explanatory variables.

data_lm_beta_q0 <- data_gdm_q0 %>% 
  dplyr::mutate(., delta_rH_200 = abs(s1.rH_200 - s2.rH_200),
                delta_Ta_200 = abs(s1.Ta_200 - s2.Ta_200),
                delta_DBH_avg = abs(s1.DBH_avg - s2.DBH_avg),
                delta_stand_density_basal_area = abs(s1.stand_density_basal_area - s2.stand_density_basal_area),
                delta_d_gini = abs(s1.d_gini - s2.d_gini),
                delta_RA_forest = abs(s1.RA_forest - s2.RA_forest),
                delta_canopy_openness_2019 = abs(s1.canopy_openness_2019 - s2.canopy_openness_2019),
                delta_dom_tot_ratio = abs(s1.dom_tot_ratio - s2.dom_tot_ratio)) %>% 
  dplyr::select(., "site1", "site2", "alg_q0_cqn", "bac_q0_cqn",
                "fun_q0_cqn", all_of(contains("delta")))

data_lm_beta_q1 <- data_gdm_q1 %>% 
  dplyr::mutate(., delta_rH_200 = abs(s1.rH_200 - s2.rH_200),
                delta_Ta_200 = abs(s1.Ta_200 - s2.Ta_200),
                delta_DBH_avg = abs(s1.DBH_avg - s2.DBH_avg),
                delta_stand_density_basal_area = abs(s1.stand_density_basal_area - s2.stand_density_basal_area),
                delta_d_gini = abs(s1.d_gini - s2.d_gini),
                delta_RA_forest = abs(s1.RA_forest - s2.RA_forest),
                delta_canopy_openness_2019 = abs(s1.canopy_openness_2019 - s2.canopy_openness_2019),
                delta_dom_tot_ratio = abs(s1.dom_tot_ratio - s2.dom_tot_ratio)) %>% 
  dplyr::select(., "site1", "site2", "alg_q1_cqn", "bac_q1_cqn",
                "fun_q1_cqn", all_of(contains("delta")))

data_lm_beta_q2 <- data_gdm_q2 %>% 
  dplyr::mutate(., delta_rH_200 = abs(s1.rH_200 - s2.rH_200),
                delta_Ta_200 = abs(s1.Ta_200 - s2.Ta_200),
                delta_DBH_avg = abs(s1.DBH_avg - s2.DBH_avg),
                delta_stand_density_basal_area = abs(s1.stand_density_basal_area - s2.stand_density_basal_area),
                delta_d_gini = abs(s1.d_gini - s2.d_gini),
                delta_RA_forest = abs(s1.RA_forest - s2.RA_forest),
                delta_canopy_openness_2019 = abs(s1.canopy_openness_2019 - s2.canopy_openness_2019),
                delta_dom_tot_ratio = abs(s1.dom_tot_ratio - s2.dom_tot_ratio)) %>% 
  dplyr::select(., "site1", "site2", "alg_q2_cqn", "bac_q2_cqn",
                "fun_q2_cqn", all_of(contains("delta")))


# Models for q = 0
data_lm_beta_q0_scaled <- data_lm_beta_q0  %>% 
  # Scale the explanatory variables since they have very different units. 
  dplyr::mutate(across(.cols = c("delta_rH_200", "delta_Ta_200",
                                 "delta_stand_density_basal_area", 
                                 "delta_DBH_avg", "delta_d_gini",
                                 "delta_RA_forest", 
                                 "delta_canopy_openness_2019",
                                 "delta_dom_tot_ratio"), scale))

lm_beta_q0_alg <- stats::lm(alg_q0_cqn ~ bac_q0_cqn +
                              fun_q0_cqn + 
                              delta_rH_200 +
                              delta_Ta_200 +
                              delta_DBH_avg +
                              delta_stand_density_basal_area +
                              delta_d_gini +
                              delta_RA_forest +
                              delta_canopy_openness_2019 +
                              delta_dom_tot_ratio,
                              data = data_lm_beta_q0_scaled)

summary(lm_beta_q0_alg)
plot(effects::allEffects(lm_beta_q0_alg))

lm_beta_q0_bac <- stats::lm(bac_q0_cqn ~ alg_q0_cqn +
                              fun_q0_cqn + 
                              delta_rH_200 +
                              delta_Ta_200 +
                              delta_DBH_avg +
                              delta_stand_density_basal_area +
                              delta_d_gini +
                              delta_RA_forest +
                              delta_canopy_openness_2019 +
                              delta_dom_tot_ratio,
                            data = data_lm_beta_q0_scaled)

summary(lm_beta_q0_bac)
plot(effects::allEffects(lm_beta_q0_bac))

lm_beta_q0_fun <- stats::lm(fun_q0_cqn ~ alg_q0_cqn + 
                              bac_q0_cqn + 
                              delta_rH_200 +
                              delta_Ta_200 +
                              delta_DBH_avg +
                              delta_stand_density_basal_area +
                              delta_d_gini +
                              delta_RA_forest +
                              delta_canopy_openness_2019 +
                              delta_dom_tot_ratio,
                            data = data_lm_beta_q0_scaled)

summary(lm_beta_q0_fun)
plot(effects::allEffects(lm_beta_q0_fun))


# Models for q = 1

data_lm_beta_q1_scaled <- data_lm_beta_q1  %>% 
  # Scale the explanatory variables since they have very different units. 
  dplyr::mutate(across(.cols = c("delta_rH_200", "delta_Ta_200",
                                 "delta_stand_density_basal_area", 
                                 "delta_DBH_avg", "delta_d_gini",
                                 "delta_RA_forest", 
                                 "delta_canopy_openness_2019",
                                 "delta_dom_tot_ratio"), scale))

lm_beta_q1_alg <- stats::lm(alg_q1_cqn ~ bac_q1_cqn +
                              fun_q1_cqn + 
                              delta_rH_200 +
                              delta_Ta_200 +
                              delta_DBH_avg +
                              delta_stand_density_basal_area +
                              delta_d_gini +
                              delta_RA_forest +
                              delta_canopy_openness_2019 +
                              delta_dom_tot_ratio,
                            data = data_lm_beta_q1_scaled)

summary(lm_beta_q1_alg)
plot(effects::allEffects(lm_beta_q1_alg))

lm_beta_q1_bac <- stats::lm(bac_q1_cqn ~ alg_q1_cqn +
                              fun_q1_cqn + 
                              delta_rH_200 +
                              delta_Ta_200 +
                              delta_DBH_avg +
                              delta_stand_density_basal_area +
                              delta_d_gini +
                              delta_RA_forest +
                              delta_canopy_openness_2019 +
                              delta_dom_tot_ratio,
                            data = data_lm_beta_q1_scaled)

summary(lm_beta_q1_bac)
plot(effects::allEffects(lm_beta_q1_bac))

lm_beta_q1_fun <- stats::lm(fun_q1_cqn ~ alg_q1_cqn + 
                              bac_q1_cqn + 
                              delta_rH_200 +
                              delta_Ta_200 +
                              delta_DBH_avg +
                              delta_stand_density_basal_area +
                              delta_d_gini +
                              delta_RA_forest +
                              delta_canopy_openness_2019 +
                              delta_dom_tot_ratio,
                            data = data_lm_beta_q1_scaled)

summary(lm_beta_q1_fun)
plot(effects::allEffects(lm_beta_q1_fun))

# Models q = 2

data_lm_beta_q2_scaled <- data_lm_beta_q2  %>% 
  # Scale the explanatory variables since they have very different units. 
  dplyr::mutate(across(.cols = c("delta_rH_200", "delta_Ta_200",
                                 "delta_stand_density_basal_area", 
                                 "delta_DBH_avg", "delta_d_gini",
                                 "delta_RA_forest", 
                                 "delta_canopy_openness_2019",
                                 "delta_dom_tot_ratio"), scale))

lm_beta_q2_alg <- stats::lm(alg_q2_cqn ~ bac_q2_cqn +
                              fun_q2_cqn + 
                              delta_rH_200 +
                              delta_Ta_200 +
                              delta_DBH_avg +
                              delta_stand_density_basal_area +
                              delta_d_gini +
                              delta_RA_forest +
                              delta_canopy_openness_2019 +
                              delta_dom_tot_ratio,
                            data = data_lm_beta_q2_scaled)

summary(lm_beta_q2_alg)
plot(effects::allEffects(lm_beta_q2_alg))

lm_beta_q2_bac <- stats::lm(bac_q2_cqn ~ alg_q2_cqn +
                              fun_q2_cqn + 
                              delta_rH_200 +
                              delta_Ta_200 +
                              delta_DBH_avg +
                              delta_stand_density_basal_area +
                              delta_d_gini +
                              delta_RA_forest +
                              delta_canopy_openness_2019 +
                              delta_dom_tot_ratio,
                            data = data_lm_beta_q2_scaled)

summary(lm_beta_q2_bac)
plot(effects::allEffects(lm_beta_q2_bac))

lm_beta_q2_fun <- stats::lm(fun_q2_cqn ~ alg_q2_cqn + 
                              bac_q2_cqn + 
                              delta_rH_200 +
                              delta_Ta_200 +
                              delta_DBH_avg +
                              delta_stand_density_basal_area +
                              delta_d_gini +
                              delta_RA_forest +
                              delta_canopy_openness_2019 +
                              delta_dom_tot_ratio,
                            data = data_lm_beta_q2_scaled)

summary(lm_beta_q2_fun)
plot(effects::allEffects(lm_beta_q2_fun))

##---------------------------------------------------------------
##                Plotting the beta diversity drivers           -
##---------------------------------------------------------------

# Colors

algae_col <- "#49BEAA"
bacteria_col <- "#e2ca20"
fungi_col <- "#ea594e"

# Extract the splines
gdm_alg_q0_splineDat <- gdm::isplineExtract(gdm_alg_q0)
gdm_bac_q0_splineDat <- gdm::isplineExtract(gdm_bac_q0)
gdm_fun_q0_splineDat <- gdm::isplineExtract(gdm_fun_q0)

gdm_alg_q1_splineDat <- gdm::isplineExtract(gdm_alg_q1)
gdm_bac_q1_splineDat <- gdm::isplineExtract(gdm_bac_q1)
gdm_fun_q1_splineDat <- gdm::isplineExtract(gdm_fun_q1)

gdm_alg_q2_splineDat <- gdm::isplineExtract(gdm_alg_q2)
gdm_bac_q2_splineDat <- gdm::isplineExtract(gdm_bac_q2)
gdm_fun_q2_splineDat <- gdm::isplineExtract(gdm_fun_q2)

# Canopy Openness
gdm_alg_q0_plot_canopy <- data.frame(x_canopy = gdm_alg_q0_splineDat$x[,"canopy_openness_2019"], 
                                     y_canopy = gdm_alg_q0_splineDat$y[,"canopy_openness_2019"])
gdm_alg_q1_plot_canopy <- data.frame(x_canopy = gdm_alg_q1_splineDat$x[,"canopy_openness_2019"], 
                                     y_canopy = gdm_alg_q1_splineDat$y[,"canopy_openness_2019"])
gdm_alg_q2_plot_canopy <- data.frame(x_canopy = gdm_alg_q2_splineDat$x[,"canopy_openness_2019"], 
                                     y_canopy = gdm_alg_q2_splineDat$y[,"canopy_openness_2019"])

gdm_bac_q0_plot_canopy <- data.frame(x_canopy = gdm_bac_q0_splineDat$x[,"canopy_openness_2019"],
                                     y_canopy = gdm_bac_q0_splineDat$y[,"canopy_openness_2019"])
gdm_bac_q1_plot_canopy <- data.frame(x_canopy = gdm_bac_q1_splineDat$x[,"canopy_openness_2019"],
                                     y_canopy = gdm_bac_q1_splineDat$y[,"canopy_openness_2019"])
gdm_bac_q2_plot_canopy <- data.frame(x_canopy = gdm_bac_q2_splineDat$x[,"canopy_openness_2019"],
                                     y_canopy = gdm_bac_q2_splineDat$y[,"canopy_openness_2019"])

gdm_fun_q0_plot_canopy <- data.frame(x_canopy = gdm_fun_q0_splineDat$x[,"canopy_openness_2019"],
                                     y_canopy = gdm_fun_q0_splineDat$y[,"canopy_openness_2019"])
gdm_fun_q1_plot_canopy <- data.frame(x_canopy = gdm_fun_q1_splineDat$x[,"canopy_openness_2019"],
                                     y_canopy = gdm_fun_q1_splineDat$y[,"canopy_openness_2019"])
gdm_fun_q2_plot_canopy <- data.frame(x_canopy = gdm_fun_q2_splineDat$x[,"canopy_openness_2019"],
                                     y_canopy = gdm_fun_q2_splineDat$y[,"canopy_openness_2019"])

######
# Effect on algal communities
######
algae_canopy_effect <- ggplot() +
  geom_line(data = gdm_alg_q0_plot_canopy, aes(x = x_canopy, y = y_canopy),
            color = algae_col, linewidth = 0.7) + 
  geom_line(data = gdm_alg_q1_plot_canopy, aes(x = x_canopy, y = y_canopy),
            color = algae_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_line(data = gdm_alg_q2_plot_canopy, aes(x = x_canopy, y = y_canopy),
            color = algae_col, linewidth = 0.7, linetype = "dotted") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
    labs(x= "Canopy openness (%)",    
         y= "Algal \u03B2-diversity")+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2))
algae_canopy_effect

######
# Effect on fungal communities
######
fungi_canopy_effect <- ggplot() +
  geom_line(data = gdm_fun_q0_plot_canopy, aes(x = x_canopy, y = y_canopy),
            color = fungi_col, linewidth = 0.7) + 
  geom_text(data = gdm_fun_q0_plot_canopy, aes(x = max(x_canopy) + 1,
                                        y = max(y_canopy)), label = "*",
            size = 4)+ 
  geom_line(data = gdm_fun_q1_plot_canopy, aes(x = x_canopy, y = y_canopy),
            color = fungi_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_text(data = gdm_fun_q1_plot_canopy, aes(x = max(x_canopy) + 1,
                                        y = max(y_canopy)), label = "*",
            size = 4)+ 
  geom_line(data = gdm_fun_q2_plot_canopy, aes(x = x_canopy, y = y_canopy),
            color = fungi_col, linewidth = 0.7, linetype = "dotted") + 
  geom_text(data = gdm_fun_q2_plot_canopy, aes(x = max(x_canopy) + 1,
                                        y = max(y_canopy)), label = "*",
            size = 4)+ 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(x= "Canopy openness (%)",    
       y= "Fungal \u03B2-diversity")+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2))
fungi_canopy_effect

######
# Effect on bacterial communities
######
bacteria_canopy_effect <- ggplot() +
  geom_line(data = gdm_bac_q0_plot_canopy, aes(x = x_canopy, y = y_canopy),
            color = bacteria_col, linewidth = 0.7) + 
  geom_text(data = gdm_bac_q0_plot_canopy, aes(x = max(x_canopy) + 1,
                                               y = max(y_canopy)), label = "*",
            size = 4) + 
  geom_line(data = gdm_bac_q1_plot_canopy, aes(x = x_canopy, y = y_canopy),
            color = bacteria_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_text(data = gdm_bac_q1_plot_canopy, aes(x = max(x_canopy) + 1,
                                               y = max(y_canopy)), label = "*",
            size = 4) + 
  geom_line(data = gdm_bac_q2_plot_canopy, aes(x = x_canopy, y = y_canopy),
            color = bacteria_col, linewidth = 0.7, linetype = "dotted") + 
  geom_text(data = gdm_bac_q2_plot_canopy, aes(x = max(x_canopy) + 1,
                                               y = max(y_canopy)), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(x= "Canopy openness (%)",    
       y= "Bacterial \u03B2-diversity")+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2))
bacteria_canopy_effect


# Relative humidity
gdm_alg_q0_plot_humidity <- data.frame(x_humidity = gdm_alg_q0_splineDat$x[,"rH_200"], 
                                     y_humidity = gdm_alg_q0_splineDat$y[,"rH_200"])
gdm_bac_q0_plot_humidity <- data.frame(x_humidity = gdm_bac_q0_splineDat$x[,"rH_200"],
                                     y_humidity = gdm_bac_q0_splineDat$y[,"rH_200"])
gdm_fun_q0_plot_humidity <- data.frame(x_humidity = gdm_fun_q0_splineDat$x[,"rH_200"],
                                     y_humidity = gdm_fun_q0_splineDat$y[,"rH_200"])

gdm_alg_q1_plot_humidity <- data.frame(x_humidity = gdm_alg_q1_splineDat$x[,"rH_200"], 
                                     y_humidity = gdm_alg_q1_splineDat$y[,"rH_200"])
gdm_bac_q1_plot_humidity <- data.frame(x_humidity = gdm_bac_q1_splineDat$x[,"rH_200"],
                                     y_humidity = gdm_bac_q1_splineDat$y[,"rH_200"])
gdm_fun_q1_plot_humidity <- data.frame(x_humidity = gdm_fun_q1_splineDat$x[,"rH_200"],
                                     y_humidity = gdm_fun_q1_splineDat$y[,"rH_200"])


gdm_alg_q2_plot_humidity <- data.frame(x_humidity = gdm_alg_q2_splineDat$x[,"rH_200"], 
                                     y_humidity = gdm_alg_q2_splineDat$y[,"rH_200"])
gdm_bac_q2_plot_humidity <- data.frame(x_humidity = gdm_bac_q2_splineDat$x[,"rH_200"],
                                     y_humidity = gdm_bac_q2_splineDat$y[,"rH_200"])
gdm_fun_q2_plot_humidity <- data.frame(x_humidity = gdm_fun_q2_splineDat$x[,"rH_200"],
                                     y_humidity = gdm_fun_q2_splineDat$y[,"rH_200"])

######
# Effect on algal communities
######
algae_humidity_effect <- ggplot() +
  geom_line(data = gdm_alg_q0_plot_humidity, aes(x = x_humidity, y = y_humidity),
            color = algae_col, linewidth = 0.7) + 
  geom_text(data = gdm_alg_q0_plot_humidity, aes(x = max(x_humidity) + 0.5,
                                               y = max(y_humidity)), label = "*",
            size = 4) +
  geom_line(data = gdm_alg_q1_plot_humidity, aes(x = x_humidity, y = y_humidity),
            color = algae_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_line(data = gdm_alg_q2_plot_humidity, aes(x = x_humidity, y = y_humidity),
            color = algae_col, linewidth = 0.7, linetype = "dotted") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Relative Humidity",    
       y= "Algal \u03B2-diversity")
algae_humidity_effect

######
# Effect on fungal communities
######
fungi_humidity_effect <- ggplot() +
  geom_line(data = gdm_fun_q0_plot_humidity, aes(x = x_humidity, y = y_humidity),
            color = fungi_col, linewidth = 0.7) + 
  geom_line(data = gdm_fun_q1_plot_humidity, aes(x = x_humidity, y = y_humidity),
            color = fungi_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_text(data = gdm_fun_q1_plot_humidity, aes(x = max(x_humidity) + 0.5,
                                                 y = max(y_humidity)), label = "*",
            size = 4) + 
  geom_line(data = gdm_fun_q2_plot_humidity, aes(x = x_humidity, y = y_humidity),
            color = fungi_col, linewidth = 0.7, linetype = "dotted") + 
  geom_text(data = gdm_fun_q2_plot_humidity, aes(x = max(x_humidity) + 0.5,
                                                 y = max(y_humidity)), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Relative Humidity",    
       y= "Fungal \u03B2-diversity")
fungi_humidity_effect

######
# Effect on bacterial communities
######
bacteria_humidity_effect <- ggplot() +
  geom_line(data = gdm_bac_q0_plot_humidity, aes(x = x_humidity, y = y_humidity),
            color = bacteria_col, linewidth = 0.7) + 
  geom_line(data = gdm_bac_q1_plot_humidity, aes(x = x_humidity, y = y_humidity),
            color = bacteria_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_line(data = gdm_bac_q2_plot_humidity, aes(x = x_humidity, y = y_humidity),
            color = bacteria_col, linewidth = 0.7, linetype = "dotted") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Relative Humidity",    
       y= "Bacterial \u03B2-diversity")
bacteria_humidity_effect

# Temperature
gdm_alg_q0_plot_temperature <- data.frame(x_temperature = gdm_alg_q0_splineDat$x[,"Ta_200"], 
                                          y_temperature = gdm_alg_q0_splineDat$y[,"Ta_200"])
gdm_bac_q0_plot_temperature <- data.frame(x_temperature = gdm_bac_q0_splineDat$x[,"Ta_200"],
                                          y_temperature = gdm_bac_q0_splineDat$y[,"Ta_200"])
gdm_fun_q0_plot_temperature <- data.frame(x_temperature = gdm_fun_q0_splineDat$x[,"Ta_200"],
                                          y_temperature = gdm_fun_q0_splineDat$y[,"Ta_200"])

gdm_alg_q1_plot_temperature <- data.frame(x_temperature = gdm_alg_q1_splineDat$x[,"Ta_200"], 
                                          y_temperature = gdm_alg_q1_splineDat$y[,"Ta_200"])
gdm_bac_q1_plot_temperature <- data.frame(x_temperature = gdm_bac_q1_splineDat$x[,"Ta_200"],
                                          y_temperature = gdm_bac_q1_splineDat$y[,"Ta_200"])
gdm_fun_q1_plot_temperature <- data.frame(x_temperature = gdm_fun_q1_splineDat$x[,"Ta_200"],
                                          y_temperature = gdm_fun_q1_splineDat$y[,"Ta_200"])


gdm_alg_q2_plot_temperature <- data.frame(x_temperature = gdm_alg_q2_splineDat$x[,"Ta_200"], 
                                          y_temperature = gdm_alg_q2_splineDat$y[,"Ta_200"])
gdm_bac_q2_plot_temperature <- data.frame(x_temperature = gdm_bac_q2_splineDat$x[,"Ta_200"],
                                          y_temperature = gdm_bac_q2_splineDat$y[,"Ta_200"])
gdm_fun_q2_plot_temperature <- data.frame(x_temperature = gdm_fun_q2_splineDat$x[,"Ta_200"],
                                          y_temperature = gdm_fun_q2_splineDat$y[,"Ta_200"])


######
# Effect on algal communities
######
algae_temperature_effect <- ggplot() +
  geom_line(data = gdm_alg_q0_plot_temperature, aes(x = x_temperature, y = y_temperature),
            color = algae_col, linewidth = 0.7)  +
  geom_line(data = gdm_alg_q1_plot_temperature, aes(x = x_temperature, y = y_temperature),
            color = algae_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_line(data = gdm_alg_q2_plot_temperature, aes(x = x_temperature, y = y_temperature),
            color = algae_col, linewidth = 0.7, linetype = "dotted") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(x= "Temperature",    
       y= "Algal \u03B2-diversity")+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2))
algae_temperature_effect

######
# Effect on fungal communities
######
fungi_temperature_effect <- ggplot() +
  geom_line(data = gdm_fun_q0_plot_temperature, aes(x = x_temperature, y = y_temperature),
            color = fungi_col, linewidth = 0.7) + 
  geom_text(data = gdm_fun_q0_plot_temperature, aes(x = max(x_temperature) + 0.1,
                                                 y = max(y_temperature)), label = "*",
            size = 4) + 
  geom_line(data = gdm_fun_q1_plot_temperature, aes(x = x_temperature, y = y_temperature),
            color = fungi_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_text(data = gdm_fun_q1_plot_temperature, aes(x = max(x_temperature) + 0.1,
                                                 y = max(y_temperature)), label = "*",
            size = 4) + 
  geom_line(data = gdm_fun_q2_plot_temperature, aes(x = x_temperature, y = y_temperature),
            color = fungi_col, linewidth = 0.7, linetype = "dotted") + 
  geom_text(data = gdm_fun_q2_plot_temperature, aes(x = max(x_temperature) + 0.1,
                                                 y = max(y_temperature)), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(x= "Temperature",    
       y= "Fungal \u03B2-diversity")+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2))
fungi_temperature_effect

######
# Effect on bacterial communities
######
bacteria_temperature_effect <- ggplot() +
  geom_line(data = gdm_bac_q0_plot_temperature, aes(x = x_temperature, y = y_temperature),
            color = bacteria_col, linewidth = 0.7) + 
  geom_line(data = gdm_bac_q1_plot_temperature, aes(x = x_temperature, y = y_temperature),
            color = bacteria_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_text(data = gdm_bac_q1_plot_temperature, aes(x = max(x_temperature) + 0.1,
                                                    y = max(y_temperature)), label = "*",
            size = 4) + 
  geom_line(data = gdm_bac_q2_plot_temperature, aes(x = x_temperature, y = y_temperature),
            color = bacteria_col, linewidth = 0.7, linetype = "dotted") + 
  geom_text(data = gdm_bac_q2_plot_temperature, aes(x = max(x_temperature) + 0.1,
                                                    y = max(y_temperature)), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(x= "Temperature",    
       y= "Bacterial \u03B2-diversity")+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2))
bacteria_temperature_effect

# Gini coefficient
gdm_alg_q0_plot_gini <- data.frame(x_gini = gdm_alg_q0_splineDat$x[,"d_gini"], 
                                       y_gini = gdm_alg_q0_splineDat$y[,"d_gini"])
gdm_bac_q0_plot_gini <- data.frame(x_gini = gdm_bac_q0_splineDat$x[,"d_gini"],
                                       y_gini = gdm_bac_q0_splineDat$y[,"d_gini"])
gdm_fun_q0_plot_gini <- data.frame(x_gini = gdm_fun_q0_splineDat$x[,"d_gini"],
                                       y_gini = gdm_fun_q0_splineDat$y[,"d_gini"])

gdm_alg_q1_plot_gini <- data.frame(x_gini = gdm_alg_q1_splineDat$x[,"d_gini"], 
                                       y_gini = gdm_alg_q1_splineDat$y[,"d_gini"])
gdm_bac_q1_plot_gini <- data.frame(x_gini = gdm_bac_q1_splineDat$x[,"d_gini"],
                                       y_gini = gdm_bac_q1_splineDat$y[,"d_gini"])
gdm_fun_q1_plot_gini <- data.frame(x_gini = gdm_fun_q1_splineDat$x[,"d_gini"],
                                       y_gini = gdm_fun_q1_splineDat$y[,"d_gini"])


gdm_alg_q2_plot_gini <- data.frame(x_gini = gdm_alg_q2_splineDat$x[,"d_gini"], 
                                       y_gini = gdm_alg_q2_splineDat$y[,"d_gini"])
gdm_bac_q2_plot_gini <- data.frame(x_gini = gdm_bac_q2_splineDat$x[,"d_gini"],
                                       y_gini = gdm_bac_q2_splineDat$y[,"d_gini"])
gdm_fun_q2_plot_gini <- data.frame(x_gini = gdm_fun_q2_splineDat$x[,"d_gini"],
                                       y_gini = gdm_fun_q2_splineDat$y[,"d_gini"])


######
# Effect on algal communities
######
algae_gini_effect <- ggplot() +
  geom_line(data = gdm_alg_q0_plot_gini, aes(x = x_gini, y = y_gini),
            color = algae_col, linewidth = 0.7)  +
  geom_line(data = gdm_alg_q1_plot_gini, aes(x = x_gini, y = y_gini),
            color = algae_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_line(data = gdm_alg_q2_plot_gini, aes(x = x_gini, y = y_gini),
            color = algae_col, linewidth = 0.7, linetype = "dotted") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Gini Coefficient",    
       y= "Algal \u03B2-diversity")
algae_gini_effect

######
# Effect on fungal communities
######
fungi_gini_effect <- ggplot() +
  geom_line(data = gdm_fun_q0_plot_gini, aes(x = x_gini, y = y_gini),
            color = fungi_col, linewidth = 0.7) + 
  geom_line(data = gdm_fun_q1_plot_gini, aes(x = x_gini, y = y_gini),
            color = fungi_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_line(data = gdm_fun_q2_plot_gini, aes(x = x_gini, y = y_gini),
            color = fungi_col, linewidth = 0.7, linetype = "dotted") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Gini Coefficient",    
       y= "Fungal \u03B2-diversity")
fungi_gini_effect

######
# Effect on bacterial communities
######
bacteria_gini_effect <- ggplot() +
  geom_line(data = gdm_bac_q0_plot_gini, aes(x = x_gini, y = y_gini),
            color = bacteria_col, linewidth = 0.7) + 
  geom_text(data = gdm_bac_q0_plot_gini, aes(x = max(x_gini) + 0.01,
                                                    y = max(y_gini)), label = "*",
            size = 4) + 
  geom_line(data = gdm_bac_q1_plot_gini, aes(x = x_gini, y = y_gini),
            color = bacteria_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_text(data = gdm_bac_q1_plot_gini, aes(x = max(x_gini) + 0.01,
                                                    y = max(y_gini)), label = "*",
            size = 4) + 
  geom_line(data = gdm_bac_q2_plot_gini, aes(x = x_gini, y = y_gini),
            color = bacteria_col, linewidth = 0.7, linetype = "dotted") + 
  geom_text(data = gdm_bac_q2_plot_gini, aes(x = max(x_gini) + 0.01,
                                                    y = max(y_gini)), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Gini Coefficient",    
       y= "Bacterial \u03B2-diversity")
bacteria_gini_effect

# Geographic Distance
gdm_alg_q0_plot_geo <- data.frame(x_geo = gdm_alg_q0_splineDat$x[,"Geographic"], 
                                   y_geo = gdm_alg_q0_splineDat$y[,"Geographic"])
gdm_bac_q0_plot_geo <- data.frame(x_geo = gdm_bac_q0_splineDat$x[,"Geographic"],
                                   y_geo = gdm_bac_q0_splineDat$y[,"Geographic"])
gdm_fun_q0_plot_geo <- data.frame(x_geo = gdm_fun_q0_splineDat$x[,"Geographic"],
                                   y_geo = gdm_fun_q0_splineDat$y[,"Geographic"])

gdm_alg_q1_plot_geo <- data.frame(x_geo = gdm_alg_q1_splineDat$x[,"Geographic"], 
                                  y_geo = gdm_alg_q1_splineDat$y[,"Geographic"])
gdm_bac_q1_plot_geo <- data.frame(x_geo = gdm_bac_q1_splineDat$x[,"Geographic"],
                                  y_geo = gdm_bac_q1_splineDat$y[,"Geographic"])
gdm_fun_q1_plot_geo <- data.frame(x_geo = gdm_fun_q1_splineDat$x[,"Geographic"],
                                  y_geo = gdm_fun_q1_splineDat$y[,"Geographic"])


gdm_alg_q2_plot_geo <- data.frame(x_geo = gdm_alg_q2_splineDat$x[,"Geographic"], 
                                  y_geo = gdm_alg_q2_splineDat$y[,"Geographic"])
gdm_bac_q2_plot_geo <- data.frame(x_geo = gdm_bac_q2_splineDat$x[,"Geographic"],
                                  y_geo = gdm_bac_q2_splineDat$y[,"Geographic"])
gdm_fun_q2_plot_geo <- data.frame(x_geo = gdm_fun_q2_splineDat$x[,"Geographic"],
                                  y_geo = gdm_fun_q2_splineDat$y[,"Geographic"])

######
# Effect on algal communities
######
algae_geo_effect <- ggplot() +
  geom_line(data = gdm_alg_q0_plot_geo, aes(x = x_geo, y = y_geo),
            color = algae_col, linewidth = 0.7) + 
  geom_text(data = gdm_alg_q0_plot_geo, aes(x = max(x_geo) + 0.2,
                                             y = max(y_geo)), label = "*",
            size = 4)  +
  geom_line(data = gdm_alg_q1_plot_geo, aes(x = x_geo, y = y_geo),
            color = algae_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_text(data = gdm_alg_q1_plot_geo, aes(x = max(x_geo) + 0.2,
                                             y = max(y_geo)), label = "*",
            size = 4)  + 
  geom_line(data = gdm_alg_q2_plot_geo, aes(x = x_geo, y = y_geo),
            color = algae_col, linewidth = 0.7, linetype = "dotted") + 
  geom_text(data = gdm_alg_q2_plot_geo, aes(x = max(x_geo) + 0.2,
                                             y = max(y_geo)), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(x= "Geographic Distance",    
       y= "Algal \u03B2-diversity")
algae_geo_effect

######
# Effect on fungal communities
######
fungi_geo_effect <- ggplot() +
  geom_line(data = gdm_fun_q0_plot_geo, aes(x = x_geo, y = y_geo),
            color = fungi_col, linewidth = 0.7) + 
  geom_text(data = gdm_fun_q0_plot_geo, aes(x = max(x_geo) + 0.2,
                                            y = max(y_geo)), label = "*",
            size = 4) + 
  geom_line(data = gdm_fun_q1_plot_geo, aes(x = x_geo, y = y_geo),
            color = fungi_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_text(data = gdm_fun_q1_plot_geo, aes(x = max(x_geo) + 0.2,
                                            y = max(y_geo)), label = "*",
            size = 4) + 
  geom_line(data = gdm_fun_q2_plot_geo, aes(x = x_geo, y = y_geo),
            color = fungi_col, linewidth = 0.7, linetype = "dotted") + 
  geom_text(data = gdm_fun_q2_plot_geo, aes(x = max(x_geo) + 0.2,
                                            y = max(y_geo)), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(x= "Geographic Distance",    
       y= "Fungal \u03B2-diversity")
fungi_geo_effect

######
# Effect on bacterial communities
######
bacteria_geo_effect <- ggplot() +
  geom_line(data = gdm_bac_q0_plot_geo, aes(x = x_geo, y = y_geo),
            color = bacteria_col, linewidth = 0.7) + 
  geom_line(data = gdm_bac_q1_plot_geo, aes(x = x_geo, y = y_geo),
            color = bacteria_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_line(data = gdm_bac_q2_plot_geo, aes(x = x_geo, y = y_geo),
            color = bacteria_col, linewidth = 0.7, linetype = "dotted") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(x= "Geographic Distance",    
       y= "Bacterial \u03B2-diversity")
bacteria_geo_effect

# Fungal beta effect
gdm_alg_q0_plot_fungi <- data.frame(x_fungi = gdm_alg_q0_splineDat$x[,"fun_q0_cqn"], 
                                    y_fungi = gdm_alg_q0_splineDat$y[,"fun_q0_cqn"])
gdm_bac_q0_plot_fungi <- data.frame(x_fungi = gdm_bac_q0_splineDat$x[,"fun_q0_cqn"],
                                    y_fungi = gdm_bac_q0_splineDat$y[,"fun_q0_cqn"])

gdm_alg_q1_plot_fungi <- data.frame(x_fungi = gdm_alg_q1_splineDat$x[,"fun_q1_cqn"], 
                                    y_fungi = gdm_alg_q1_splineDat$y[,"fun_q1_cqn"])
gdm_bac_q1_plot_fungi <- data.frame(x_fungi = gdm_bac_q1_splineDat$x[,"fun_q1_cqn"],
                                    y_fungi = gdm_bac_q1_splineDat$y[,"fun_q1_cqn"])

gdm_alg_q2_plot_fungi <- data.frame(x_fungi = gdm_alg_q2_splineDat$x[,"fun_q2_cqn"], 
                                    y_fungi = gdm_alg_q2_splineDat$y[,"fun_q2_cqn"])
gdm_bac_q2_plot_fungi <- data.frame(x_fungi = gdm_bac_q2_splineDat$x[,"fun_q2_cqn"],
                                    y_fungi = gdm_bac_q2_splineDat$y[,"fun_q2_cqn"])


######
# Effect on algal communities
######
algae_fungi_effect <- ggplot() +
  geom_line(data = gdm_alg_q0_plot_fungi, aes(x = x_fungi, y = y_fungi),
            color = algae_col, linewidth = 0.7) + 
  geom_text(data = gdm_alg_q0_plot_fungi, aes(x = max(x_fungi) + 0.01,
                                            y = max(y_fungi)), label = "*",
            size = 4)  +
  geom_line(data = gdm_alg_q1_plot_fungi, aes(x = x_fungi, y = y_fungi),
            color = algae_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_text(data = gdm_alg_q1_plot_fungi, aes(x = max(x_fungi) + 0.01,
                                            y = max(y_fungi)), label = "*",
            size = 4)  + 
  geom_line(data = gdm_alg_q2_plot_fungi, aes(x = x_fungi, y = y_fungi),
            color = algae_col, linewidth = 0.7, linetype = "dotted") + 
  geom_text(data = gdm_alg_q2_plot_fungi, aes(x = max(x_fungi) + 0.01,
                                            y = max(y_fungi)), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(x= "Fungal \u03B2-diversity",    
       y= "Algal \u03B2-diversity")+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2))
algae_fungi_effect

######
# Effect on bacterial communities
######
bacteria_fungi_effect <- ggplot() +
  geom_line(data = gdm_bac_q0_plot_fungi, aes(x = x_fungi, y = y_fungi),
            color = bacteria_col, linewidth = 0.7) + 
  geom_text(data = gdm_bac_q0_plot_fungi, aes(x = max(x_fungi) + 0.02,
                                              y = max(y_fungi)), label = "*",
            size = 4) + 
  geom_line(data = gdm_bac_q1_plot_fungi, aes(x = x_fungi, y = y_fungi),
            color = bacteria_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_text(data = gdm_bac_q1_plot_fungi, aes(x = max(x_fungi) + 0.02,
                                              y = max(y_fungi)), label = "*",
            size = 4) + 
  geom_line(data = gdm_bac_q2_plot_fungi, aes(x = x_fungi, y = y_fungi),
            color = bacteria_col, linewidth = 0.7, linetype = "dotted") + 
  geom_text(data = gdm_bac_q2_plot_fungi, aes(x = max(x_fungi) + 0.02,
                                              y = max(y_fungi)), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Fungal \u03B2-diversity",    
       y= "Bacterial \u03B2-diversity")+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2))
bacteria_fungi_effect

# algae beta effect
gdm_fun_q0_plot_algae <- data.frame(x_algae = gdm_fun_q0_splineDat$x[,"alg_q0_cqn"], 
                                    y_algae = gdm_fun_q0_splineDat$y[,"alg_q0_cqn"])
gdm_bac_q0_plot_algae <- data.frame(x_algae = gdm_bac_q0_splineDat$x[,"alg_q0_cqn"],
                                    y_algae = gdm_bac_q0_splineDat$y[,"alg_q0_cqn"])

gdm_fun_q1_plot_algae <- data.frame(x_algae = gdm_fun_q1_splineDat$x[,"alg_q1_cqn"], 
                                    y_algae = gdm_fun_q1_splineDat$y[,"alg_q1_cqn"])
gdm_bac_q1_plot_algae <- data.frame(x_algae = gdm_bac_q1_splineDat$x[,"alg_q1_cqn"],
                                    y_algae = gdm_bac_q1_splineDat$y[,"alg_q1_cqn"])

gdm_fun_q2_plot_algae <- data.frame(x_algae = gdm_fun_q2_splineDat$x[,"alg_q2_cqn"], 
                                    y_algae = gdm_fun_q2_splineDat$y[,"alg_q2_cqn"])
gdm_bac_q2_plot_algae <- data.frame(x_algae = gdm_bac_q2_splineDat$x[,"alg_q2_cqn"],
                                    y_algae = gdm_bac_q2_splineDat$y[,"alg_q2_cqn"])

######
# Effect on fungal communities
######
fungi_algae_effect <- ggplot() +
  geom_line(data = gdm_fun_q0_plot_algae, aes(x = x_algae, y = y_algae),
            color = fungi_col, linewidth = 0.7) + 
  geom_text(data = gdm_fun_q0_plot_algae, aes(x = max(x_algae) + 0.02,
                                              y = max(y_algae)), label = "*",
            size = 4) + 
  geom_line(data = gdm_fun_q1_plot_algae, aes(x = x_algae, y = y_algae),
            color = fungi_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_text(data = gdm_fun_q1_plot_algae, aes(x = max(x_algae) + 0.02,
                                              y = max(y_algae)), label = "*",
            size = 4) + 
  geom_line(data = gdm_fun_q2_plot_algae, aes(x = x_algae, y = y_algae),
            color = fungi_col, linewidth = 0.7, linetype = "dotted") + 
  geom_text(data = gdm_fun_q2_plot_algae, aes(x = max(x_algae) + 0.02,
                                              y = max(y_algae)), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(x= "Algal \u03B2-diversity",    
       y= "Fungal \u03B2-diversity")+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2))
fungi_algae_effect

######
# Effect on bacterial communities
######
bacteria_algae_effect <- ggplot() +
  geom_line(data = gdm_bac_q0_plot_algae, aes(x = x_algae, y = y_algae),
            color = bacteria_col, linewidth = 0.7) + 
  geom_text(data = gdm_bac_q0_plot_algae, aes(x = max(x_algae) + 0.02,
                                              y = max(y_algae)), label = "*",
            size = 4) + 
  geom_line(data = gdm_bac_q1_plot_algae, aes(x = x_algae, y = y_algae),
            color = bacteria_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_text(data = gdm_bac_q1_plot_algae, aes(x = max(x_algae) + 0.02,
                                              y = max(y_algae)), label = "*",
            size = 4) + 
  geom_line(data = gdm_bac_q2_plot_algae, aes(x = x_algae, y = y_algae),
            color = bacteria_col, linewidth = 0.7, linetype = "dotted") + 
  geom_text(data = gdm_bac_q2_plot_algae, aes(x = max(x_algae) + 0.02,
                                              y = max(y_algae)), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(x= "Algal \u03B2-diversity")+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2))
bacteria_algae_effect

# bacteria beta effect
gdm_fun_q0_plot_bacteria <- data.frame(x_bacteria = gdm_fun_q0_splineDat$x[,"bac_q0_cqn"], 
                                       y_bacteria = gdm_fun_q0_splineDat$y[,"bac_q0_cqn"])
gdm_alg_q0_plot_bacteria <- data.frame(x_bacteria = gdm_alg_q0_splineDat$x[,"bac_q0_cqn"],
                                       y_bacteria = gdm_alg_q0_splineDat$y[,"bac_q0_cqn"])

gdm_fun_q1_plot_bacteria <- data.frame(x_bacteria = gdm_fun_q1_splineDat$x[,"bac_q1_cqn"], 
                                       y_bacteria = gdm_fun_q1_splineDat$y[,"bac_q1_cqn"])
gdm_alg_q1_plot_bacteria <- data.frame(x_bacteria = gdm_alg_q1_splineDat$x[,"bac_q1_cqn"],
                                       y_bacteria = gdm_alg_q1_splineDat$y[,"bac_q1_cqn"])

gdm_fun_q2_plot_bacteria <- data.frame(x_bacteria = gdm_fun_q2_splineDat$x[,"bac_q2_cqn"], 
                                       y_bacteria = gdm_fun_q2_splineDat$y[,"bac_q2_cqn"])
gdm_alg_q2_plot_bacteria <- data.frame(x_bacteria = gdm_alg_q2_splineDat$x[,"bac_q2_cqn"],
                                       y_bacteria = gdm_alg_q2_splineDat$y[,"bac_q2_cqn"])

######
# Effect on algal communities
######
algae_bacteria_effect <- ggplot() +
  geom_line(data = gdm_alg_q0_plot_bacteria, aes(x = x_bacteria, y = y_bacteria),
            color = algae_col, linewidth = 0.7)  + 
  geom_text(data = gdm_alg_q0_plot_bacteria, aes(x = max(x_bacteria) + 0.02,
                                              y = max(y_bacteria)), label = "*",
            size = 4) +
  geom_line(data = gdm_alg_q1_plot_bacteria, aes(x = x_bacteria, y = y_bacteria),
            color = algae_col, linetype = "dotdash", linewidth = 0.7) + 
  geom_text(data = gdm_alg_q1_plot_bacteria, aes(x = max(x_bacteria) + 0.02,
                                              y = max(y_bacteria)), label = "*",
            size = 4) + 
  geom_line(data = gdm_alg_q2_plot_bacteria, aes(x = x_bacteria, y = y_bacteria),
            color = algae_col, linetype = "dotted", linewidth = 0.7) + 
  geom_text(data = gdm_alg_q2_plot_bacteria, aes(x = max(x_bacteria) + 0.02,
                                              y = max(y_bacteria)), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"),legend.key.width=unit(1.5,"cm")) +
  labs(x= "Bacterial \u03B2-diversity",    
       y= "Algal \u03B2-diversity") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 2))
algae_bacteria_effect

######
# Effect on fungal communities
######
fungi_bacteria_effect <- ggplot() +
  geom_line(data = gdm_fun_q0_plot_bacteria, aes(x = x_bacteria, y = y_bacteria),
            color = fungi_col, linewidth = 0.7) + 
  geom_text(data = gdm_fun_q0_plot_bacteria, aes(x = max(x_bacteria) + 0.02,
                                                 y = max(y_bacteria)), label = "*",
            size = 4) + 
  geom_line(data = gdm_fun_q1_plot_bacteria, aes(x = x_bacteria, y = y_bacteria),
            color = fungi_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_text(data = gdm_fun_q1_plot_bacteria, aes(x = max(x_bacteria) + 0.02,
                                                 y = max(y_bacteria)), label = "*",
            size = 4) + 
  geom_line(data = gdm_fun_q2_plot_bacteria, aes(x = x_bacteria, y = y_bacteria),
            color = fungi_col, linewidth = 0.7, linetype = "dotted") + 
  geom_text(data = gdm_fun_q2_plot_bacteria, aes(x = max(x_bacteria) + 0.02,
                                                 y = max(y_bacteria)), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Bacterial \u03B2-diversity",    
       y= "Fungal \u03B2-diversity") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 2))
fungi_bacteria_effect

# Stand density effect
gdm_fun_q0_plot_density <- data.frame(x_density = gdm_fun_q0_splineDat$x[,"stand_density_basal_area"], 
                                    y_density = gdm_fun_q0_splineDat$y[,"stand_density_basal_area"])
gdm_bac_q0_plot_density <- data.frame(x_density = gdm_bac_q0_splineDat$x[,"stand_density_basal_area"],
                                    y_density = gdm_bac_q0_splineDat$y[,"stand_density_basal_area"])
gdm_alg_q0_plot_density <- data.frame(x_density = gdm_alg_q0_splineDat$x[,"stand_density_basal_area"],
                                    y_density = gdm_alg_q0_splineDat$y[,"stand_density_basal_area"])

gdm_fun_q1_plot_density <- data.frame(x_density = gdm_fun_q1_splineDat$x[,"stand_density_basal_area"], 
                                    y_density = gdm_fun_q1_splineDat$y[,"stand_density_basal_area"])
gdm_bac_q1_plot_density <- data.frame(x_density = gdm_bac_q1_splineDat$x[,"stand_density_basal_area"],
                                    y_density = gdm_bac_q1_splineDat$y[,"stand_density_basal_area"])
gdm_alg_q1_plot_density <- data.frame(x_density = gdm_alg_q1_splineDat$x[,"stand_density_basal_area"],
                                    y_density = gdm_alg_q1_splineDat$y[,"stand_density_basal_area"])

gdm_fun_q2_plot_density <- data.frame(x_density = gdm_fun_q2_splineDat$x[,"stand_density_basal_area"], 
                                    y_density = gdm_fun_q2_splineDat$y[,"stand_density_basal_area"])
gdm_bac_q2_plot_density <- data.frame(x_density = gdm_bac_q2_splineDat$x[,"stand_density_basal_area"],
                                    y_density = gdm_bac_q2_splineDat$y[,"stand_density_basal_area"])
gdm_alg_q2_plot_density <- data.frame(x_density = gdm_alg_q2_splineDat$x[,"stand_density_basal_area"],
                                    y_density = gdm_alg_q2_splineDat$y[,"stand_density_basal_area"])

######
# Effect on algal communities
######
algae_density_effect <- ggplot() +
  geom_line(data = gdm_alg_q0_plot_density, aes(x = x_density, y = y_density),
            color = algae_col, linewidth = 0.7)  +
  geom_line(data = gdm_alg_q1_plot_density, aes(x = x_density, y = y_density),
            color = algae_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_line(data = gdm_alg_q2_plot_density, aes(x = x_density, y = y_density),
            color = algae_col, linewidth = 0.7, linetype = "dotted") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Stand density",    
       y= "Algal \u03B2-diversity")
algae_density_effect

######
# Effect on fungal communities
######
fungi_density_effect <- ggplot() +
  geom_line(data = gdm_fun_q0_plot_density, aes(x = x_density, y = y_density),
            color = fungi_col, linewidth = 0.7) + 
  geom_text(data = gdm_fun_q0_plot_density, aes(x = max(x_density) + 0.02,
                                                 y = max(y_density)), label = "*",
            size = 4) + 
  geom_line(data = gdm_fun_q1_plot_density, aes(x = x_density, y = y_density),
            color = fungi_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_line(data = gdm_fun_q2_plot_density, aes(x = x_density, y = y_density),
            color = fungi_col, linewidth = 0.7, linetype = "dotted") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Stand density",    
       y= "Fungal \u03B2-diversity")
fungi_density_effect

######
# Effect on bacterial communities
######
bacteria_density_effect <- ggplot() +
  geom_line(data = gdm_bac_q0_plot_density, aes(x = x_density, y = y_density),
            color = bacteria_col, linewidth = 0.7) + 
  geom_line(data = gdm_bac_q1_plot_density, aes(x = x_density, y = y_density),
            color = bacteria_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_line(data = gdm_bac_q2_plot_density, aes(x = x_density, y = y_density),
            color = bacteria_col, linewidth = 0.7, linetype = "dotted") + 
  geom_text(data = gdm_bac_q2_plot_density, aes(x = max(x_density) + 0.02,
                                                y = max(y_density)), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Stand density",    
       y= "Bacterial \u03B2-diversity")
bacteria_density_effect

# Surrounding Forest effect
gdm_fun_q0_plot_forest <- data.frame(x_forest = gdm_fun_q0_splineDat$x[,"RA_forest"], 
                                    y_forest = gdm_fun_q0_splineDat$y[,"RA_forest"])
gdm_bac_q0_plot_forest <- data.frame(x_forest = gdm_bac_q0_splineDat$x[,"RA_forest"],
                                    y_forest = gdm_bac_q0_splineDat$y[,"RA_forest"])
gdm_alg_q0_plot_forest <- data.frame(x_forest = gdm_alg_q0_splineDat$x[,"RA_forest"],
                                    y_forest = gdm_alg_q0_splineDat$y[,"RA_forest"])

gdm_fun_q1_plot_forest <- data.frame(x_forest = gdm_fun_q1_splineDat$x[,"RA_forest"], 
                                    y_forest = gdm_fun_q1_splineDat$y[,"RA_forest"])
gdm_bac_q1_plot_forest <- data.frame(x_forest = gdm_bac_q1_splineDat$x[,"RA_forest"],
                                    y_forest = gdm_bac_q1_splineDat$y[,"RA_forest"])
gdm_alg_q1_plot_forest <- data.frame(x_forest = gdm_alg_q1_splineDat$x[,"RA_forest"],
                                    y_forest = gdm_alg_q1_splineDat$y[,"RA_forest"])

gdm_fun_q2_plot_forest <- data.frame(x_forest = gdm_fun_q2_splineDat$x[,"RA_forest"], 
                                    y_forest = gdm_fun_q2_splineDat$y[,"RA_forest"])
gdm_bac_q2_plot_forest <- data.frame(x_forest = gdm_bac_q2_splineDat$x[,"RA_forest"],
                                    y_forest = gdm_bac_q2_splineDat$y[,"RA_forest"])
gdm_alg_q2_plot_forest <- data.frame(x_forest = gdm_alg_q2_splineDat$x[,"RA_forest"],
                                    y_forest = gdm_alg_q2_splineDat$y[,"RA_forest"])

######
# Effect on algal communities
######
algae_forest_effect <- ggplot() +
  geom_line(data = gdm_alg_q0_plot_forest, aes(x = x_forest, y = y_forest),
            color = algae_col, linewidth = 0.7)  +
  geom_line(data = gdm_alg_q1_plot_forest, aes(x = x_forest, y = y_forest),
            color = algae_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_line(data = gdm_alg_q2_plot_forest, aes(x = x_forest, y = y_forest),
            color = algae_col, linewidth = 0.7, linetype = "dotted") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Surrounding Forest",    
       y= "Algal \u03B2-diversity")
algae_forest_effect

######
# Effect on fungal communities
######
fungi_forest_effect <- ggplot() +
  geom_line(data = gdm_fun_q0_plot_forest, aes(x = x_forest, y = y_forest),
            color = fungi_col, linewidth = 0.7) + 
  geom_line(data = gdm_fun_q1_plot_forest, aes(x = x_forest, y = y_forest),
            color = fungi_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_line(data = gdm_fun_q2_plot_forest, aes(x = x_forest, y = y_forest),
            color = fungi_col, linewidth = 0.7, linetype = "dotted") + 
  geom_text(data = gdm_fun_q2_plot_forest, aes(x = max(x_forest) + 0.02,
                                                y = max(y_forest)), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Surrounding Forest",    
       y= "Fungal \u03B2-diversity")
fungi_forest_effect

######
# Effect on bacterial communities
######
bacteria_forest_effect <- ggplot() +
  geom_line(data = gdm_bac_q0_plot_forest, aes(x = x_forest, y = y_forest),
            color = bacteria_col, linewidth = 0.7) + 
  geom_line(data = gdm_bac_q1_plot_forest, aes(x = x_forest, y = y_forest),
            color = bacteria_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_text(data = gdm_bac_q1_plot_forest, aes(x = max(x_forest) + 0.02,
                                               y = max(y_forest)), label = "*",
            size = 4) + 
  geom_line(data = gdm_bac_q2_plot_forest, aes(x = x_forest, y = y_forest),
            color = bacteria_col, linewidth = 0.7, linetype = "dotted") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Surrounding Forest",    
       y= "Bacterial \u03B2-diversity")
bacteria_forest_effect

# Dominant Tree effect
gdm_fun_q0_plot_dominant <- data.frame(x_dominant = gdm_fun_q0_splineDat$x[,"dom_tot_ratio"], 
                                     y_dominant = gdm_fun_q0_splineDat$y[,"dom_tot_ratio"])
gdm_bac_q0_plot_dominant <- data.frame(x_dominant = gdm_bac_q0_splineDat$x[,"dom_tot_ratio"],
                                     y_dominant = gdm_bac_q0_splineDat$y[,"dom_tot_ratio"])
gdm_alg_q0_plot_dominant <- data.frame(x_dominant = gdm_alg_q0_splineDat$x[,"dom_tot_ratio"],
                                     y_dominant = gdm_alg_q0_splineDat$y[,"dom_tot_ratio"])

gdm_fun_q1_plot_dominant <- data.frame(x_dominant = gdm_fun_q1_splineDat$x[,"dom_tot_ratio"], 
                                     y_dominant = gdm_fun_q1_splineDat$y[,"dom_tot_ratio"])
gdm_bac_q1_plot_dominant <- data.frame(x_dominant = gdm_bac_q1_splineDat$x[,"dom_tot_ratio"],
                                     y_dominant = gdm_bac_q1_splineDat$y[,"dom_tot_ratio"])
gdm_alg_q1_plot_dominant <- data.frame(x_dominant = gdm_alg_q1_splineDat$x[,"dom_tot_ratio"],
                                     y_dominant = gdm_alg_q1_splineDat$y[,"dom_tot_ratio"])

gdm_fun_q2_plot_dominant <- data.frame(x_dominant = gdm_fun_q2_splineDat$x[,"dom_tot_ratio"], 
                                     y_dominant = gdm_fun_q2_splineDat$y[,"dom_tot_ratio"])
gdm_bac_q2_plot_dominant <- data.frame(x_dominant = gdm_bac_q2_splineDat$x[,"dom_tot_ratio"],
                                     y_dominant = gdm_bac_q2_splineDat$y[,"dom_tot_ratio"])
gdm_alg_q2_plot_dominant <- data.frame(x_dominant = gdm_alg_q2_splineDat$x[,"dom_tot_ratio"],
                                     y_dominant = gdm_alg_q2_splineDat$y[,"dom_tot_ratio"])

######
# Effect on algal communities
######
algae_dominant_effect <- ggplot() +
  geom_line(data = gdm_alg_q0_plot_dominant, aes(x = x_dominant, y = y_dominant),
            color = algae_col, linewidth = 0.7)  +
  geom_line(data = gdm_alg_q1_plot_dominant, aes(x = x_dominant, y = y_dominant),
            color = algae_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_line(data = gdm_alg_q2_plot_dominant, aes(x = x_dominant, y = y_dominant),
            color = algae_col, linewidth = 0.7, linetype = "dotted") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Dominant Trees:Total Trees",    
       y= "Algal \u03B2-diversity")
algae_dominant_effect

######
# Effect on fungal communities
######
fungi_dominant_effect <- ggplot() +
  geom_line(data = gdm_fun_q0_plot_dominant, aes(x = x_dominant, y = y_dominant),
            color = fungi_col, linewidth = 0.7) + 
  geom_line(data = gdm_fun_q1_plot_dominant, aes(x = x_dominant, y = y_dominant),
            color = fungi_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_line(data = gdm_fun_q2_plot_dominant, aes(x = x_dominant, y = y_dominant),
            color = fungi_col, linewidth = 0.7, linetype = "dotted") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Dominant Trees:Total Trees",    
       y= "Fungal \u03B2-diversity")
fungi_dominant_effect

######
# Effect on bacterial communities
######
bacteria_dominant_effect <- ggplot() +
  geom_line(data = gdm_bac_q0_plot_dominant, aes(x = x_dominant, y = y_dominant),
            color = bacteria_col, linewidth = 0.7) + 
  geom_line(data = gdm_bac_q1_plot_dominant, aes(x = x_dominant, y = y_dominant),
            color = bacteria_col, linewidth = 0.7, linetype = "dotdash") + 
  geom_line(data = gdm_bac_q2_plot_dominant, aes(x = x_dominant, y = y_dominant),
            color = bacteria_col, linewidth = 0.7, linetype = "dotted") + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(x= "Dominant Trees:Total Trees",    
       y= "Bacterial \u03B2-diversity")
bacteria_dominant_effect

#####
# Plotting
#####
blank <- ggplot() + theme_bw() + theme(panel.border = element_blank())

beta_legend <- ggpubr::get_legend(ggplot() +
  geom_line(data = gdm_alg_q0_plot_bacteria, aes(x = x_bacteria, y = y_bacteria, linetype = "solid"),
            color = "black", linewidth = 0.7)  +
  geom_line(data = gdm_alg_q1_plot_bacteria, aes(x = x_bacteria, y = y_bacteria, linetype = "dotdash"),
            color = "black", linewidth = 0.7) + 
  geom_line(data = gdm_alg_q2_plot_bacteria, aes(x = x_bacteria, y = y_bacteria, linetype = "dotted"),
            color = "black", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), legend.position = "bottom") +
  labs(x= "Bacterial \u03B2-diversity",    
       y= "Algal \u03B2-diversity") +
  scale_linetype_manual(name = element_text("Diversity Level"), 
                        values =c('solid',
                                  'dotdash', 
                                  'dotted'),
                        labels = c('q = 0 (all species)', 'q = 1 (typical)', 'q = 2 (dominant)'))+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2)))


beta_effects_algae1 <- ggpubr::ggarrange(algae_bacteria_effect, 
                                        algae_fungi_effect,
                                        nrow = 1, ncol = 2)  %>% 
  ggpubr::annotate_figure(top =  ggpubr::text_grob("Biotic Effects", size = 7, family = "sans", face = "bold"))

beta_effects_algae2 <- ggpubr::ggarrange(algae_canopy_effect,
                                         algae_temperature_effect,
                                         nrow = 1, ncol = 2) %>% 
  ggpubr::annotate_figure(top = ggpubr::text_grob("Abiotic Effects", size = 7, family = "sans", face = "bold"))

beta_effects_algae3 <- ggpubr::ggarrange(beta_effects_algae1,
                                         beta_effects_algae2,
                                      nrow = 1, ncol = 2) 

beta_effects_fungi <- ggpubr::ggarrange(fungi_bacteria_effect,
                                        fungi_algae_effect, 
                                        fungi_canopy_effect, 
                                        fungi_temperature_effect,
                                        nrow = 1, ncol = 4)  %>% 
  ggpubr::annotate_figure(top = ggpubr::text_grob("",  color = fungi_col, face = "bold", size = 7))
                                    
beta_effects_bacteria <- ggpubr::ggarrange(bacteria_fungi_effect, 
                                           bacteria_algae_effect, 
                                           bacteria_canopy_effect,
                                           bacteria_temperature_effect,
                                           nrow = 1, ncol = 4)  %>% 
  ggpubr::annotate_figure(top = ggpubr::text_grob("",  color = bacteria_col, face = "bold", size = 7))

effects_beta <- ggpubr::ggarrange(beta_effects_algae3, beta_effects_fungi, beta_effects_bacteria,
                                  nrow = 3, ncol = 1, legend.grob = beta_legend, legend = "bottom",
                                  heights = c(1, 0.925, 0.925)) 


effects_beta
ggplot2::ggsave("effects_beta.pdf", plot = effects_beta, device = cairo_pdf,
       width = 175, height = 150, units = "mm")

ggplot2::ggsave("effects_beta.png", plot = effects_beta, device = png,
       width = 175, height = 150, units = "mm")


###
#Community Barplots
###
# Setting the colors for the community bar graphs. 
my_cols <- paletteer::paletteer_d('ggsci::default_igv')

#ALGAE

# Aggregate the taxa at order level. This is to allow better visualization, 
# while still retaining the most information.
phy_order_alg <- phy_algae_bark %>%
  microbiome::aggregate_taxa(level = "Order")

# Subset the phyloseq object to the top 24 orders and put the rest in 
# a category "Others", based on relative abundance. 
phy_alg_ord_top25 <- fantaxtic::top_taxa(phy_algae_bark,
                    tax_level = 'Order',
                    n_taxa =  24,
                    by_proportion = TRUE,
                    merged_label = "Other",
                    include_na_taxa = T) 

phy_alg_ord_top25_named <- fantaxtic::name_na_taxa(phy_alg_ord_top25$ps_obj, include_rank = T)

# Transform the subset dataset to compositional (relative) abundances.
phy_alg_ord_top25_plot <-  phy_alg_ord_top25_named %>%
  microbiome::aggregate_taxa(level = "Order") %>%   
  microbiome::transform(transform = "compositional") 

# Extract the names of the Orders.
taxa_names(phy_alg_ord_top25_plot) <- tax_table(phy_alg_ord_top25_plot)[, 4]

# Sort the taxa names alphabetically. 
taxa_names_alg_ord <- sort(taxa_names(phy_alg_ord_top25_plot))

# To get our desired plotting order and group names we need to change 
# the exploratory names and order them as factors.
sampledata_algae <- data.frame(sample_data(phy_alg_ord_top25_plot))
sampledata_algae <- sampledata_algae %>% 
  mutate(across("exploratory", stringr::str_replace, "Alb", "Swabian Alb")) %>% 
  mutate(across("exploratory", stringr::str_replace, "Hainich", "Hainich-Dün")) %>% 
  mutate(across("exploratory", stringr::str_replace, "Schorfheide", "Schorfheide-Chorin")) 

sampledata_algae$exploratory <- factor(sampledata_algae$exploratory, 
                                          levels = c("Swabian Alb", "Hainich-Dün", "Schorfheide-Chorin"))  

sample_data(phy_alg_ord_top25_plot) <- sample_data(sampledata_algae)

# Custom plotting to make a nice stacked barplot. 
alg_ord_plots <- phy_alg_ord_top25_plot %>%
  microbiome::plot_composition(group_by =  'exploratory', otu.sort = taxa_names_alg_ord) +
  scale_fill_manual(values = ggplot2::alpha(my_cols, 0.9), name = 'Order') +
  guides(fill = guide_legend(title.position = 'top')) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', linewidth = 0.25),
        axis.text.x =  element_blank(),
        axis.text.y =  element_text(colour = "black"),
        axis.title = element_text(colour = "black"),
        legend.position = 'bottom', 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 7),
        legend.title =  element_text(size = 7),
        legend.key.size = unit(2.5, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        legend.background = element_rect(fill = 'transparent'),
        text = element_text(colour = 'black', size = 7)) + 
  xlab('Sample') +
  ylab('Relative abundance')  +
  scale_y_continuous(label = scales::percent)  + 
  labs( subtitle = 'Algae')

alg_ord_plots

#BACTERIA

# Aggregate the taxa at order level. This is to allow better visualization, 
# while still retaining the most information.

# Subset the phyloseq object to the top 24 orders and put the rest in 
# a category "Others", based on relative abundance. 
phy_bac_ord_top25 <- fantaxtic::top_taxa(phy_bacteria_bark,
                                         tax_level = 'Order',
                                         n_taxa =  24,
                                         by_proportion = TRUE,
                                         merged_label = "Other",
                                         include_na_taxa = T) 

phy_bac_ord_top25_named <- fantaxtic::name_na_taxa(phy_bac_ord_top25$ps_obj, include_rank = T)

# Transform the subset dataset to compositional (relative) abundances.
phy_bac_ord_top25_plot <-  phy_bac_ord_top25_named %>%
  microbiome::aggregate_taxa(level = "Order") %>%  
  microbiome::transform(transform = "compositional") 

# Extract the names of the Orders.
taxa_names(phy_bac_ord_top25_plot) <- tax_table(phy_bac_ord_top25_plot)[, 4]

# Sort the taxa names alphabetically. 
taxa_names_bac_ord <- sort(taxa_names(phy_bac_ord_top25_plot))

# To get our desired plotting order and group names we need to change 
# the exploratory names and order them as factors.
sampledata_bacteria <- data.frame(sample_data(phy_bac_ord_top25_plot))
sampledata_bacteria <- sampledata_bacteria %>% 
  mutate(across("exploratory", stringr::str_replace, "Alb", "Swabian Alb")) %>% 
  mutate(across("exploratory", stringr::str_replace, "Hainich", "Hainich-Dün")) %>% 
  mutate(across("exploratory", stringr::str_replace, "Schorfheide", "Schorfheide-Chorin")) 

sampledata_bacteria$exploratory <- factor(sampledata_bacteria$exploratory, 
                                       levels = c("Swabian Alb", "Hainich-Dün", "Schorfheide-Chorin"))  

sample_data(phy_bac_ord_top25_plot) <- sample_data(sampledata_bacteria)


# Custom plotting to make a nice stacked barplot. 
bac_ord_plots <- phy_bac_ord_top25_plot %>%
  microbiome::plot_composition( group_by =  'exploratory', otu.sort = taxa_names_bac_ord) +
  scale_fill_manual(values = ggplot2::alpha(my_cols, 0.9), name = 'Order') +
  guides(fill = guide_legend(title.position = 'top')) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', linewidth = 0.25),
        axis.text.x =  element_blank(),
        axis.text.y =  element_text(colour = "black"),
        axis.title = element_text(colour = "black"),
        legend.position = 'bottom', 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 7),
        legend.title =  element_text(size = 7),
        legend.key.size = unit(2.5, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        legend.background = element_rect(fill = 'transparent'),
        text = element_text(colour = 'black', size = 7)) + 
  xlab('Sample') +
  ylab('Relative abundance')  +
  scale_y_continuous(label = scales::percent)  + 
  labs( subtitle = 'Bacteria') 

bac_ord_plots

#FUNGI

# Aggregate the taxa at order level. This is to allow better visualization, 
# while still retaining the most information.


# Subset the phyloseq object to the top 24 orders and put the rest in 
# a category "Others", based on relative abundance. 
phy_fun_ord_top25 <- fantaxtic::top_taxa(phy_fungi_bark,
                                         tax_level = 'Order',
                                         n_taxa =  24,
                                         by_proportion = TRUE,
                                         merged_label = "Other",
                                         include_na_taxa = T) 

phy_fun_ord_top25_named <- fantaxtic::name_na_taxa(phy_fun_ord_top25$ps_obj, include_rank = T)


# Transform the subset dataset to compositional (relative) abundances.
phy_fun_ord_top25_plot <-  phy_fun_ord_top25_named %>%
  microbiome::aggregate_taxa(level = "Order") %>%  
  microbiome::transform(transform = "compositional") 

# Extract the names of the Orders.
taxa_names(phy_fun_ord_top25_plot) <- tax_table(phy_fun_ord_top25_plot)[, 4]

# Sort the taxa names alphabetically. 
taxa_names_fun_ord <- sort(taxa_names(phy_fun_ord_top25_plot))

# To get our desired plotting order and group names we need to change 
# the exploratory names and order them as factors.
sampledata_fungi <- data.frame(sample_data(phy_fun_ord_top25_plot))
sampledata_fungi <- sampledata_fungi %>% 
  mutate(across("exploratory", stringr::str_replace, "Alb", "Swabian Alb")) %>% 
  mutate(across("exploratory", stringr::str_replace, "Hainich", "Hainich-Dün")) %>% 
  mutate(across("exploratory", stringr::str_replace, "Schorfheide", "Schorfheide-Chorin")) 

sampledata_fungi$exploratory <- factor(sampledata_fungi$exploratory, 
                                       levels = c("Swabian Alb", "Hainich-Dün", "Schorfheide-Chorin"))  

sample_data(phy_fun_ord_top25_plot) <- sample_data(sampledata_fungi)

# Custom plotting to make a nice stacked barplot. 
fun_ord_plots <- phy_fun_ord_top25_plot %>%
  microbiome::plot_composition(group_by =  'exploratory', otu.sort = taxa_names_fun_ord) +
  scale_fill_manual(values = ggplot2::alpha(my_cols, 0.9), name = 'Order') +
  guides(fill = guide_legend(title.position = 'top')) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = 'black', linewidth = 0.25),
        axis.text.x =  element_blank(),
        axis.text.y =  element_text(colour = "black"),
        axis.title = element_text(colour = "black"),
        legend.position = 'bottom', 
        plot.title = element_text(vjust = -4, hjust = 0.03), 
        legend.text = element_text(colour = 'black', size = 7),
        legend.title =  element_text(size = 7),
        legend.key.size = unit(2.5, 'mm'),
        axis.ticks.length.x = unit(-0.2, "cm"), 
        legend.box.spacing = unit(-4, 'mm'),
        legend.background = element_rect(fill = 'transparent'),
        text = element_text(colour = 'black', size = 7)) + 
  xlab('Sample') +
  ylab('Relative abundance')  +
  scale_y_continuous(label = scales::percent)  + 
  labs( subtitle = 'Fungi')
fun_ord_plots


# List the plots and combine them to make the final figure.

finished_barplots <- ggpubr::ggarrange(alg_ord_plots, fun_ord_plots, bac_ord_plots,
                                       nrow = 3, ncol = 1)
finished_barplots

ggplot2::ggsave("finished_barplots.pdf", plot = finished_barplots, device = cairo_pdf,
       width = 175, height = 250, units = "mm")

ggplot2::ggsave("finished_barplots.png", plot = finished_barplots, device = png,
                width = 175, height = 250, units = "mm")

###
#Variance partitioning plot
###

# Import the results of the Alpha diversity variation partioning
alg_variance_lm <- readRDS(here::here("alg_variance_lm.rds"))
fun_variance_lm <- readRDS(here::here("fun_variance_lm.rds"))
bac_variance_lm <- readRDS(here::here("bac_variance_lm.rds"))

# Combine alpha into one dataframe 

alg_variance_lm <- alg_variance_lm  %>% 
  cbind(div_lev = rep("\u03B1-diversity", nrow(alg_variance_lm))) %>% 
  cbind(organism = rep("Algae", nrow(alg_variance_lm)))

fun_variance_lm <- fun_variance_lm  %>% 
  cbind(div_lev = rep("\u03B1-diversity", nrow(fun_variance_lm))) %>% 
  cbind(organism = rep("Fungi", nrow(fun_variance_lm)))

bac_variance_lm <- bac_variance_lm  %>% 
  cbind(div_lev = rep("\u03B1-diversity", nrow(bac_variance_lm))) %>% 
  cbind(organism = rep("Bacteria", nrow(bac_variance_lm)))

variance_lm <- rbind(alg_variance_lm, fun_variance_lm, bac_variance_lm) %>% 
  dplyr::rename(variance = Proportion) %>% 
  dplyr::mutate(variance = (variance * 100))

# Combine beta into one dataframe 

alg_variance_gdm <- alg_variance  %>% 
  cbind(div_lev = rep("\u03B2-diversity", nrow(alg_variance))) %>% 
  cbind(organism = rep("Algae", nrow(alg_variance)))

fun_variance_gdm <- fun_variance  %>% 
  cbind(div_lev = rep("\u03B2-diversity", nrow(fun_variance))) %>% 
  cbind(organism = rep("Fungi", nrow(alg_variance)))

bac_variance_gdm <- bac_variance %>% 
  cbind(div_lev = rep("\u03B2-diversity", nrow(bac_variance))) %>% 
  cbind(organism = rep("Bacteria", nrow(alg_variance)))

variance_beta <- rbind(alg_variance_gdm, fun_variance_gdm, bac_variance_gdm) %>% 
  dplyr::rename(variable = VARIABLE_SET,
                variance = DEVIANCE)

variance_full <- rbind(variance_lm, variance_beta)  %>% 
  mutate(organism = factor(organism, levels = c("Algae", "Fungi", "Bacteria")))
  
# Create faceted multipanel plot with ggpubr and ggplot2

div_labels <- c("\u03B1-diversity", "\u03B2-diversity")

variance_barplot <- ggpubr::ggbarplot(variance_full, x = "q_lev", y = "variance", fill = "variable",
                                      color = NA, 
                                      palette = c("#999999", "#99876d", "#78AB7A", "#E5914A",
                                                        "#6D5882", "#F0E442", "#0072B2", "#DA3E52")) +
  ggplot2::facet_grid(organism ~ div_lev, space="free", scales="free", switch = "y") +
  ggplot2::scale_x_discrete(position = "top") +
  ggplot2::geom_hline(aes(yintercept = 0), linetype = "dotted", linewidth = 0.3) +
  ggplot2::coord_flip() + 
  ggplot2::ylab("Explained variance") +
  ggplot2::theme(text = element_text(size = 7),
                 legend.title = element_blank(),
                 legend.position = "right",
                 axis.title.y = element_blank(),
                 axis.line.y = element_blank(), 
                 axis.ticks.y = element_blank(),
                 axis.line.x = element_line(linewidth = 0.3),
                 strip.background = element_rect(color = "grey")) 
variance_barplot

ggplot2::ggsave("variance_barplot.pdf", plot = variance_barplot, device = cairo_pdf,
       width = 175, height = 87.5, units = "mm")

ggplot2::ggsave("variance_barplot.png", plot = variance_barplot, device = png,
                width = 175, height = 87.5, units = "mm")

# Miscallenous 

phybarkalb_alg <- subset_samples(phy_algae_bark, exploratory == "Alb")
phybarkalb_alg <- phyloseq::prune_taxa(phyloseq::taxa_sums(phybarkalb_alg) != 0, phybarkalb_alg)
phybarkalb_alg

phybarkhai_alg <- subset_samples(phy_algae_bark, exploratory =="Hainich")
phybarkhai_alg <- phyloseq::prune_taxa(phyloseq::taxa_sums(phybarkhai_alg) != 0, phybarkhai_alg)
phybarkhai_alg

phybarksch_alg <- subset_samples(phy_algae_bark, exploratory == "Schorfheide")
phybarksch_alg <- phyloseq::prune_taxa(phyloseq::taxa_sums(phybarksch_alg) != 0, phybarksch_alg)
phybarksch_alg

phybarkalb_fun <- subset_samples(phy_fungi_bark, exploratory == "Alb")
phybarkalb_fun <- phyloseq::prune_taxa(phyloseq::taxa_sums(phybarkalb_fun) != 0, phybarkalb_fun)
phybarkalb_fun

phybarkhai_fun <- subset_samples(phy_fungi_bark, exploratory == "Hainich")
phybarkhai_fun <- phyloseq::prune_taxa(phyloseq::taxa_sums(phybarkhai_fun) != 0, phybarkhai_fun)
phybarkhai_fun

phybarksch_fun <- subset_samples(phy_fungi_bark, exploratory == "Schorfheide")
phybarksch_fun <- phyloseq::prune_taxa(phyloseq::taxa_sums(phybarksch_fun) != 0, phybarksch_fun)
phybarksch_fun

phybarkalb_bac <- subset_samples(phy_bacteria_bark, exploratory == "Alb")
phybarkalb_bac <- phyloseq::prune_taxa(phyloseq::taxa_sums(phybarkalb_bac) != 0, phybarkalb_bac)
phybarkalb_bac

phybarkhai_bac <- subset_samples(phy_bacteria_bark, exploratory == "Hainich")
phybarkhai_bac <- phyloseq::prune_taxa(phyloseq::taxa_sums(phybarkhai_bac) != 0, phybarkhai_bac)
phybarkhai_bac

phybarksch_bac <- subset_samples(phy_bacteria_bark, exploratory == "Schorfheide")
phybarksch_bac <- phyloseq::prune_taxa(phyloseq::taxa_sums(phybarksch_bac) != 0, phybarksch_bac)
phybarksch_bac

