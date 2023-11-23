#################################################################
##                          Section 1                          ##
##                       Package Loading                       ##
#################################################################

library(here); packageVersion("here")
#1.0.1
library(effects); packageVersion("effects")
#4.2.2
library(hillR); packageVersion("hillR")
#0.5.1
library(tidyverse); packageVersion("tidyverse")
#2.0.0
library(phyloseq); packageVersion("phyloseq")
#1.44.0
library(ggeffects); packageVersion("ggeffects")
#1.2.3

#################################################################
##                          Section 2                          ##
##                        Data Loading                         ##
#################################################################
#Read in the phyloseq files from the data cleaning section. 
phy_algae_bark <- base::readRDS(here("Data", "phy_algae_bark.rds"))

phy_bacteria_bark <- base::readRDS(here("Data", "phy_bacteria_bark_new.rds"))

phy_fungi_bark <- base::readRDS(here("Data", "phy_fungi_bark.rds"))

metadata_bark <- base::data.frame(phyloseq::sample_data(phy_fungi_bark)) %>% 
  tibble::rownames_to_column(var = "Sample_ID") %>%  
  dplyr::select(-one_of("tree_type", "stand_density_abundance",
                        "precipitation_radolan", "enl_2019",
                        "stand_evenness_basal_area", "d_SD",
                        "PAR", "Plot_ID", "substrate"))


#################################################################
##                          Section 3                          ##
##               Linear Models of Alpha Diversity              ##
#################################################################

####
# First Step 
####

# Simple linear model of y ~ Exploratory + host tree + tree independent variables
# and the other organismal groups. The response variables are Hill Numbers 
# q = 0, q = 1, q = 2.  

# Calculate the hill numbers.
# First we need to extract the ASV table from the phyloseq object.

asv_algae_hill_bark <- base::data.frame(phyloseq::otu_table(phy_algae_bark)) 

alg_q0_bark_alpha <- base::data.frame("alg_q0" = hillR::hill_taxa(asv_algae_hill_bark, q = 0, MARGIN = 2)) %>% 
  tibble::rownames_to_column(var = "Sample_ID")
alg_q1_bark_alpha <- base::data.frame("alg_q1" = hillR::hill_taxa(asv_algae_hill_bark, q = 1, MARGIN = 2)) %>% 
  tibble::rownames_to_column(var = "Sample_ID")
alg_q2_bark_alpha <- base::data.frame("alg_q2" = hillR::hill_taxa(asv_algae_hill_bark, q = 2, MARGIN = 2)) %>% 
  tibble::rownames_to_column(var = "Sample_ID")

asv_bacteria_hill_bark <- base::data.frame(phyloseq::otu_table(phy_bacteria_bark)) 

bac_q0_bark_alpha <- base::data.frame("bac_q0" = hillR::hill_taxa(asv_bacteria_hill_bark, q = 0, MARGIN = 2)) %>% 
  tibble::rownames_to_column(var = "Sample_ID")
bac_q1_bark_alpha <- base::data.frame("bac_q1" = hillR::hill_taxa(asv_bacteria_hill_bark, q = 1, MARGIN = 2)) %>% 
  tibble::rownames_to_column(var = "Sample_ID")
bac_q2_bark_alpha <- base::data.frame("bac_q2" = hillR::hill_taxa(asv_bacteria_hill_bark, q = 2, MARGIN = 2)) %>% 
  tibble::rownames_to_column(var = "Sample_ID")

asv_fungi_hill_bark <- base::data.frame(phyloseq::otu_table(phy_fungi_bark)) 

fun_q0_bark_alpha <- base::data.frame("fun_q0" = hillR::hill_taxa(asv_fungi_hill_bark, q = 0, MARGIN = 2)) %>% 
  tibble::rownames_to_column(var = "Sample_ID")
fun_q1_bark_alpha <- base::data.frame("fun_q1" = hillR::hill_taxa(asv_fungi_hill_bark, q = 1, MARGIN = 2)) %>% 
  tibble::rownames_to_column(var = "Sample_ID")
fun_q2_bark_alpha <- base::data.frame("fun_q2" = hillR::hill_taxa(asv_fungi_hill_bark, q = 2, MARGIN = 2)) %>% 
  tibble::rownames_to_column(var = "Sample_ID")

# Join with the rest of the metadata.
metadata_alpha <- metadata_bark  %>% 
  dplyr::inner_join(.,alg_q0_bark_alpha) %>% 
  dplyr::inner_join(.,alg_q1_bark_alpha) %>% 
  dplyr::inner_join(.,alg_q2_bark_alpha) %>% 
  dplyr::inner_join(.,bac_q0_bark_alpha) %>% 
  dplyr::inner_join(.,bac_q1_bark_alpha) %>% 
  dplyr::inner_join(.,bac_q2_bark_alpha) %>% 
  dplyr::inner_join(.,fun_q0_bark_alpha) %>% 
  dplyr::inner_join(.,fun_q1_bark_alpha) %>% 
  dplyr::inner_join(.,fun_q2_bark_alpha) 

# Calculate the library size and append it to the metadata. 
alg_library <- data.frame(library_size_alg = phyloseq::sample_sums(phy_algae_bark)) %>% 
  tibble::rownames_to_column(var = "Sample_ID")

fun_library <- data.frame(library_size_fun = phyloseq::sample_sums(phy_fungi_bark)) %>% 
  tibble::rownames_to_column(var = "Sample_ID")

bac_library <- data.frame(library_size_bac = phyloseq::sample_sums(phy_bacteria_bark)) %>% 
  tibble::rownames_to_column(var = "Sample_ID")

metadata_alpha <- metadata_alpha %>% 
  dplyr::inner_join(., alg_library) %>% 
  dplyr::inner_join(., fun_library) %>% 
  dplyr::inner_join(., bac_library)

# Scale the variables to make effects comparable between the models. 
metadata_alpha_scaled <- metadata_alpha %>% 
  dplyr::mutate(across(.cols = c("alg_q0", "alg_q1", "alg_q2",
                                 "bac_q0", "bac_q1", "bac_q2",
                                 "fun_q0", "fun_q1", "fun_q2"),~(scale(.) %>% as.numeric))) %>% 
  # Furthermore, scale the explanatory variables since they have very different units. 
  dplyr::mutate(across(.cols = c("rH_200", "Ta_200", "stand_density_basal_area", 
                                 "DBH_avg", "d_gini", "RA_forest", 
                                 "canopy_openness_2019", "dom_tot_ratio"),~(scale(.) %>% as.numeric))) %>% 
  dplyr::mutate(across(.cols = c("library_size_alg", "library_size_fun",
                                 "library_size_bac"), ~(log(.) %>% as.numeric)))

metadata_alpha_scaled$exploratory <- as.factor(metadata_alpha_scaled$exploratory)

############
# Linear Models of alpha diversity
###########
# Create linear models including all of the previously chosen variables. 
# We do not employ any model selection to make effect sizes etc comparable between organism groups
# and levels of alpha diversity. 
#ALGAE
lm_q0_alg_second <- stats::lm(alg_q0 ~ exploratory + 
                         rH_200 + 
                         Ta_200 +
                         stand_density_basal_area +
                         DBH_avg +
                         d_gini +
                         canopy_openness_2019 +
                         dom_tot_ratio + 
                         RA_forest + 
                         bac_q0 +
                         fun_q0 + 
                         offset(library_size_alg),
                       data = metadata_alpha_scaled,
                       na.action = "na.fail")
base::summary(lm_q0_alg_second)  
base::plot(allEffects(lm_q0_alg_second))
base::plot(stats::residuals(lm_q0_alg_second), stats::fitted(lm_q0_alg_second))
stats::qqnorm(stats::residuals(lm_q0_alg_second))
stats::qqline(stats::residuals(lm_q0_alg_second))

lm_q1_alg_second <- stats::lm(alg_q1 ~ exploratory + 
                         rH_200 + 
                         Ta_200 +
                         stand_density_basal_area +
                         DBH_avg +
                         d_gini +
                         canopy_openness_2019 +
                         dom_tot_ratio + 
                         RA_forest + 
                         bac_q1 +
                         fun_q1 + 
                           offset(library_size_alg),
                       data = metadata_alpha_scaled,
                       na.action = "na.fail")
base::summary(lm_q1_alg_second)  
base::plot(allEffects(lm_q1_alg_second))
base::plot(stats::residuals(lm_q1_alg_second), stats::fitted(lm_q1_alg_second))
stats::qqnorm(stats::residuals(lm_q1_alg_second))
stats::qqline(stats::residuals(lm_q1_alg_second))

lm_q2_alg_second <- stats::lm(alg_q2 ~ exploratory + 
                         rH_200 + 
                         Ta_200 +
                         stand_density_basal_area +
                         DBH_avg +
                         d_gini +
                         canopy_openness_2019 +
                         dom_tot_ratio + 
                         RA_forest + 
                         bac_q2 +
                         fun_q2 + 
                           offset(library_size_alg),
                       data = metadata_alpha_scaled,
                       na.action = "na.fail")
base::summary(lm_q2_alg_second)  
base::plot(allEffects(lm_q2_alg_second))
base::plot(stats::residuals(lm_q2_alg_second), stats::fitted(lm_q2_alg_second))
stats::qqnorm(stats::residuals(lm_q2_alg_second))
stats::qqline(stats::residuals(lm_q2_alg_second))

# Bacteria
lm_q0_bac_second <- stats::lm(bac_q0 ~ exploratory + 
                         rH_200 + 
                         Ta_200 +
                         stand_density_basal_area +
                         DBH_avg +
                         d_gini +
                         canopy_openness_2019 +
                         dom_tot_ratio + 
                         RA_forest + 
                         fun_q0 +
                         alg_q0 + 
                           offset(library_size_bac),
                       data = metadata_alpha_scaled,
                       na.action = "na.fail")
base::summary(lm_q0_bac_second)  
base::plot(allEffects(lm_q0_bac_second))
base::plot(stats::residuals(lm_q0_bac_second), stats::fitted(lm_q0_bac_second))
stats::qqnorm(stats::residuals(lm_q0_bac_second))
stats::qqline(stats::residuals(lm_q0_bac_second))

lm_q1_bac_second <- stats::lm(bac_q1 ~ exploratory + 
                         rH_200 + 
                         Ta_200 +
                         stand_density_basal_area +
                         DBH_avg +
                         d_gini +
                         canopy_openness_2019 +
                         dom_tot_ratio + 
                         RA_forest + 
                         fun_q1 +
                         alg_q1 + 
                           offset(library_size_bac),
                       data = metadata_alpha_scaled,
                       na.action = "na.fail")
base::summary(lm_q1_bac_second)  
base::plot(allEffects(lm_q1_bac_second))
base::plot(stats::residuals(lm_q1_bac_second), stats::fitted(lm_q1_bac_second))
stats::qqnorm(stats::residuals(lm_q1_bac_second))
stats::qqline(stats::residuals(lm_q1_bac_second))

lm_q2_bac_second <- stats::lm(bac_q2 ~ exploratory + 
                         rH_200 + 
                         Ta_200 +
                         stand_density_basal_area +
                         DBH_avg +
                         d_gini +
                         canopy_openness_2019 +
                         dom_tot_ratio + 
                         RA_forest + 
                         fun_q2 +
                         alg_q2 + 
                           offset(library_size_bac),
                       data = metadata_alpha_scaled,
                       na.action = "na.fail")
base::summary(lm_q2_bac_second)  
base::plot(allEffects(lm_q2_bac_second))
base::plot(stats::residuals(lm_q2_bac_second), stats::fitted(lm_q2_bac_second))
stats::qqnorm(stats::residuals(lm_q2_bac_second))
stats::qqline(stats::residuals(lm_q2_bac_second))

# Fungi
lm_q0_fun_second <- stats::lm(fun_q0 ~ exploratory + 
                         rH_200 + 
                         Ta_200 +
                         stand_density_basal_area +
                         DBH_avg +
                         d_gini +
                         canopy_openness_2019 +
                         dom_tot_ratio + 
                         RA_forest + 
                         bac_q0 +
                         alg_q0 + 
                           offset(library_size_fun),
                       data = metadata_alpha_scaled,
                       na.action = "na.fail")
base::summary(lm_q0_fun_second)  
base::plot(allEffects(lm_q0_fun_second))
base::plot(stats::residuals(lm_q0_fun_second), stats::fitted(lm_q0_fun_second))
stats::qqnorm(stats::residuals(lm_q0_fun_second))
stats::qqline(stats::residuals(lm_q0_fun_second))

lm_q1_fun_second <- stats::lm(fun_q1 ~ exploratory + 
                         rH_200 + 
                         Ta_200 +
                         stand_density_basal_area +
                         DBH_avg +
                         d_gini +
                         canopy_openness_2019 +
                         dom_tot_ratio + 
                         RA_forest + 
                         bac_q1 +
                         alg_q1 + 
                           offset(library_size_fun),
                       data = metadata_alpha_scaled,
                       na.action = "na.fail")
base::summary(lm_q1_fun_second)  
base::plot(allEffects(lm_q1_fun_second))
base::plot(stats::residuals(lm_q1_fun_second), stats::fitted(lm_q1_fun_second))
stats::qqnorm(stats::residuals(lm_q1_fun_second))
stats::qqline(stats::residuals(lm_q1_fun_second))

lm_q2_fun_second <- stats::lm(fun_q2 ~ exploratory + 
                         rH_200 + 
                         Ta_200 +
                         stand_density_basal_area +
                         DBH_avg +
                         d_gini +
                         canopy_openness_2019 +
                         dom_tot_ratio + 
                         RA_forest + 
                         bac_q2 +
                         alg_q2 + 
                           offset(library_size_fun),
                       data = metadata_alpha_scaled,
                       na.action = "na.fail")
base::summary(lm_q2_fun_second)  
base::plot(allEffects(lm_q2_fun_second))
base::plot(stats::residuals(lm_q2_fun_second), stats::fitted(lm_q2_fun_second))
stats::qqnorm(stats::residuals(lm_q2_fun_second))
stats::qqline(stats::residuals(lm_q2_fun_second))

# Calculate Benjamini Hochberg corrected p-values of all of the linear models. 
p_step2 <- (summary(lm_q0_alg_second)$coefficients[-1,4])
p_step2 <- c(p_step2, (summary(lm_q1_alg_second)$coefficients[-1,4]))
p_step2 <- c(p_step2, (summary(lm_q2_alg_second)$coefficients[-1,4]))

p_step2 <- c(p_step2, summary(lm_q0_bac_second)$coefficients[-1,4])
p_step2 <- c(p_step2, (summary(lm_q1_bac_second)$coefficients[-1,4]))
p_step2 <- c(p_step2, (summary(lm_q2_bac_second)$coefficients[-1,4]))

p_step2 <- c(p_step2, summary(lm_q0_fun_second)$coefficients[-1,4])
p_step2 <- c(p_step2, (summary(lm_q1_fun_second)$coefficients[-1,4]))
p_step2 <- c(p_step2, (summary(lm_q2_fun_second)$coefficients[-1,4]))

p_adj_step2 <- p.adjust(p_step2, method = "fdr")

p_vals_step2 <- data.frame(level = c(rep("algae_q0", 12),
                                     rep("algae_q1", 12),
                                     rep("algae_q2", 12),
                                     rep("bacteria_q0", 12),
                                     rep("bacteria_q1", 12),
                                     rep("bacteria_q2", 12),
                                     rep("fungi_q0", 12),
                                     rep("fungi_q1", 12),
                                     rep("fungi_q2", 12)),
                           variable = names(p_step2) , 
                           p_val = round(p_step2, 3), 
                           p_val_adj = round(p_adj_step2,3))

##---------------------------------------------------------------
##                   Variance Partitioning                      -
##---------------------------------------------------------------
# For the ModEVA function we need to create models for all of the different combinations of 
# abiotic, biotic and geographic variables. 
# A = biotic; B = abiotic, C = exploratory, AB = biotic + abiotic
# AC = biotic + exploratory, BC = abiotic + exploratory 

###
# Algae
###
#q0 
lm_q0_alg_second_biotic <- stats::lm(alg_q0 ~ bac_q0 +
                                       fun_q0 + 
                                       offset(library_size_alg),
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q0_alg_second_abiotic <- stats::lm(alg_q0 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest + 
                                        offset(library_size_alg),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q0_alg_second_geo <- stats::lm(alg_q0 ~ exploratory + 
                                    offset(library_size_alg),
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q0_alg_second_geo_bio <- stats::lm(alg_q0 ~ exploratory +
                                        bac_q0 +
                                        fun_q0 + 
                                        offset(library_size_alg),
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q0_alg_second_geo_abio <- stats::lm(alg_q0 ~ exploratory +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_alg),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q0_alg_second_bio_abio <- stats::lm(alg_q0 ~ bac_q0 +
                                         fun_q0 +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_alg),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

alg_var_lm_q0 <- modEvA::varPart(A = summary(lm_q0_alg_second_biotic)$r.squared,
                B = summary(lm_q0_alg_second_abiotic)$r.squared,
                C = summary(lm_q0_alg_second_geo)$r.squared,
                AB = summary(lm_q0_alg_second_bio_abio)$r.squared,
                AC = summary(lm_q0_alg_second_geo_bio)$r.squared,
                BC = summary(lm_q0_alg_second_geo_abio)$r.squared,
                ABC = summary(lm_q0_alg_second)$r.squared,
                A.name = "Biotic",
                B.name = "Abiotic",
                C.name = "Exploratory") %>% 
  mutate(across(Proportion, round, 3)) %>% 
  tibble::rownames_to_column(var = "variable")

#q1 
lm_q1_alg_second_biotic <- stats::lm(alg_q1 ~ bac_q1 +
                                       fun_q1 + 
                                       offset(library_size_alg),
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q1_alg_second_abiotic <- stats::lm(alg_q1 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest + 
                                        offset(library_size_alg),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q1_alg_second_geo <- stats::lm(alg_q1 ~ exploratory + 
                                    offset(library_size_alg),
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q1_alg_second_geo_bio <- stats::lm(alg_q1 ~ exploratory +
                                        bac_q1 +
                                        fun_q1 + 
                                        offset(library_size_alg),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q1_alg_second_geo_abio <- stats::lm(alg_q1 ~ exploratory +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_alg),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

lm_q1_alg_second_bio_abio <- stats::lm(alg_q1 ~ bac_q1 +
                                         fun_q1 +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_alg),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

alg_var_lm_q1 <- modEvA::varPart(A = summary(lm_q1_alg_second_biotic)$r.squared,
                B = summary(lm_q1_alg_second_abiotic)$r.squared,
                C = summary(lm_q1_alg_second_geo)$r.squared,
                AB = summary(lm_q1_alg_second_bio_abio)$r.squared,
                AC = summary(lm_q1_alg_second_geo_bio)$r.squared,
                BC = summary(lm_q1_alg_second_geo_abio)$r.squared,
                ABC = summary(lm_q1_alg_second)$r.squared,
                A.name = "Biotic",
                B.name = "Abiotic",
                C.name = "Exploratory") %>% 
  mutate(across(Proportion, round, 3))%>% 
  tibble::rownames_to_column(var = "variable")

#q2 
lm_q2_alg_second_biotic <- stats::lm(alg_q2 ~ bac_q2 +
                                       fun_q2 + 
                                       offset(library_size_alg),
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q2_alg_second_abiotic <- stats::lm(alg_q2 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest + 
                                        offset(library_size_alg),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q2_alg_second_geo <- stats::lm(alg_q2 ~ exploratory + 
                                    offset(library_size_alg),
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q2_alg_second_geo_bio <- stats::lm(alg_q2 ~ exploratory +
                                        bac_q2 +
                                        fun_q2 + 
                                        offset(library_size_alg),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q2_alg_second_geo_abio <- stats::lm(alg_q2 ~ exploratory +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_alg),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

lm_q2_alg_second_bio_abio <- stats::lm(alg_q2 ~ bac_q2 +
                                         fun_q2 +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_alg),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

alg_var_lm_q2 <- modEvA::varPart(A = summary(lm_q2_alg_second_biotic)$r.squared,
                B = summary(lm_q2_alg_second_abiotic)$r.squared,
                C = summary(lm_q2_alg_second_geo)$r.squared,
                AB = summary(lm_q2_alg_second_bio_abio)$r.squared,
                AC = summary(lm_q2_alg_second_geo_bio)$r.squared,
                BC = summary(lm_q2_alg_second_geo_abio)$r.squared,
                ABC = summary(lm_q2_alg_second)$r.squared,
                A.name = "Biotic",
                B.name = "Abiotic",
                C.name = "Exploratory")  %>% 
  dplyr::mutate(across(Proportion, round, 3)) %>% 
  tibble::rownames_to_column(var = "variable")

# Combine in one table 
alg_variance_lm <- rbind(alg_var_lm_q0, alg_var_lm_q1, alg_var_lm_q2) %>% 
  cbind(q_lev = c(rep("q0", 8), rep("q1", 8), rep("q2", 8))) %>% 
  mutate(variable = replace(variable, variable == "Abiotic", "abiotic (a)")) %>% 
  mutate(variable = replace(variable, variable == "Biotic", "biotic (b)")) %>% 
  mutate(variable = replace(variable, variable == "Exploratory", "geographic (g)"))  %>% 
  mutate(variable = replace(variable, variable == "Unexplained", "unexplained"))  %>% 
  mutate(variable = replace(variable, variable == "Biotic_Abiotic", "b+a")) %>% 
  mutate(variable = replace(variable, variable == "Biotic_Exploratory", "b+g")) %>% 
  mutate(variable = replace(variable, variable == "Abiotic_Exploratory", "a+g")) %>% 
  mutate(variable = replace(variable, variable == "Biotic_Abiotic_Exploratory", "a+b+g")) %>% 
  mutate(variable = factor(variable, levels = c("unexplained","a+b+g", "a+g", "b+g", "b+a",
                                                        "geographic (g)", "abiotic (a)", "biotic (b)")))  %>% 
  mutate(q_lev = factor(q_lev, levels = c("q2", "q1", "q0")))

###
# Fungi
###
#q0 
lm_q0_fun_second_biotic <- stats::lm(fun_q0 ~ bac_q0 +
                                       alg_q0 + 
                                       offset(library_size_fun),
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q0_fun_second_abiotic <- stats::lm(fun_q0 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest + 
                                        offset(library_size_fun),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q0_fun_second_geo <- stats::lm(fun_q0 ~ exploratory + 
                                    offset(library_size_fun),
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q0_fun_second_geo_bio <- stats::lm(fun_q0 ~ exploratory +
                                        bac_q0 +
                                        alg_q0 + 
                                        offset(library_size_fun),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q0_fun_second_geo_abio <- stats::lm(fun_q0 ~ exploratory +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_fun),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

lm_q0_fun_second_bio_abio <- stats::lm(fun_q0 ~ bac_q0 +
                                         alg_q0 +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_fun),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

fun_var_lm_q0 <- modEvA::varPart(A = summary(lm_q0_fun_second_biotic)$r.squared,
                B = summary(lm_q0_fun_second_abiotic)$r.squared,
                C = summary(lm_q0_fun_second_geo)$r.squared,
                AB = summary(lm_q0_fun_second_bio_abio)$r.squared,
                AC = summary(lm_q0_fun_second_geo_bio)$r.squared,
                BC = summary(lm_q0_fun_second_geo_abio)$r.squared,
                ABC = summary(lm_q0_fun_second)$r.squared,
                A.name = "Biotic",
                B.name = "Abiotic",
                C.name = "Exploratory") %>% 
  mutate(across(Proportion, round, 3)) %>% 
  tibble::rownames_to_column(var = "variable")

#q1 
lm_q1_fun_second_biotic <- stats::lm(fun_q1 ~ bac_q1 +
                                       alg_q1 + 
                                       offset(library_size_fun),
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q1_fun_second_abiotic <- stats::lm(fun_q1 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest + 
                                        offset(library_size_fun),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q1_fun_second_geo <- stats::lm(fun_q1 ~ exploratory + 
                                    offset(library_size_fun),
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q1_fun_second_geo_bio <- stats::lm(fun_q1 ~ exploratory +
                                        bac_q1 +
                                        alg_q1 + 
                                        offset(library_size_fun),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q1_fun_second_geo_abio <- stats::lm(fun_q1 ~ exploratory +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_fun),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

lm_q1_fun_second_bio_abio <- stats::lm(fun_q1 ~ bac_q1 +
                                         alg_q1 +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_fun),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

fun_var_lm_q1 <- modEvA::varPart(A = summary(lm_q1_fun_second_biotic)$r.squared,
                B = summary(lm_q1_fun_second_abiotic)$r.squared,
                C = summary(lm_q1_fun_second_geo)$r.squared,
                AB = summary(lm_q1_fun_second_bio_abio)$r.squared,
                AC = summary(lm_q1_fun_second_geo_bio)$r.squared,
                BC = summary(lm_q1_fun_second_geo_abio)$r.squared,
                ABC = summary(lm_q1_fun_second)$r.squared,
                A.name = "Biotic",
                B.name = "Abiotic",
                C.name = "Exploratory") %>% 
  mutate(across(Proportion, round, 3)) %>% 
  tibble::rownames_to_column(var = "variable")

#q2 
lm_q2_fun_second_biotic <- stats::lm(fun_q2 ~ bac_q2 +
                                       alg_q2 + 
                                       offset(library_size_fun),
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q2_fun_second_abiotic <- stats::lm(fun_q2 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest + 
                                        offset(library_size_fun),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q2_fun_second_geo <- stats::lm(fun_q2 ~ exploratory + 
                                    offset(library_size_fun),
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q2_fun_second_geo_bio <- stats::lm(fun_q2 ~ exploratory +
                                        bac_q2 +
                                        alg_q2 + 
                                        offset(library_size_fun),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q2_fun_second_geo_abio <- stats::lm(fun_q2 ~ exploratory +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_fun),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

lm_q2_fun_second_bio_abio <- stats::lm(fun_q2 ~ bac_q2 +
                                         alg_q2 +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_fun),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

fun_var_lm_q2 <- modEvA::varPart(A = summary(lm_q2_fun_second_biotic)$r.squared,
                B = summary(lm_q2_fun_second_abiotic)$r.squared,
                C = summary(lm_q2_fun_second_geo)$r.squared,
                AB = summary(lm_q2_fun_second_bio_abio)$r.squared,
                AC = summary(lm_q2_fun_second_geo_bio)$r.squared,
                BC = summary(lm_q2_fun_second_geo_abio)$r.squared,
                ABC = summary(lm_q2_fun_second)$r.squared,
                A.name = "Biotic",
                B.name = "Abiotic",
                C.name = "Exploratory") %>% 
  mutate(across(Proportion, round, 3)) %>% 
  tibble::rownames_to_column(var = "variable")

fun_variance_lm <- rbind(fun_var_lm_q0, fun_var_lm_q1, fun_var_lm_q2) %>% 
  cbind(q_lev = c(rep("q0", 8), rep("q1", 8), rep("q2", 8))) %>% 
  mutate(variable = replace(variable, variable == "Abiotic", "abiotic (a)")) %>% 
  mutate(variable = replace(variable, variable == "Biotic", "biotic (b)")) %>% 
  mutate(variable = replace(variable, variable == "Exploratory", "geographic (g)"))  %>% 
  mutate(variable = replace(variable, variable == "Unexplained", "unexplained"))  %>% 
  mutate(variable = replace(variable, variable == "Biotic_Abiotic", "b+a")) %>% 
  mutate(variable = replace(variable, variable == "Biotic_Exploratory", "b+g")) %>% 
  mutate(variable = replace(variable, variable == "Abiotic_Exploratory", "a+g")) %>% 
  mutate(variable = replace(variable, variable == "Biotic_Abiotic_Exploratory", "a+b+g")) %>% 
  mutate(variable = factor(variable, levels = c("unexplained","a+b+g", "a+g", "b+g", "b+a",
                                                "geographic (g)", "abiotic (a)", "biotic (b)")))  %>% 
  mutate(q_lev = factor(q_lev, levels = c("q2", "q1", "q0")))

###
# Bacteria
###
#q0 
lm_q0_bac_second_biotic <- stats::lm(bac_q0 ~ fun_q0 +
                                       alg_q0 + 
                                       offset(library_size_bac),
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q0_bac_second_abiotic <- stats::lm(bac_q0 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest + 
                                        offset(library_size_bac),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q0_bac_second_geo <- stats::lm(bac_q0 ~ exploratory + 
                                    offset(library_size_bac),
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q0_bac_second_geo_bio <- stats::lm(bac_q0 ~ exploratory +
                                        fun_q0 +
                                        alg_q0 + 
                                        offset(library_size_bac),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q0_bac_second_geo_abio <- stats::lm(bac_q0 ~ exploratory +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_bac),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

lm_q0_bac_second_bio_abio <- stats::lm(bac_q0 ~ fun_q0 +
                                         alg_q0 +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_bac),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

bac_var_lm_q0 <- modEvA::varPart(A = summary(lm_q0_bac_second_biotic)$r.squared,
                B = summary(lm_q0_bac_second_abiotic)$r.squared,
                C = summary(lm_q0_bac_second_geo)$r.squared,
                AB = summary(lm_q0_bac_second_bio_abio)$r.squared,
                AC = summary(lm_q0_bac_second_geo_bio)$r.squared,
                BC = summary(lm_q0_bac_second_geo_abio)$r.squared,
                ABC = summary(lm_q0_bac_second)$r.squared,
                A.name = "Biotic",
                B.name = "Abiotic",
                C.name = "Exploratory") %>% 
  mutate(across(Proportion, round, 3)) %>% 
  tibble::rownames_to_column(var = "variable")

#q1 
lm_q1_bac_second_biotic <- stats::lm(bac_q1 ~ fun_q1 +
                                       alg_q1 + 
                                       offset(library_size_bac),
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q1_bac_second_abiotic <- stats::lm(bac_q1 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest + 
                                        offset(library_size_bac),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q1_bac_second_geo <- stats::lm(bac_q1 ~ exploratory + 
                                    offset(library_size_bac),
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q1_bac_second_geo_bio <- stats::lm(bac_q1 ~ exploratory +
                                        fun_q1 +
                                        alg_q1 + 
                                        offset(library_size_bac),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q1_bac_second_geo_abio <- stats::lm(bac_q1 ~ exploratory +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_bac),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

lm_q1_bac_second_bio_abio <- stats::lm(bac_q1 ~ fun_q1 +
                                         alg_q1 +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_bac),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

bac_var_lm_q1 <- modEvA::varPart(A = summary(lm_q1_bac_second_biotic)$r.squared,
                B = summary(lm_q1_bac_second_abiotic)$r.squared,
                C = summary(lm_q1_bac_second_geo)$r.squared,
                AB = summary(lm_q1_bac_second_bio_abio)$r.squared,
                AC = summary(lm_q1_bac_second_geo_bio)$r.squared,
                BC = summary(lm_q1_bac_second_geo_abio)$r.squared,
                ABC = summary(lm_q1_bac_second)$r.squared,
                A.name = "Biotic",
                B.name = "Abiotic",
                C.name = "Exploratory") %>% 
  mutate(across(Proportion, round, 3)) %>% 
  tibble::rownames_to_column(var = "variable")

#q2 
lm_q2_bac_second_biotic <- stats::lm(bac_q2 ~ fun_q2 +
                                       alg_q2 + 
                                       offset(library_size_bac),
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q2_bac_second_abiotic <- stats::lm(bac_q2 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest + 
                                        offset(library_size_bac),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q2_bac_second_geo <- stats::lm(bac_q2 ~ exploratory + 
                                    offset(library_size_bac),
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q2_bac_second_geo_bio <- stats::lm(bac_q2 ~ exploratory +
                                        fun_q2 +
                                        alg_q2 + 
                                        offset(library_size_bac),
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q2_bac_second_geo_abio <- stats::lm(bac_q2 ~ exploratory +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_bac),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

lm_q2_bac_second_bio_abio <- stats::lm(bac_q2 ~ fun_q2 +
                                         alg_q2 +
                                         rH_200 + 
                                         Ta_200 +
                                         stand_density_basal_area +
                                         DBH_avg +
                                         d_gini +
                                         canopy_openness_2019 +
                                         dom_tot_ratio + 
                                         RA_forest + 
                                         offset(library_size_bac),
                                       data = metadata_alpha_scaled,
                                       na.action = "na.fail")

bac_var_lm_q2 <- modEvA::varPart(A = summary(lm_q2_bac_second_biotic)$r.squared,
                B = summary(lm_q2_bac_second_abiotic)$r.squared,
                C = summary(lm_q2_bac_second_geo)$r.squared,
                AB = summary(lm_q2_bac_second_bio_abio)$r.squared,
                AC = summary(lm_q2_bac_second_geo_bio)$r.squared,
                BC = summary(lm_q2_bac_second_geo_abio)$r.squared,
                ABC = summary(lm_q2_bac_second)$r.squared,
                A.name = "Biotic",
                B.name = "Abiotic",
                C.name = "Exploratory") %>% 
  mutate(across(Proportion, round, 3)) %>% 
  tibble::rownames_to_column(var = "variable")

bac_variance_lm <- rbind(bac_var_lm_q0, bac_var_lm_q1, bac_var_lm_q2) %>% 
  cbind(q_lev = c(rep("q0", 8), rep("q1", 8), rep("q2", 8))) %>% 
  mutate(variable = replace(variable, variable == "Abiotic", "abiotic (a)")) %>% 
  mutate(variable = replace(variable, variable == "Biotic", "biotic (b)")) %>% 
  mutate(variable = replace(variable, variable == "Exploratory", "geographic (g)"))  %>% 
  mutate(variable = replace(variable, variable == "Unexplained", "unexplained"))  %>% 
  mutate(variable = replace(variable, variable == "Biotic_Abiotic", "b+a")) %>% 
  mutate(variable = replace(variable, variable == "Biotic_Exploratory", "b+g")) %>% 
  mutate(variable = replace(variable, variable == "Abiotic_Exploratory", "a+g")) %>% 
  mutate(variable = replace(variable, variable == "Biotic_Abiotic_Exploratory", "a+b+g")) %>% 
  mutate(variable = factor(variable, levels = c("unexplained","a+b+g", "a+g", "b+g", "b+a",
                                                "geographic (g)", "abiotic (a)", "biotic (b)")))  %>% 
  mutate(q_lev = factor(q_lev, levels = c("q2", "q1", "q0")))

# Save the tables because we will create a combined figure with the variance partitioning 
# of beta diversity. 
base::saveRDS(bac_variance_lm, here("Data", "bac_variance_lm.rds"))
base::saveRDS(alg_variance_lm, here("Data", "alg_variance_lm.rds"))
base::saveRDS(fun_variance_lm, here("Data", "fun_variance_lm.rds"))

#####
# Effect plots of alpha diversity
#####

# Colors

algae_col <- "#49BEAA"
bacteria_col <- "#e2ca20"
fungi_col <- "#ea594e"
      
# Effects on algae
alg_q0_pred_fun <- ggeffects::ggpredict(lm_q0_alg_second, terms = "fun_q0")
alg_q1_pred_fun <- ggeffects::ggpredict(lm_q1_alg_second, terms = "fun_q1")
alg_q2_pred_fun <- ggeffects::ggpredict(lm_q2_alg_second, terms = "fun_q2")

algae_fungal_effect <- ggplot() +
  geom_line(data = alg_q0_pred_fun,mapping = aes(x = x, y = predicted),
            color = algae_col, linewidth = 0.7) + 
  geom_text(data = alg_q0_pred_fun, aes(x = max(x) - 0.08,
                                        y = max(predicted) + 0.05), label = "***",
            size = 4) +
  geom_line(data = alg_q1_pred_fun,mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = alg_q2_pred_fun,mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(   
       x= "Fungal \u03B1-diversity")  +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
algae_fungal_effect

alg_q0_pred_bac <- ggeffects::ggpredict(lm_q0_alg_second, terms = "bac_q0")
alg_q1_pred_bac <- ggeffects::ggpredict(lm_q1_alg_second, terms = "bac_q1")
alg_q2_pred_bac <- ggeffects::ggpredict(lm_q2_alg_second, terms = "bac_q2")

algae_bacterial_effect <- ggplot() +
  geom_line(data = alg_q0_pred_bac, mapping = aes(x = x, y = predicted),
            color = algae_col, linewidth = 0.7) +
  geom_line(data = alg_q1_pred_bac, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = alg_q2_pred_bac, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(y= "Algal \u03B1-diversity",    
       x= "Bacterial \u03B1-diversity") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
algae_bacterial_effect

alg_q0_pred_rH_200 <- ggeffects::ggpredict(lm_q0_alg_second, terms = "rH_200")
alg_q1_pred_rH_200 <- ggeffects::ggpredict(lm_q1_alg_second, terms = "rH_200")
alg_q2_pred_rH_200 <- ggeffects::ggpredict(lm_q2_alg_second, terms = "rH_200")

algae_humidity_effect <- ggplot() +
  geom_line(data = alg_q0_pred_rH_200, mapping = aes(x = x, y = predicted),
            color = algae_col, linewidth = 0.7) +
  geom_line(data = alg_q1_pred_rH_200, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = alg_q2_pred_rH_200, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Relative humidity") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
algae_humidity_effect

alg_q0_pred_Ta_200 <- ggeffects::ggpredict(lm_q0_alg_second, terms = "Ta_200")
alg_q1_pred_Ta_200 <- ggeffects::ggpredict(lm_q1_alg_second, terms = "Ta_200")
alg_q2_pred_Ta_200 <- ggeffects::ggpredict(lm_q2_alg_second, terms = "Ta_200")

algae_temperature_effect <- ggplot() +
  geom_line(data = alg_q0_pred_Ta_200, mapping = aes(x = x, y = predicted),
            color = algae_col, linewidth = 0.7) +
  geom_line(data = alg_q1_pred_Ta_200, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = alg_q2_pred_Ta_200, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(y = "\u03B1-diversity",    
       x= "Temperature") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
algae_temperature_effect

alg_q0_pred_stand_density_basal_area <- ggeffects::ggpredict(lm_q0_alg_second, terms = "stand_density_basal_area")
alg_q1_pred_stand_density_basal_area <- ggeffects::ggpredict(lm_q1_alg_second, terms = "stand_density_basal_area")
alg_q2_pred_stand_density_basal_area <- ggeffects::ggpredict(lm_q2_alg_second, terms = "stand_density_basal_area")

algae_stand_density_effect <- ggplot() +
  geom_line(data = alg_q0_pred_stand_density_basal_area, mapping = aes(x = x, y = predicted),
            color = algae_col, linewidth = 0.7) +
  geom_line(data = alg_q1_pred_stand_density_basal_area, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = alg_q2_pred_stand_density_basal_area, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotted", linewidth = 0.7)+ 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(y = "\u03B1-diversity",    
       x= "Stand Density")  +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
algae_stand_density_effect

alg_q0_pred_DBH_avg <- ggeffects::ggpredict(lm_q0_alg_second, terms = "DBH_avg")
alg_q1_pred_DBH_avg <- ggeffects::ggpredict(lm_q1_alg_second, terms = "DBH_avg")
alg_q2_pred_DBH_avg <- ggeffects::ggpredict(lm_q2_alg_second, terms = "DBH_avg")

algae_DBH_effect <- ggplot() +
  geom_line(data = alg_q0_pred_DBH_avg, mapping = aes(x = x, y = predicted),
            color = algae_col, linewidth = 0.7) +
  geom_line(data = alg_q1_pred_DBH_avg, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = alg_q2_pred_DBH_avg, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(y = "\u03B1-diversity",    
       x= "Average DBH") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
algae_DBH_effect

alg_q0_pred_d_gini <- ggeffects::ggpredict(lm_q0_alg_second, terms = "d_gini")
alg_q1_pred_d_gini <- ggeffects::ggpredict(lm_q1_alg_second, terms = "d_gini")
alg_q2_pred_d_gini <- ggeffects::ggpredict(lm_q2_alg_second, terms = "d_gini")

algae_gini_effect <- ggplot() +
  geom_line(data = alg_q0_pred_d_gini, mapping = aes(x = x, y = predicted),
            color = algae_col, linewidth = 0.7) +
  geom_line(data = alg_q1_pred_d_gini, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = alg_q2_pred_d_gini, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(y = "\u03B1-diversity",    
       x= "Gini Coefficient") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
algae_gini_effect

alg_q0_pred_canopy_openness_2019 <- ggeffects::ggpredict(lm_q0_alg_second, terms = "canopy_openness_2019")
alg_q1_pred_canopy_openness_2019 <- ggeffects::ggpredict(lm_q1_alg_second, terms = "canopy_openness_2019")
alg_q2_pred_canopy_openness_2019 <- ggeffects::ggpredict(lm_q2_alg_second, terms = "canopy_openness_2019")

algae_canopy_effect <- ggplot() +
  geom_line(data = alg_q0_pred_canopy_openness_2019, mapping = aes(x = x, y = predicted),
            color = algae_col, linewidth = 0.7) +
  geom_line(data = alg_q1_pred_canopy_openness_2019, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = alg_q2_pred_canopy_openness_2019, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Canopy openness") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
algae_canopy_effect

alg_q0_pred_dom_tot_ratio <- ggeffects::ggpredict(lm_q0_alg_second, terms = "dom_tot_ratio")
alg_q1_pred_dom_tot_ratio <- ggeffects::ggpredict(lm_q1_alg_second, terms = "dom_tot_ratio")
alg_q2_pred_dom_tot_ratio <- ggeffects::ggpredict(lm_q2_alg_second, terms = "dom_tot_ratio")

algae_dom_tot_ratio_effect <- ggplot() +
  geom_line(data = alg_q0_pred_dom_tot_ratio, mapping = aes(x = x, y = predicted),
            color = algae_col, linewidth = 0.7) +
  geom_line(data = alg_q1_pred_dom_tot_ratio, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = alg_q2_pred_dom_tot_ratio, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(y = "\u03B1-diversity",    
       x= "Dominant Trees:Total Trees") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
algae_dom_tot_ratio_effect

alg_q0_pred_RA_forest <- ggeffects::ggpredict(lm_q0_alg_second, terms = "RA_forest")
alg_q1_pred_RA_forest <- ggeffects::ggpredict(lm_q1_alg_second, terms = "RA_forest")
alg_q2_pred_RA_forest <- ggeffects::ggpredict(lm_q2_alg_second, terms = "RA_forest")

algae_forest_effect <- ggplot() +
  geom_line(data = alg_q0_pred_RA_forest, mapping = aes(x = x, y = predicted),
            color = algae_col, linewidth = 0.7) +
  geom_line(data = alg_q1_pred_RA_forest, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = alg_q2_pred_RA_forest, mapping = aes(x = x, y = predicted),
            color = algae_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(y = "\u03B1-diversity",   
       x= "Surrounding Forest Area") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
algae_forest_effect

###
#Region effect
###

alg_q0_pred_region <- ggeffects::ggpredict(lm_q0_alg_second, terms = "exploratory")
alg_q1_pred_region <- ggeffects::ggpredict(lm_q1_alg_second, terms = "exploratory")
alg_q2_pred_region <- ggeffects::ggpredict(lm_q2_alg_second, terms = "exploratory")

algae_region_effect <- ggplot2::ggplot() +
  ggplot2::geom_boxplot(data = alg_q0_pred_region, mapping = aes(x = x, y = predicted),
                        color = algae_col) +
  ggplot2::geom_boxplot(data = alg_q1_pred_region, mapping = aes(x = x, y = predicted),
                        color = algae_col, linetype = "dotdash", linewidth = 0.7) +
  ggplot2:: geom_boxplot(data = alg_q2_pred_region, mapping = aes(x = x, y = predicted),
                         color = algae_col, linetype = "dotted", linewidth = 0.7) +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                 text = element_text(size = 7, family = "sans")) +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1)) +
  scale_x_discrete(labels = c("South-West", "Central", "North-East")) +
  labs(y = "\u03B1-diversity",
       x= "Region") 

algae_region_effect


###########
# Effects on bacteria
###########

bac_q0_pred_fun <- ggeffects::ggpredict(lm_q0_bac_second, terms = "fun_q0")
bac_q1_pred_fun <- ggeffects::ggpredict(lm_q1_bac_second, terms = "fun_q1")
bac_q2_pred_fun <- ggeffects::ggpredict(lm_q2_bac_second, terms = "fun_q2")

bacteria_fungal_effect <- ggplot() +
  geom_line(data = bac_q0_pred_fun,mapping = aes(x = x, y = predicted),
            color = bacteria_col, linewidth = 0.7) + 
  geom_text(data = bac_q0_pred_fun, aes(x = max(x) - 0.2 ,
                                        y = max(predicted) + 0.05), label = "***",
            size = 4) +
  geom_line(data = bac_q1_pred_fun,mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotdash", linewidth = 0.7) + 
  geom_text(data = bac_q1_pred_fun, aes(x = max(x) + 0.8,
                                        y = max(predicted) + 0.05), label = "***",
            size = 4) +
  geom_line(data = bac_q2_pred_fun,mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans")) +
  labs(y = "Bacterial \u03B1-diversity",    
    x= "Fungal \u03B1-diversity")+
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
bacteria_fungal_effect

bac_q0_pred_alg <- ggeffects::ggpredict(lm_q0_bac_second, terms = "alg_q0")
bac_q1_pred_alg <- ggeffects::ggpredict(lm_q1_bac_second, terms = "alg_q1")
bac_q2_pred_alg <- ggeffects::ggpredict(lm_q2_bac_second, terms = "alg_q2")

bacteria_algal_effect <- ggplot() +
  geom_line(data = bac_q0_pred_alg, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linewidth = 0.7) +
  geom_line(data = bac_q1_pred_alg, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = bac_q2_pred_alg, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotted", linewidth = 0.7) + 
  geom_text(data = bac_q2_pred_alg, aes(x = max(x),
                                        y = max(predicted)+ 0.02), label = "**",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Algal \u03B1-diversity") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
bacteria_algal_effect

bac_q0_pred_rH_200 <- ggeffects::ggpredict(lm_q0_bac_second, terms = "rH_200")
bac_q1_pred_rH_200 <- ggeffects::ggpredict(lm_q1_bac_second, terms = "rH_200")
bac_q2_pred_rH_200 <- ggeffects::ggpredict(lm_q2_bac_second, terms = "rH_200")

bacteria_humidity_effect <- ggplot() +
  geom_line(data = bac_q0_pred_rH_200, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linewidth = 0.7) +
  geom_line(data = bac_q1_pred_rH_200, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = bac_q2_pred_rH_200, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(   
       x= "Relative humidity") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
bacteria_humidity_effect

bac_q0_pred_Ta_200 <- ggeffects::ggpredict(lm_q0_bac_second, terms = "Ta_200")
bac_q1_pred_Ta_200 <- ggeffects::ggpredict(lm_q1_bac_second, terms = "Ta_200")
bac_q2_pred_Ta_200 <- ggeffects::ggpredict(lm_q2_bac_second, terms = "Ta_200")

bacteria_temperature_effect <- ggplot() +
  geom_line(data = bac_q0_pred_Ta_200, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linewidth = 0.7) +
  geom_line(data = bac_q1_pred_Ta_200, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = bac_q2_pred_Ta_200, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotted", linewidth = 0.7)  + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Temperature") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
bacteria_temperature_effect

bac_q0_pred_stand_density_basal_area <- ggeffects::ggpredict(lm_q0_bac_second, terms = "stand_density_basal_area")
bac_q1_pred_stand_density_basal_area <- ggeffects::ggpredict(lm_q1_bac_second, terms = "stand_density_basal_area")
bac_q2_pred_stand_density_basal_area <- ggeffects::ggpredict(lm_q2_bac_second, terms = "stand_density_basal_area")

bacteria_stand_density_effect <- ggplot() +
  geom_line(data = bac_q0_pred_stand_density_basal_area, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linewidth = 0.7) +
  geom_line(data = bac_q1_pred_stand_density_basal_area, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = bac_q2_pred_stand_density_basal_area, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotted", linewidth = 0.7)+ 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Stand Density")  +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
bacteria_stand_density_effect

bac_q0_pred_DBH_avg <- ggeffects::ggpredict(lm_q0_bac_second, terms = "DBH_avg")
bac_q1_pred_DBH_avg <- ggeffects::ggpredict(lm_q1_bac_second, terms = "DBH_avg")
bac_q2_pred_DBH_avg <- ggeffects::ggpredict(lm_q2_bac_second, terms = "DBH_avg")

bacteria_DBH_effect <- ggplot() +
  geom_line(data = bac_q0_pred_DBH_avg, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linewidth = 0.7) +
  geom_line(data = bac_q1_pred_DBH_avg, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = bac_q2_pred_DBH_avg, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotted", linewidth = 0.7) + 
  geom_text(data = bac_q2_pred_DBH_avg, aes(x = max(x) + 0.1,
                                           y = max(predicted) - 0.04), label = "*",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Average DBH") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
bacteria_DBH_effect

bac_q0_pred_d_gini <- ggeffects::ggpredict(lm_q0_bac_second, terms = "d_gini")
bac_q1_pred_d_gini <- ggeffects::ggpredict(lm_q1_bac_second, terms = "d_gini")
bac_q2_pred_d_gini <- ggeffects::ggpredict(lm_q2_bac_second, terms = "d_gini")

bacteria_gini_effect <- ggplot() +
  geom_line(data = bac_q0_pred_d_gini, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linewidth = 0.7) + 
  geom_text(data = bac_q0_pred_d_gini, aes(x = max(x) + 0.2,
                                        y = max(predicted) - 0.02), label = "***",
            size = 4)+
  geom_line(data = bac_q1_pred_d_gini, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotdash", linewidth = 0.7) + 
  geom_text(data = bac_q1_pred_d_gini, aes(x = max(x) + 0.2,
                                        y = max(predicted) - 0.02), label = "***",
            size = 4)+
  geom_line(data = bac_q2_pred_d_gini, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotted", linewidth = 0.7) + 
  geom_text(data = bac_q2_pred_d_gini, aes(x = max(x) + 0.2,
                                        y = max(predicted) - 0.05), label = "***",
            size = 4)+ 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Gini Coefficient") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
bacteria_gini_effect

bac_q0_pred_canopy_openness_2019 <- ggeffects::ggpredict(lm_q0_bac_second, terms = "canopy_openness_2019")
bac_q1_pred_canopy_openness_2019 <- ggeffects::ggpredict(lm_q1_bac_second, terms = "canopy_openness_2019")
bac_q2_pred_canopy_openness_2019 <- ggeffects::ggpredict(lm_q2_bac_second, terms = "canopy_openness_2019")

bacteria_canopy_effect <- ggplot() +
  geom_line(data = bac_q0_pred_canopy_openness_2019, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linewidth = 0.7) + 
  geom_text(data = bac_q0_pred_canopy_openness_2019, aes(x = max(x) + 0.15,
                                           y = min(predicted) + 0.1), label = "***",
            size = 4) +
  geom_line(data = bac_q1_pred_canopy_openness_2019, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotdash", linewidth = 0.7) + 
  geom_text(data = bac_q1_pred_canopy_openness_2019, aes(x = max(x) + 0.3,
                                           y = min(predicted) + 0.07), label = "*",
            size = 4) +
  geom_line(data = bac_q2_pred_canopy_openness_2019, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Canopy openness") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
bacteria_canopy_effect

bac_q0_pred_dom_tot_ratio <- ggeffects::ggpredict(lm_q0_bac_second, terms = "dom_tot_ratio")
bac_q1_pred_dom_tot_ratio <- ggeffects::ggpredict(lm_q1_bac_second, terms = "dom_tot_ratio")
bac_q2_pred_dom_tot_ratio <- ggeffects::ggpredict(lm_q2_bac_second, terms = "dom_tot_ratio")

bacteria_dom_tot_ratio_effect <- ggplot() +
  geom_line(data = bac_q0_pred_dom_tot_ratio, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linewidth = 0.7) +
  geom_line(data = bac_q1_pred_dom_tot_ratio, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = bac_q2_pred_dom_tot_ratio, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Dominant Trees:Total Trees") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
bacteria_dom_tot_ratio_effect

bac_q0_pred_RA_forest <- ggeffects::ggpredict(lm_q0_bac_second, terms = "RA_forest")
bac_q1_pred_RA_forest <- ggeffects::ggpredict(lm_q1_bac_second, terms = "RA_forest")
bac_q2_pred_RA_forest <- ggeffects::ggpredict(lm_q2_bac_second, terms = "RA_forest")

bacteria_forest_effect <- ggplot() +
  geom_line(data = bac_q0_pred_RA_forest, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linewidth = 0.7) +
  geom_line(data = bac_q1_pred_RA_forest, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = bac_q2_pred_RA_forest, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Surrounding Forest Area") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
bacteria_forest_effect

###
#Region effect
###

bac_q0_pred_region <- ggeffects::ggpredict(lm_q0_bac_second, terms = "exploratory")
bac_q1_pred_region <- ggeffects::ggpredict(lm_q1_bac_second, terms = "exploratory")
bac_q2_pred_region <- ggeffects::ggpredict(lm_q2_bac_second, terms = "exploratory")

bacteria_region_effect <- ggplot2::ggplot() +
  ggplot2::geom_boxplot(data = bac_q0_pred_region, mapping = aes(x = x, y = predicted),
                        color = bacteria_col) +
  ggplot2::geom_boxplot(data = bac_q1_pred_region, mapping = aes(x = x, y = predicted),
                        color = bacteria_col, linetype = "dotdash", linewidth = 0.7) +
  ggplot2:: geom_boxplot(data = bac_q2_pred_region, mapping = aes(x = x, y = predicted),
                         color = bacteria_col, linetype = "dotted", linewidth = 0.7) +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                 text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1)) +
  scale_x_discrete(labels = c("South-West", "Central", "North-East")) +
  labs(x= "Region") 

bacteria_region_effect

###########
# Effects on Fungi
###########

fun_q0_pred_bac <- ggeffects::ggpredict(lm_q0_fun_second, terms = "bac_q0")
fun_q1_pred_bac <- ggeffects::ggpredict(lm_q1_fun_second, terms = "bac_q1")
fun_q2_pred_bac <- ggeffects::ggpredict(lm_q2_fun_second, terms = "bac_q2")

fungi_bacterial_effect <- ggplot() +
  geom_line(data = fun_q0_pred_bac,mapping = aes(x = x, y = predicted),
            color = fungi_col, linewidth = 0.7) + 
  geom_text(data = fun_q0_pred_bac, aes(x = max(x) - 0.9,
                                        y = max(predicted)), label = "***",
            size = 4) +
  geom_line(data = fun_q1_pred_bac,mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = fun_q2_pred_bac,mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"),legend.key.width=unit(1.5,"cm")) +
  labs(y= "Fungal \u03B1-diversity",    
       x= "Bacterial \u03B1-diversity")   +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
fungi_bacterial_effect

fun_q0_pred_alg <- ggeffects::ggpredict(lm_q0_fun_second, terms = "alg_q0")
fun_q1_pred_alg <- ggeffects::ggpredict(lm_q1_fun_second, terms = "alg_q1")
fun_q2_pred_alg <- ggeffects::ggpredict(lm_q2_fun_second, terms = "alg_q2")

fungi_algal_effect <- ggplot() +
  geom_line(data = fun_q0_pred_alg, mapping = aes(x = x, y = predicted),
            color = fungi_col, linewidth = 0.7) + 
  geom_text(data = fun_q0_pred_alg, aes(x = max(x) + 0.35,
                                        y = max(predicted) - 0.1), label = "**",
            size = 4) +
  geom_line(data = fun_q1_pred_alg, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = fun_q2_pred_alg, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Algal \u03B1-diversity") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
fungi_algal_effect

fun_q0_pred_rH_200 <- ggeffects::ggpredict(lm_q0_fun_second, terms = "rH_200")
fun_q1_pred_rH_200 <- ggeffects::ggpredict(lm_q1_fun_second, terms = "rH_200")
fun_q2_pred_rH_200 <- ggeffects::ggpredict(lm_q2_fun_second, terms = "rH_200")

fungi_humidity_effect <- ggplot() +
  geom_line(data = fun_q0_pred_rH_200, mapping = aes(x = x, y = predicted),
            color = fungi_col, linewidth = 0.7) + 
  geom_text(data = fun_q0_pred_rH_200, aes(x = max(x) + 0.2,
                                        y = max(predicted) + 0.02), label = "**",
            size = 4) +
  geom_line(data = fun_q1_pred_rH_200, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotdash", linewidth = 0.7) + 
  geom_text(data = fun_q1_pred_rH_200, aes(x = max(x) + 0.2,
                                        y = max(predicted)), label = "**",
            size = 4) +
  geom_line(data = fun_q2_pred_rH_200, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotted", linewidth = 0.7) + 
  geom_text(data = fun_q2_pred_rH_200, aes(x = max(x) + 0.2,
                                        y = max(predicted) - 0.2), label = "**",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Relative humidity") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
fungi_humidity_effect

fun_q0_pred_Ta_200 <- ggeffects::ggpredict(lm_q0_fun_second, terms = "Ta_200")
fun_q1_pred_Ta_200 <- ggeffects::ggpredict(lm_q1_fun_second, terms = "Ta_200")
fun_q2_pred_Ta_200 <- ggeffects::ggpredict(lm_q2_fun_second, terms = "Ta_200")

fungi_temperature_effect <- ggplot() +
  geom_line(data = fun_q0_pred_Ta_200, mapping = aes(x = x, y = predicted),
            color = fungi_col, linewidth = 0.7) +
  geom_line(data = fun_q1_pred_Ta_200, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = fun_q2_pred_Ta_200, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Temperature") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
fungi_temperature_effect

fun_q0_pred_stand_density_basal_area <- ggeffects::ggpredict(lm_q0_fun_second, terms = "stand_density_basal_area")
fun_q1_pred_stand_density_basal_area <- ggeffects::ggpredict(lm_q1_fun_second, terms = "stand_density_basal_area")
fun_q2_pred_stand_density_basal_area <- ggeffects::ggpredict(lm_q2_fun_second, terms = "stand_density_basal_area")

fungi_stand_density_effect <- ggplot() +
  geom_line(data = fun_q0_pred_stand_density_basal_area, mapping = aes(x = x, y = predicted),
            color = fungi_col, linewidth = 0.7) +
  geom_line(data = fun_q1_pred_stand_density_basal_area, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = fun_q2_pred_stand_density_basal_area, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotted", linewidth = 0.7)+ 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Stand Density")  +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
fungi_stand_density_effect

fun_q0_pred_DBH_avg <- ggeffects::ggpredict(lm_q0_fun_second, terms = "DBH_avg")
fun_q1_pred_DBH_avg <- ggeffects::ggpredict(lm_q1_fun_second, terms = "DBH_avg")
fun_q2_pred_DBH_avg <- ggeffects::ggpredict(lm_q2_fun_second, terms = "DBH_avg")

fungi_DBH_effect <- ggplot() +
  geom_line(data = fun_q0_pred_DBH_avg, mapping = aes(x = x, y = predicted),
            color = fungi_col, linewidth = 0.7) +
  geom_line(data = fun_q1_pred_DBH_avg, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = fun_q2_pred_DBH_avg, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Average DBH") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
fungi_DBH_effect

fun_q0_pred_d_gini <- ggeffects::ggpredict(lm_q0_fun_second, terms = "d_gini")
fun_q1_pred_d_gini <- ggeffects::ggpredict(lm_q1_fun_second, terms = "d_gini")
fun_q2_pred_d_gini <- ggeffects::ggpredict(lm_q2_fun_second, terms = "d_gini")

fungi_gini_effect <- ggplot() +
  geom_line(data = fun_q0_pred_d_gini, mapping = aes(x = x, y = predicted),
            color = fungi_col, linewidth = 0.7) +
  geom_line(data = fun_q1_pred_d_gini, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = fun_q2_pred_d_gini, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Gini Coefficient") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
fungi_gini_effect

fun_q0_pred_canopy_openness_2019 <- ggeffects::ggpredict(lm_q0_fun_second, terms = "canopy_openness_2019")
fun_q1_pred_canopy_openness_2019 <- ggeffects::ggpredict(lm_q1_fun_second, terms = "canopy_openness_2019")
fun_q2_pred_canopy_openness_2019 <- ggeffects::ggpredict(lm_q2_fun_second, terms = "canopy_openness_2019")

fungi_canopy_effect <- ggplot() +
  geom_line(data = fun_q0_pred_canopy_openness_2019, mapping = aes(x = x, y = predicted),
            color = fungi_col, linewidth = 0.7) + 
  geom_text(data = fun_q0_pred_canopy_openness_2019, aes(x = max(x) + 0.2,
                                                         y = max(predicted) + 0.01), label = "**",
            size = 4) +
  geom_line(data = fun_q1_pred_canopy_openness_2019, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotdash", linewidth = 0.7) + 
  geom_text(data = fun_q1_pred_canopy_openness_2019, aes(x = max(x) + 0.2,
                                                         y = max(predicted) + 0.01), label = "**",
            size = 4) +
  geom_line(data = fun_q2_pred_canopy_openness_2019, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotted", linewidth = 0.7) + 
  geom_text(data = fun_q2_pred_canopy_openness_2019, aes(x = max(x) + 0.2,
                                                         y = max(predicted) - 0.2), label = "**",
            size = 4) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Canopy openness") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
fungi_canopy_effect

fun_q0_pred_dom_tot_ratio <- ggeffects::ggpredict(lm_q0_fun_second, terms = "dom_tot_ratio")
fun_q1_pred_dom_tot_ratio <- ggeffects::ggpredict(lm_q1_fun_second, terms = "dom_tot_ratio")
fun_q2_pred_dom_tot_ratio <- ggeffects::ggpredict(lm_q2_fun_second, terms = "dom_tot_ratio")

fungi_dom_tot_ratio_effect <- ggplot() +
  geom_line(data = fun_q0_pred_dom_tot_ratio, mapping = aes(x = x, y = predicted),
            color = fungi_col, linewidth = 0.7) +
  geom_line(data = fun_q1_pred_dom_tot_ratio, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = fun_q2_pred_dom_tot_ratio, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Dominant Trees:Total Trees") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
fungi_dom_tot_ratio_effect

fun_q0_pred_RA_forest <- ggeffects::ggpredict(lm_q0_fun_second, terms = "RA_forest")
fun_q1_pred_RA_forest <- ggeffects::ggpredict(lm_q1_fun_second, terms = "RA_forest")
fun_q2_pred_RA_forest <- ggeffects::ggpredict(lm_q2_fun_second, terms = "RA_forest")

fungi_forest_effect <- ggplot() +
  geom_line(data = fun_q0_pred_RA_forest, mapping = aes(x = x, y = predicted),
            color = fungi_col, linewidth = 0.7) +
  geom_line(data = fun_q1_pred_RA_forest, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotdash", linewidth = 0.7) +
  geom_line(data = fun_q2_pred_RA_forest, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotted", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Surrounding Forest Area") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
fungi_forest_effect

###
#Region effect
###

fun_q0_pred_region <- ggeffects::ggpredict(lm_q0_fun_second, terms = "exploratory")
fun_q1_pred_region <- ggeffects::ggpredict(lm_q1_fun_second, terms = "exploratory")
fun_q2_pred_region <- ggeffects::ggpredict(lm_q2_fun_second, terms = "exploratory")

fungi_region_effect <- ggplot2::ggplot() +
  ggplot2::geom_boxplot(data = fun_q0_pred_region, mapping = aes(x = x, y = predicted),
                        color = fungi_col) +
  ggplot2::geom_boxplot(data = fun_q1_pred_region, mapping = aes(x = x, y = predicted),
                        color = fungi_col, linetype = "dotdash", linewidth = 0.7) +
  ggplot2:: geom_boxplot(data = fun_q2_pred_region, mapping = aes(x = x, y = predicted),
                         color = fungi_col, linetype = "dotted", linewidth = 0.7) +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                 text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1)) +
  scale_x_discrete(labels = c("South-West", "Central", "North-East")) +
  labs(x= "Region") 

fungi_region_effect

############
# Plotting of the main text figures
############
# Grep a legend to append to the figure of the combined effect plots.
alpha_legend <- ggpubr::get_legend(ggplot() +
  geom_line(data = bac_q0_pred_fun,mapping = aes(x = x, y = predicted, linetype = "solid"),
            color = "black", linewidth = 0.7) + 
  geom_text(data = bac_q0_pred_fun, aes(x = max(x) + 0.1,
                                        y = max(predicted)), label = "*",
            size = 4) +
  geom_line(data = bac_q1_pred_fun,mapping = aes(x = x, y = predicted, linetype = "dotdash"),
            color = "black", linewidth = 0.7) +
  geom_line(data = bac_q2_pred_fun,mapping = aes(x = x, y = predicted, linetype = "dotted"),
            color = "black", linewidth = 0.7) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 7, family = "sans"), legend.position = "bottom") +
  scale_linetype_manual(name = element_text("Diversity Level"), 
                        values =c('solid',
                                  'dotdash', 
                                  'dotted'),
                        labels = c('q = 0 (all species)', 'q = 1 (typical)', 'q = 2 (dominant)')) +
  labs(y= "Bacterial \u03B1-diversity",    
       x= "Fungal \u03B1-diversity")  +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1)))

# Put a title for the biotic effects side. 
effects_algae1 <- ggpubr::ggarrange(algae_bacterial_effect, algae_fungal_effect,
                                    nrow = 1, ncol = 2) %>% 
  ggpubr::annotate_figure(top = ggpubr::text_grob("Biotic effects", size = 7, family = "sans", face = "bold"))

# Put a title for the abiotic effects side. 
effects_algae2 <- ggpubr::ggarrange(algae_canopy_effect, algae_humidity_effect, 
                                   nrow = 1, ncol = 2) %>% 
  ggpubr::annotate_figure(top = ggpubr::text_grob("Abiotic effects", size = 7, family = "sans", face = "bold"))

# Combine the algae abiotic and biotic sides. 
effects_algae3 <- ggpubr::ggarrange(effects_algae1,
                                    effects_algae2,
                                    nrow = 1, ncol = 2) 

# Fungal plots.
effects_fungi <- ggpubr::ggarrange(fungi_bacterial_effect, fungi_algal_effect, fungi_canopy_effect, fungi_humidity_effect, 
                                   nrow = 1, ncol = 4) %>% 
  ggpubr::annotate_figure(top = ggpubr::text_grob("",  color = fungi_col, face = "bold", size = 7)) 

# Bacterial plots.
effects_bacteria <- ggpubr::ggarrange(bacteria_fungal_effect, bacteria_algal_effect, bacteria_canopy_effect, bacteria_humidity_effect,
                                      nrow = 1, ncol = 4) %>% 
  ggpubr::annotate_figure(top = ggpubr::text_grob("",  color = bacteria_col, face = "bold", size = 7))

# Combine them in one big plot for the paper.
effects_alpha <- ggpubr::ggarrange(effects_algae3, effects_fungi, effects_bacteria,
                                   nrow = 3, common.legend = T, legend = "bottom",
                                   legend.grob = alpha_legend)
effects_alpha

# Save the plots as EPS and PNG. 
ggplot2::ggsave(here("Figures", "effects_alpha.eps"), plot = effects_alpha, device = cairo_ps,
       width = 175, height = 150, units = "mm", dpi = 600, bg = "white")

ggplot2::ggsave(here("Figures", "effects_alpha.png"), plot = effects_alpha, device = png,
       width = 175, height = 150, units = "mm")


#####
# Plotting supplementary alpha effect figures
#####
# Same procedure as above. 

###
#Algae
###

supplementary_alpha_algae <- ggpubr::ggarrange(algae_temperature_effect,
                                              algae_DBH_effect,
                                              algae_gini_effect,
                                              algae_stand_density_effect,
                                              algae_dom_tot_ratio_effect,
                                              algae_forest_effect,
                                              algae_region_effect,
                                              nrow = 7, ncol = 1)  %>% 
  ggpubr::annotate_figure(top = ggpubr::text_grob("Algae",  color = "black", face = "bold", size = 7))

supplementary_alpha_algae

###
#Fungi
###

supplementary_alpha_fungi <- ggpubr::ggarrange(fungi_temperature_effect,
                                              fungi_DBH_effect,
                                              fungi_gini_effect,
                                              fungi_stand_density_effect,
                                              fungi_dom_tot_ratio_effect,
                                              fungi_forest_effect,
                                              fungi_region_effect,
                                              nrow = 7, ncol = 1) %>% 
  ggpubr::annotate_figure(top = ggpubr::text_grob("Fungi",  color = "black", face = "bold", size = 7))

supplementary_alpha_fungi

###
#Bacteria
###

supplementary_alpha_bacteria <- ggpubr::ggarrange(bacteria_temperature_effect,
                                                 bacteria_DBH_effect,
                                                 bacteria_gini_effect,
                                                 bacteria_stand_density_effect,
                                                 bacteria_dom_tot_ratio_effect,
                                                 bacteria_forest_effect,
                                                 bacteria_region_effect,
                                                 nrow = 7, ncol = 1) %>% 
  ggpubr::annotate_figure(top = ggpubr::text_grob("Bacteria",  color = "black", face = "bold", size = 7))

supplementary_alpha_bacteria

# Final arrangement 

effects_alpha_supplementary <- ggpubr::ggarrange(supplementary_alpha_algae,
                                                supplementary_alpha_fungi,
                                                supplementary_alpha_bacteria,
                                                nrow = 1, ncol = 3, legend.grob = alpha_legend, legend = "bottom") 

effects_alpha_supplementary


ggplot2::ggsave(here("Figures", "effects_alpha_supplementary.eps"), plot = effects_alpha_supplementary, device = cairo_ps,
                width = 175, height = 250, units = "mm", dpi = 600, bg = "white")

ggplot2::ggsave(here("Figures","effects_alpha_supplementary.png"), plot = effects_alpha_supplementary, device = png,
                width = 175, height = 250, units = "mm")


#################################################################
##                          Section 4                          ##
##                   Miscallenous Data Export                  ##
#################################################################

# Export the metadata for the supplementary material. 
metadata_final <- metadata_alpha %>% 
  dplyr::mutate(exploratory = stringr::str_replace(exploratory, "Alb", "Swabian Alb"),
                exploratory = stringr::str_replace(exploratory, "Hainich", "Hainich-Dün"),
                exploratory = stringr::str_replace(exploratory, "Schorfheide", "Schorfheide-Chorin")) %>%
  dplyr::mutate(dominant_tree = stringr::str_replace(dominant_tree, "Picea_abies", "Picea abies"),
                dominant_tree = stringr::str_replace(dominant_tree, "Fagus_sylvatica", "Fagus sylvatica"),
                dominant_tree = stringr::str_replace(dominant_tree, "Pinus_sylvestris", "Pinus sylvestris")) %>%
  dplyr::rename(relative_humidity = rH_200,
                temperature = Ta_200,
                stand_density = stand_density_basal_area,
                average_DBH = DBH_avg,
                gini_coefficient = d_gini,
                forest_area = RA_forest,
                canopy_openness = canopy_openness_2019,
                ratio_of_dominant_trees = dom_tot_ratio,
                library_size_algae = library_size_alg,
                library_size_fungi = library_size_fun,
                library_size_bacteria = library_size_bac) %>% 
  dplyr::mutate(average_DBH = base::round(average_DBH, 2),
                canopy_openness = base::round(canopy_openness, 2),
                gini_coefficient = base::round(gini_coefficient, 2),
                stand_density = base::round(stand_density, 2),
                forest_area = base::round(forest_area, 2)) %>% 
  dplyr::select("Sample_ID", "exploratory", "dominant_tree","relative_humidity", "temperature",
                "average_DBH", "canopy_openness", "gini_coefficient", "stand_density",
                "ratio_of_dominant_trees", "forest_area", 
                "library_size_algae", "library_size_fungi", "library_size_bacteria")

write.csv(metadata_final, "metadata_final.csv", row.names = F)

# Export the table with coefficient estimates. 

estimates <- (summary(lm_q0_alg_second)$coefficients[-1, 1]) %>% 
  c(., (summary(lm_q1_alg_second)$coefficients[-1, 1])) %>% 
  c(., (summary(lm_q2_alg_second)$coefficients[-1, 1])) %>% 
  c(., (summary(lm_q0_bac_second)$coefficients[-1, 1])) %>%
  c(., (summary(lm_q1_bac_second)$coefficients[-1, 1])) %>%
  c(., (summary(lm_q2_bac_second)$coefficients[-1, 1])) %>%
  c(., (summary(lm_q0_fun_second)$coefficients[-1, 1])) %>%
  c(., (summary(lm_q1_fun_second)$coefficients[-1, 1])) %>%
  c(., (summary(lm_q2_fun_second)$coefficients[-1, 1]))

sta_err <- (summary(lm_q0_alg_second)$coefficients[-1, 2]) %>% 
  c(., (summary(lm_q1_alg_second)$coefficients[-1, 2])) %>% 
  c(., (summary(lm_q2_alg_second)$coefficients[-1, 2])) %>% 
  c(., (summary(lm_q0_bac_second)$coefficients[-1, 2])) %>%
  c(., (summary(lm_q1_bac_second)$coefficients[-1, 2])) %>%
  c(., (summary(lm_q2_bac_second)$coefficients[-1, 2])) %>%
  c(., (summary(lm_q0_fun_second)$coefficients[-1, 2])) %>%
  c(., (summary(lm_q1_fun_second)$coefficients[-1, 2])) %>%
  c(., (summary(lm_q2_fun_second)$coefficients[-1, 2]))
  
coefs <- data.frame(level = c(rep("algae_q0", 12),
                              rep("algae_q1", 12),
                              rep("algae_q2", 12),
                              rep("bacteria_q0", 12),
                              rep("bacteria_q1", 12),
                              rep("bacteria_q2", 12),
                              rep("fungi_q0", 12),
                              rep("fungi_q1", 12),
                              rep("fungi_q2", 12)),
                    variable = names(estimates),
                    estimate = round(estimates, 3),
                    standard_error = round(sta_err, 3)) %>% 
  dplyr::left_join(p_vals_step2, by = c("variable", "level"))

write.csv(coefs, "coef_table.csv", row.names = F)
