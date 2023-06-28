#################################################################
##                          Section 1                          ##
##                       Package Loading                       ##
#################################################################

library(here); packageVersion("here")
#1.0.1
library(effects); packageVersion("effects")
#4.2.1
library(hillR); packageVersion("hillR")
#0.5.1
library(tidyverse); packageVersion("tidyverse")
#1.3.2
library(phyloseq); packageVersion("phyloseq")
#1.40.0
library(ggdist); packageVersion("ggdist")
#3.1.1
library(gghalves); packageVersion("gghalves")
#0.1.3
library(ggeffects); packageVersion("ggeffects")
#1.1.2

#################################################################
##                          Section 2                          ##
##                        Data Loading                         ##
#################################################################

phy_algae_bark <- base::readRDS("phy_algae_bark.rds")

phy_bacteria_bark <- base::readRDS("phy_bacteria_bark.rds")

phy_fungi_bark <- base::readRDS("phy_fungi_bark.rds")

metadata_bark <- base::readRDS("metadata_bark.rds")

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

# Join with the rest of the metadata and rename the column to have distinct columns for the organisms.
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

# Scale the response variables to make effects comparable between our three organismal groups. 
metadata_alpha_scaled <- metadata_alpha %>% 
  dplyr::mutate(across(.cols = c("alg_q0", "alg_q1", "alg_q2",
                                 "bac_q0", "bac_q1", "bac_q2",
                                 "fun_q0", "fun_q1", "fun_q2"),~(scale(.) %>% as.vector))) %>% 
  # Furthermore, scale the explanatory variables since they have very different units. 
  dplyr::mutate(across(.cols = c("rH_200", "Ta_200", "stand_density_basal_area", 
                                 "DBH_avg", "d_gini", "RA_forest", 
                                 "canopy_openness_2019", "dom_tot_ratio"),~(scale(.) %>% as.vector)))

#####
# Model the alpha diversity in a first step.
#####

#ALGAE
lm_q0_alg <- stats::lm(alg_q0 ~ exploratory + 
                         dominant_tree + 
                         RA_forest + 
                         bac_q0 +
                         fun_q0,
                       data = metadata_alpha_scaled)
base::summary(lm_q0_alg)  
base::plot(effects::allEffects(lm_q0_alg))
base::plot(stats::residuals(lm_q0_alg), stats::fitted(lm_q0_alg))
stats::qqnorm(stats::residuals(lm_q0_alg))
stats::qqline(stats::residuals(lm_q0_alg))

lm_q1_alg <- stats::lm(alg_q1 ~ exploratory + 
                         dominant_tree + 
                         RA_forest + 
                         bac_q1 +
                         fun_q1,
                       data = metadata_alpha_scaled)
base::summary(lm_q1_alg)  
base::plot(allEffects(lm_q1_alg))
base::plot(stats::residuals(lm_q1_alg), stats::fitted(lm_q1_alg))
stats::qqnorm(stats::residuals(lm_q1_alg))
stats::qqline(stats::residuals(lm_q1_alg))

lm_q2_alg <- stats::lm(alg_q2 ~ exploratory + 
                         dominant_tree + 
                         RA_forest + 
                         bac_q2 +
                         fun_q2,
                       data = metadata_alpha_scaled)
base::summary(lm_q2_alg)  
base::plot(allEffects(lm_q2_alg))
base::plot(stats::residuals(lm_q2_alg), stats::fitted(lm_q2_alg))
stats::qqnorm(stats::residuals(lm_q2_alg))
stats::qqline(stats::residuals(lm_q2_alg))

# Bacteria
lm_q0_bac <- stats::lm(bac_q0 ~ exploratory + 
                         dominant_tree + 
                         RA_forest + 
                         alg_q0 +
                         fun_q0,
                       data = metadata_alpha_scaled)
base::summary(lm_q0_bac)  
base::plot(allEffects(lm_q0_bac))
base::plot(stats::residuals(lm_q0_bac), stats::fitted(lm_q0_bac))
stats::qqnorm(stats::residuals(lm_q0_bac))
stats::qqline(stats::residuals(lm_q0_bac))

lm_q1_bac <- stats::lm(bac_q1 ~ exploratory + 
                         dominant_tree + 
                         RA_forest + 
                         alg_q1 +
                         fun_q1,
                       data = metadata_alpha_scaled)
base::summary(lm_q1_bac)  
base::plot(allEffects(lm_q1_bac))
base::plot(stats::residuals(lm_q1_bac), stats::fitted(lm_q1_bac))
stats::qqnorm(stats::residuals(lm_q1_bac))
stats::qqline(stats::residuals(lm_q1_bac))

lm_q2_bac <- stats::lm(bac_q2 ~ exploratory + 
                         dominant_tree + 
                         RA_forest + 
                         alg_q2 +
                         fun_q2,
                       data = metadata_alpha_scaled)
base::summary(lm_q2_bac)  
base::plot(allEffects(lm_q2_bac))
base::plot(stats::residuals(lm_q2_bac), stats::fitted(lm_q2_bac))
stats::qqnorm(stats::residuals(lm_q2_bac))
stats::qqline(stats::residuals(lm_q2_bac))

# Fungi
lm_q0_fun <- stats::lm(fun_q0 ~ exploratory + 
                         dominant_tree + 
                         RA_forest + 
                         bac_q0 +
                         alg_q0,
                       data = metadata_alpha_scaled)
base::summary(lm_q0_fun)  
base::plot(allEffects(lm_q0_fun))
base::plot(stats::residuals(lm_q0_fun), stats::fitted(lm_q0_fun))
stats::qqnorm(stats::residuals(lm_q0_fun))
stats::qqline(stats::residuals(lm_q0_fun))

lm_q1_fun <- stats::lm(fun_q1 ~ exploratory + 
                         dominant_tree + 
                         RA_forest + 
                         bac_q1 +
                         alg_q1,
                       data = metadata_alpha_scaled)
base::summary(lm_q1_fun)  
base::plot(allEffects(lm_q1_fun))
base::plot(stats::residuals(lm_q1_fun), stats::fitted(lm_q1_fun))
stats::qqnorm(stats::residuals(lm_q1_fun))
stats::qqline(stats::residuals(lm_q1_fun))

lm_q2_fun <- stats::lm(fun_q2 ~ exploratory + 
                         dominant_tree + 
                         RA_forest + 
                         bac_q2 +
                         alg_q2,
                       data = metadata_alpha_scaled)
base::summary(lm_q2_fun)  
base::plot(allEffects(lm_q2_fun))
base::plot(stats::residuals(lm_q2_fun), stats::fitted(lm_q2_fun))
stats::qqnorm(stats::residuals(lm_q2_fun))
stats::qqline(stats::residuals(lm_q2_fun))

# Correct the p-values using the Benjamini-Hochberg procedure. 

p_step1 <- (summary(lm_q0_alg)$coefficients[-1,4])
p_step1 <- c(p_step1, summary(lm_q0_bac)$coefficients[-1,4])
p_step1 <- c(p_step1, summary(lm_q0_fun)$coefficients[-1,4])

p_step1 <- c(p_step1, (summary(lm_q1_alg)$coefficients[-1,4]))
p_step1 <- c(p_step1, (summary(lm_q1_bac)$coefficients[-1,4]))
p_step1 <- c(p_step1, (summary(lm_q1_fun)$coefficients[-1,4]))

p_step1 <- c(p_step1, (summary(lm_q2_alg)$coefficients[-1,4]))
p_step1 <- c(p_step1, (summary(lm_q2_bac)$coefficients[-1,4]))
p_step1 <- c(p_step1, (summary(lm_q2_fun)$coefficients[-1,4]))

p_adj_step1 <- p.adjust(p_step1, method = "fdr")

p_vals_step1 <- data.frame(names(p_step1) ,round(p_step1, 3), round(p_adj_step1,3))

############
# Second step 
###########

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
                         fun_q0,
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
                         fun_q1,
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
                         fun_q2,
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
                         alg_q0,
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
                         alg_q1,
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
                         alg_q2,
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
                         alg_q0,
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
                         alg_q1,
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
                         alg_q2,
                       data = metadata_alpha_scaled,
                       na.action = "na.fail")
base::summary(lm_q2_fun_second)  
base::plot(allEffects(lm_q2_fun_second))
base::plot(stats::residuals(lm_q2_fun_second), stats::fitted(lm_q2_fun_second))
stats::qqnorm(stats::residuals(lm_q2_fun_second))
stats::qqline(stats::residuals(lm_q2_fun_second))

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

p_vals_step2 <- data.frame(names(p_step2) ,round(p_step2, 3), round(p_adj_step2,3))

##---------------------------------------------------------------
##                   Variance Partitioning                      -
##---------------------------------------------------------------

# A = biotic; B = abiotic, C = exploratory, AB = biotic + abiotic
# AC = biotic + exploratory, BC = abiotic + exploratory 

###
# Algae
###
#q0 
lm_q0_alg_second_biotic <- stats::lm(alg_q0 ~ bac_q0 +
                                       fun_q0,
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q0_alg_second_abiotic <- stats::lm(alg_q0 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest,
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q0_alg_second_geo <- stats::lm(alg_q0 ~ exploratory,
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q0_alg_second_geo_bio <- stats::lm(alg_q0 ~ exploratory +
                                        bac_q0 +
                                        fun_q0,
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
                                         RA_forest,
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
                                         RA_forest,
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
                                       fun_q1,
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q1_alg_second_abiotic <- stats::lm(alg_q1 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest,
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q1_alg_second_geo <- stats::lm(alg_q1 ~ exploratory,
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q1_alg_second_geo_bio <- stats::lm(alg_q1 ~ exploratory +
                                        bac_q1 +
                                        fun_q1,
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
                                         RA_forest,
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
                                         RA_forest,
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
                                       fun_q2,
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q2_alg_second_abiotic <- stats::lm(alg_q2 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest,
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q2_alg_second_geo <- stats::lm(alg_q2 ~ exploratory,
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q2_alg_second_geo_bio <- stats::lm(alg_q2 ~ exploratory +
                                        bac_q2 +
                                        fun_q2,
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
                                         RA_forest,
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
                                         RA_forest,
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
                                       alg_q0,
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q0_fun_second_abiotic <- stats::lm(fun_q0 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest,
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q0_fun_second_geo <- stats::lm(fun_q0 ~ exploratory,
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q0_fun_second_geo_bio <- stats::lm(fun_q0 ~ exploratory +
                                        bac_q0 +
                                        alg_q0,
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
                                         RA_forest,
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
                                         RA_forest,
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
                                       alg_q1,
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q1_fun_second_abiotic <- stats::lm(fun_q1 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest,
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q1_fun_second_geo <- stats::lm(fun_q1 ~ exploratory,
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q1_fun_second_geo_bio <- stats::lm(fun_q1 ~ exploratory +
                                        bac_q1 +
                                        alg_q1,
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
                                         RA_forest,
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
                                         RA_forest,
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
                                       alg_q2,
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q2_fun_second_abiotic <- stats::lm(fun_q2 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest,
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q2_fun_second_geo <- stats::lm(fun_q2 ~ exploratory,
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q2_fun_second_geo_bio <- stats::lm(fun_q2 ~ exploratory +
                                        bac_q2 +
                                        alg_q2,
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
                                         RA_forest,
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
                                         RA_forest,
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
                                       alg_q0,
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q0_bac_second_abiotic <- stats::lm(bac_q0 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest,
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q0_bac_second_geo <- stats::lm(bac_q0 ~ exploratory,
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q0_bac_second_geo_bio <- stats::lm(bac_q0 ~ exploratory +
                                        fun_q0 +
                                        alg_q0,
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
                                         RA_forest,
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
                                         RA_forest,
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
                                       alg_q1,
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q1_bac_second_abiotic <- stats::lm(bac_q1 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest,
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q1_bac_second_geo <- stats::lm(bac_q1 ~ exploratory,
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q1_bac_second_geo_bio <- stats::lm(bac_q1 ~ exploratory +
                                        fun_q1 +
                                        alg_q1,
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
                                         RA_forest,
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
                                         RA_forest,
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
                                       alg_q2,
                                     data = metadata_alpha_scaled,
                                     na.action = "na.fail")

lm_q2_bac_second_abiotic <- stats::lm(bac_q2 ~ rH_200 + 
                                        Ta_200 +
                                        stand_density_basal_area +
                                        DBH_avg +
                                        d_gini +
                                        canopy_openness_2019 +
                                        dom_tot_ratio + 
                                        RA_forest,
                                      data = metadata_alpha_scaled,
                                      na.action = "na.fail")

lm_q2_bac_second_geo <- stats::lm(bac_q2 ~ exploratory,
                                  data = metadata_alpha_scaled,
                                  na.action = "na.fail")

lm_q2_bac_second_geo_bio <- stats::lm(bac_q2 ~ exploratory +
                                        fun_q2 +
                                        alg_q2,
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
                                         RA_forest,
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
                                         RA_forest,
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

saveRDS(alg_variance_lm, here::here("alg_variance_lm.rds"))
saveRDS(fun_variance_lm, here::here("fun_variance_lm.rds"))
saveRDS(bac_variance_lm, here::here("bac_variance_lm.rds"))

#####
# Effects plots
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
  geom_text(data = alg_q0_pred_fun, aes(x = max(x) + 0.1,
                                        y = max(predicted)), label = "*",
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
  geom_text(data = alg_q2_pred_bac, aes(x = max(x) + 0.1,
                                        y = max(predicted)), label = "*",
            size = 4) + 
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
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
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
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
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
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
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
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
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
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Percentage of dominant tree") +
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
        text = element_text(size = 7, family = "sans"), axis.title.y = element_blank()) +
  labs(    
       x= "Surrounding Forest Area") +
  scale_y_continuous(labels = function(x) format(x, nsmall = 1))
algae_forest_effect

###########
# Effects on bacteria
###########

bac_q0_pred_fun <- ggeffects::ggpredict(lm_q0_bac_second, terms = "fun_q0")
bac_q1_pred_fun <- ggeffects::ggpredict(lm_q1_bac_second, terms = "fun_q1")
bac_q2_pred_fun <- ggeffects::ggpredict(lm_q2_bac_second, terms = "fun_q2")

bacteria_fungal_effect <- ggplot() +
  geom_line(data = bac_q0_pred_fun,mapping = aes(x = x, y = predicted),
            color = bacteria_col, linewidth = 0.7) + 
  geom_text(data = bac_q0_pred_fun, aes(x = max(x) + 0.1,
                                        y = max(predicted)), label = "*",
            size = 4) +
  geom_line(data = bac_q1_pred_fun,mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotdash", linewidth = 0.7) +
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
  geom_text(data = bac_q2_pred_alg, aes(x = max(x) + 0.1,
                                        y = max(predicted)), label = "*",
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
            color = bacteria_col, linetype = "dotted", linewidth = 0.7) + 
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
  geom_text(data = bac_q0_pred_d_gini, aes(x = max(x) + 0.1,
                                        y = max(predicted)), label = "*",
            size = 4)+
  geom_line(data = bac_q1_pred_d_gini, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotdash", linewidth = 0.7) + 
  geom_text(data = bac_q1_pred_d_gini, aes(x = max(x) + 0.1,
                                        y = max(predicted)), label = "*",
            size = 4)+
  geom_line(data = bac_q2_pred_d_gini, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotted", linewidth = 0.7) + 
  geom_text(data = bac_q2_pred_d_gini, aes(x = max(x) + 0.1,
                                        y = max(predicted) - 0.05), label = "*",
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
  geom_text(data = bac_q0_pred_canopy_openness_2019, aes(x = max(x) + 0.1,
                                           y = min(predicted)), label = "*",
            size = 4) +
  geom_line(data = bac_q1_pred_canopy_openness_2019, mapping = aes(x = x, y = predicted),
            color = bacteria_col, linetype = "dotdash", linewidth = 0.7) + 
  geom_text(data = bac_q1_pred_canopy_openness_2019, aes(x = max(x) + 0.1,
                                           y = min(predicted)), label = "*",
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
       x= "Percentage of dominant tree") +
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

###########
# Effects on Fungi
###########

fun_q0_pred_bac <- ggeffects::ggpredict(lm_q0_fun_second, terms = "bac_q0")
fun_q1_pred_bac <- ggeffects::ggpredict(lm_q1_fun_second, terms = "bac_q1")
fun_q2_pred_bac <- ggeffects::ggpredict(lm_q2_fun_second, terms = "bac_q2")

fungi_bacterial_effect <- ggplot() +
  geom_line(data = fun_q0_pred_bac,mapping = aes(x = x, y = predicted),
            color = fungi_col, linewidth = 0.7) + 
  geom_text(data = fun_q0_pred_bac, aes(x = max(x) + 0.1,
                                        y = max(predicted)), label = "*",
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
  geom_text(data = fun_q0_pred_alg, aes(x = max(x) + 0.1,
                                        y = max(predicted)), label = "*",
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
  geom_text(data = fun_q0_pred_rH_200, aes(x = max(x) + 0.1,
                                        y = max(predicted)), label = "*",
            size = 4) +
  geom_line(data = fun_q1_pred_rH_200, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotdash", linewidth = 0.7) + 
  geom_text(data = fun_q1_pred_rH_200, aes(x = max(x) + 0.1,
                                        y = max(predicted) - 0.1), label = "*",
            size = 4) +
  geom_line(data = fun_q2_pred_rH_200, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotted", linewidth = 0.7) + 
  geom_text(data = fun_q2_pred_rH_200, aes(x = max(x) + 0.1,
                                        y = max(predicted)), label = "*",
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
  geom_text(data = fun_q0_pred_canopy_openness_2019, aes(x = max(x) + 0.1,
                                                         y = max(predicted) ), label = "*",
            size = 4) +
  geom_line(data = fun_q1_pred_canopy_openness_2019, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotdash", linewidth = 0.7) + 
  geom_text(data = fun_q1_pred_canopy_openness_2019, aes(x = max(x) + 0.1,
                                                         y = max(predicted) - 0.1), label = "*",
            size = 4) +
  geom_line(data = fun_q2_pred_canopy_openness_2019, mapping = aes(x = x, y = predicted),
            color = fungi_col, linetype = "dotted", linewidth = 0.7) + 
  geom_text(data = fun_q2_pred_canopy_openness_2019, aes(x = max(x) + 0.1,
                                                         y = max(predicted)), label = "*",
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
       x= "Percentage of dominant tree") +
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

############
# Arrange in final figure
############

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

effects_algae1 <- ggpubr::ggarrange(algae_bacterial_effect, algae_fungal_effect,
                                    nrow = 1, ncol = 2) %>% 
  ggpubr::annotate_figure(top = ggpubr::text_grob("Biotic effects", size = 7, family = "sans", face = "bold"))

effects_algae2 <- ggpubr::ggarrange(algae_canopy_effect, algae_humidity_effect, 
                                   nrow = 1, ncol = 2) %>% 
  ggpubr::annotate_figure(top = ggpubr::text_grob("Abiotic effects", size = 7, family = "sans", face = "bold"))

effects_algae3 <- ggpubr::ggarrange(effects_algae1,
                                    effects_algae2,
                                    nrow = 1, ncol = 2) 

effects_fungi <- ggpubr::ggarrange(fungi_bacterial_effect, fungi_algal_effect, fungi_canopy_effect, fungi_humidity_effect, 
                                   nrow = 1, ncol = 4) %>% 
  ggpubr::annotate_figure(top = ggpubr::text_grob("",  color = fungi_col, face = "bold", size = 7)) 

effects_bacteria <- ggpubr::ggarrange(bacteria_fungal_effect, bacteria_algal_effect, bacteria_canopy_effect, bacteria_humidity_effect,
                                      nrow = 1, ncol = 4) %>% 
  ggpubr::annotate_figure(top = ggpubr::text_grob("",  color = bacteria_col, face = "bold", size = 7))

effects_alpha <- ggpubr::ggarrange(effects_algae3, effects_fungi, effects_bacteria,
                                   nrow = 3, common.legend = T, legend = "bottom",
                                   legend.grob = alpha_legend)
effects_alpha
ggplot2::ggsave("effects_alpha.pdf", plot = effects_alpha, device = cairo_pdf,
       width = 175, height = 150, units = "mm")

ggplot2::ggsave("effects_alpha.png", plot = effects_alpha, device = png,
       width = 175, height = 150, units = "mm")


