##%######################################################%##
#                                                          #
####       Characterising Clusters from Sequences       ####
#                       09/11/2020                         #
##%######################################################%##

# The previous script I have, cluster_sequences_summary.R, was to retrieve the
# sequences and their associated metadata for each of the clusters, but did not 
# provide a summary of what the clusters actually are as a whole. This is what this
# script does, it compares all of the sequences within a cluster and decides if it
# is a mixed or "pure" cluster with regards to each of the characteristics. 

# As of 06/11/2020, this is only for POL 1.5% GD clusters. 

setwd("~/Documents/PhD/Data/aim_1/SEARCH_clustering")

# Libraries ----

library(readr)        # for data import
library(tidyr)        # for data manipulation
library(dplyr)        # for data manipulation

# Import data ----

POL_clusters <- read_csv("POL_cluster_sequences.csv", 
                         col_types = cols(gag_seq = col_skip(), 
                                          pol_seq = col_skip(), 
                                          has_seq = col_skip(), 
                                          search_id = col_skip(), 
                                          collection_date = col_skip(), 
                                          region_name = col_skip(), 
                                          age = col_skip(), 
                                          community_name = col_skip(), 
                                          incident = col_skip(), 
                                          male = col_skip(), 
                                          occupation = col_skip(), 
                                          Gender = col_factor(levels = c("Female", "Male")), 
                                          `HIV-1 infection` = col_factor(levels = c("Chronic", "Incident")), 
                                          Region = col_factor(levels = c("Western Uganda", "Eastern Uganda", "Kenya")), 
                                          `Sample collection date` = col_date(format = "%Y-%m-%d"), 
                                          `Age category` = col_factor(levels = c("15 to 20","21 to 49", "50 or above")), 
                                          Occupation = col_factor(levels = c("Formal sector", "High-risk informal sector","Low-risk informal sector", "Other", "No job or disabled", "21")), 
                                          Group = col_factor(levels = c("Intervention","Control")), 
                                          gag_REGA_Subtype = col_skip(), 
                                          gag_trim_REGA_Subtype = col_skip()))

GAG_clusters <- read_csv("GAG_cluster_sequences.csv", 
                         col_types = cols(gag_seq = col_skip(), 
                                          pol_seq = col_skip(), 
                                          has_seq = col_skip(), 
                                          search_id = col_skip(), 
                                          collection_date = col_skip(), 
                                          region_name = col_skip(), 
                                          age = col_skip(), 
                                          community_name = col_skip(), 
                                          incident = col_skip(), 
                                          male = col_skip(), 
                                          occupation = col_skip(), 
                                          Gender = col_factor(levels = c("Female", "Male")), 
                                          `HIV-1 infection` = col_factor(levels = c("Chronic", "Incident")), 
                                          Region = col_factor(levels = c("Western Uganda", "Eastern Uganda", "Kenya")), 
                                          `Sample collection date` = col_date(format = "%Y-%m-%d"), 
                                          `Age category` = col_factor(levels = c("15 to 20","21 to 49", "50 or above")), 
                                          Occupation = col_factor(levels = c("Formal sector", "High-risk informal sector","Low-risk informal sector", "Other", "No job or disabled", "21")), 
                                          Group = col_factor(levels = c("Intervention","Control")), 
                                          pol_REGA_Subtype = col_skip(), 
                                          gag_trim_REGA_Subtype = col_skip()))

GAG_trim_clusters <- read_csv("GAG_trim_cluster_sequences.csv", 
                              col_types = cols(gag_seq = col_skip(), 
                                               pol_seq = col_skip(), 
                                               has_seq = col_skip(), 
                                               search_id = col_skip(), 
                                               collection_date = col_skip(), 
                                               region_name = col_skip(), 
                                               age = col_skip(), 
                                               community_name = col_skip(), 
                                               incident = col_skip(), 
                                               male = col_skip(), 
                                               occupation = col_skip(), 
                                               Gender = col_factor(levels = c("Female", "Male")), 
                                               `HIV-1 infection` = col_factor(levels = c("Chronic", "Incident")), 
                                               Region = col_factor(levels = c("Western Uganda", "Eastern Uganda", "Kenya")), 
                                               `Sample collection date` = col_date(format = "%Y-%m-%d"), 
                                               `Age category` = col_factor(levels = c("15 to 20","21 to 49", "50 or above")), 
                                               Occupation = col_factor(levels = c("Formal sector", "High-risk informal sector","Low-risk informal sector", "Other", "No job or disabled", "21")), 
                                               Group = col_factor(levels = c("Intervention","Control")), 
                                               pol_REGA_Subtype = col_skip(), 
                                               gag_REGA_Subtype = col_skip()))

# Select only 1.5% GD Clusters ----

# POL
POL_1.5_with_4.5 <- dplyr::filter(POL_clusters,
                         !is.na(POL_1.5_Cluster))
POL_1.5 <- dplyr::select(POL_1.5_with_4.5, -POL_4.5_Cluster)

# GAG
GAG_1.5_with_4.5 <- dplyr::filter(GAG_clusters,
                         !is.na(GAG_1.5_Cluster))
GAG_1.5 <- dplyr::select(GAG_1.5_with_4.5, -GAG_4.5_Cluster)

# GAG_trim
GAG_trim_1.5_with_4.5 <- dplyr::filter(GAG_trim_clusters,
                              !is.na(GAG_trim_1.5_Cluster))
GAG_trim_1.5 <- dplyr::select(GAG_trim_1.5_with_4.5, -GAG_trim_4.5_Cluster)

# Select only 4.5% GD Clusters ----

# POL
POL_4.5 <- dplyr::select(POL_clusters, -POL_1.5_Cluster)

# GAG
GAG_4.5 <- dplyr::select(GAG_clusters, -GAG_1.5_Cluster)

# GAG_trim
GAG_trim_4.5 <- dplyr::select(GAG_trim_clusters, -GAG_trim_1.5_Cluster)


# Function to characterise clusters, for size up to n = 5 ----

characterise.clusters <- function(x){
  
  colnames(x)[12] <- "REGA_subtype" # Rename so the loop works for POL and GAG and GAG_trim ----
  
  list <- split(x, f = x[,13]) # Split sequences into separate data frames for each cluster number ----
  
  cluster <- list() # Create an empty list to store the output ----
  
  for (i in 1:length(list)) # Loop to characterise the clusters ----
  {
    if (nrow(list[[i]]) == 2){
      
      cluster_sex <- ifelse(list[[i]]$Gender[1] == list[[i]]$Gender[2] & 
                              list[[i]]$Gender[1] =="Male", "Same-sex cluster (Male)", 
                            ifelse(list[[i]]$Gender[1] == list[[i]]$Gender[2] & 
                                     list[[i]]$Gender[1] == "Female", "Same-sex cluster (Female)",
                                   "Heterosexual cluster"))
      
      cluster_status <- ifelse(list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[2] & 
                                 list[[i]]$`HIV-1 infection`[1] =="Chronic", "Chronic cluster", 
                               ifelse(list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[2] & 
                                        list[[i]]$`HIV-1 infection`[1] == "Incident", "Incident cluster",
                                      "Mixed cluster"))
      
      cluster_region <- ifelse(list[[i]]$Region[1] == list[[i]]$Region[2], "Intra-region", 
                               "Inter-region")
      
      cluster_community <- ifelse(list[[i]]$Community[1] == list[[i]]$Community[2], "Intra-community", 
                                  "Inter-community")
      
      cluster_age <- ifelse(list[[i]]$`Age category`[1] == list[[i]]$`Age category`[2] & 
                              list[[i]]$`Age category`[1] == "15 to 20", "15 to 20 cluster", 
                            ifelse(list[[i]]$`Age category`[1] == list[[i]]$`Age category`[2] & 
                                     list[[i]]$`Age category`[1] == "21 to 49", "21 to 49 cluster",
                                   ifelse(list[[i]]$`Age category`[1] == list[[i]]$`Age category`[2] & 
                                            list[[i]]$`Age category`[1] == "50 or above", "50 or above cluster",
                                          "Mixed cluster")))
      
      cluster_occupation <- ifelse(list[[i]]$Occupation[1] == list[[i]]$Occupation[2], "Intra-occupation", 
                                   "Inter-occupation")
      
      cluster_group <- ifelse(list[[i]]$Group[1] == list[[i]]$Group[2] & 
                                list[[i]]$Group[1] == "Intervention", "Intervention cluster", 
                              ifelse(list[[i]]$Group[1] == list[[i]]$Group[2] & 
                                       list[[i]]$Group[1] == "Control", "Control cluster",
                                     "Mixed cluster"))
      
      cluster_subtype <- ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                  list[[i]]$REGA_subtype[1] == "A1", "A1 cluster",
                                ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                         list[[i]]$REGA_subtype[1] == "D", "D cluster",
                                       ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                                list[[i]]$REGA_subtype[1] == "C", "C cluster",
                                              ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                                       list[[i]]$REGA_subtype[1] == "A1,C", "A1,C cluster",
                                                     ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                                              list[[i]]$REGA_subtype[1] == "A1,D", "A1,D cluster",
                                                            "Other recombinant cluster")))))
      
      cluster[[i]] <- data.frame(cluster_sex, cluster_status, cluster_region, cluster_community, cluster_age, cluster_occupation, cluster_group, cluster_subtype)
      
    } else {if (nrow(list[[i]]) == 3){
      
      cluster_sex <- ifelse(list[[i]]$Gender[1] == list[[i]]$Gender[2] & 
                              list[[i]]$Gender[1] == list[[i]]$Gender[3] & 
                              list[[i]]$Gender[1] =="Male", "Same-sex cluster (Male)", 
                            ifelse(list[[i]]$Gender[1] == list[[i]]$Gender[2] & 
                                     list[[i]]$Gender[1] == list[[i]]$Gender[3] & 
                                     list[[i]]$Gender[1] == "Female", "Same-sex cluster (Female)",
                                   "Heterosexual cluster"))
      
      cluster_status <- ifelse(list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[2] & 
                                 list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[3] & 
                                 list[[i]]$`HIV-1 infection`[1] =="Chronic", "Chronic cluster", 
                               ifelse(list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[2] & 
                                        list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[3] & 
                                        list[[i]]$`HIV-1 infection`[1] == "Incident", "Incident cluster",
                                      "Mixed cluster"))
      
      cluster_region <- ifelse(list[[i]]$Region[1] == list[[i]]$Region[2] & 
                                 list[[i]]$Region[1] == list[[i]]$Region[3], "Intra-region", 
                               "Inter-region")
      
      cluster_community <- ifelse(list[[i]]$Community[1] == list[[i]]$Community[2]  & 
                                    list[[i]]$Community[1] == list[[i]]$Community[3], "Intra-community", 
                                  "Inter-community")
      
      cluster_age <- ifelse(list[[i]]$`Age category`[1] == list[[i]]$`Age category`[2] & 
                              list[[i]]$`Age category`[1] == list[[i]]$`Age category`[3] & 
                              list[[i]]$`Age category`[1] == "15 to 20", "15 to 20 cluster", 
                            ifelse(list[[i]]$`Age category`[1] == list[[i]]$`Age category`[2] & 
                                     list[[i]]$`Age category`[1] == list[[i]]$`Age category`[3] & 
                                     list[[i]]$`Age category`[1] == "21 to 49", "21 to 49 cluster",
                                   ifelse(list[[i]]$`Age category`[1] == list[[i]]$`Age category`[2] & 
                                            list[[i]]$`Age category`[1] == list[[i]]$`Age category`[3] & 
                                            list[[i]]$`Age category`[1] == "50 or above", "50 or above cluster",
                                          "Mixed cluster")))
      
      cluster_occupation <- ifelse(list[[i]]$Occupation[1] == list[[i]]$Occupation[2] & 
                                     list[[i]]$Occupation[1] == list[[i]]$Occupation[3], "Intra-occupation", 
                                   "Inter-occupation")
      
      cluster_group <- ifelse(list[[i]]$Group[1] == list[[i]]$Group[2] & 
                                list[[i]]$Group[1] == list[[i]]$Group[3] & 
                                list[[i]]$Group[1] == "Intervention", "Intervention cluster", 
                              ifelse(list[[i]]$Group[1] == list[[i]]$Group[2]  & 
                                       list[[i]]$Group[1] == list[[i]]$Group[3] & 
                                       list[[i]]$Group[1] == "Control", "Control cluster",
                                     "Mixed cluster"))
      
      cluster_subtype <- ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                  list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[3] & 
                                  list[[i]]$REGA_subtype[1] == "A1", "A1 cluster",
                                ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                         list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[3] & 
                                         list[[i]]$REGA_subtype[1] == "D", "D cluster",
                                       ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                                list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[3] & 
                                                list[[i]]$REGA_subtype[1] == "C", "C cluster",
                                              ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                                       list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[3] & 
                                                       list[[i]]$REGA_subtype[1] == "A1,C", "A1,C cluster",
                                                     ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                                              list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[3] & 
                                                              list[[i]]$REGA_subtype[1] == "A1,D", "A1,D cluster",
                                                            "Other recombinant cluster")))))
      
      cluster[[i]] <- data.frame(cluster_sex, cluster_status, cluster_region, cluster_community, cluster_age, cluster_occupation, cluster_group, cluster_subtype)
      
    } else {if (nrow(list[[i]]) == 4){
      
      cluster_sex <- ifelse(list[[i]]$Gender[1] == list[[i]]$Gender[2] & 
                              list[[i]]$Gender[1] == list[[i]]$Gender[3] & 
                              list[[i]]$Gender[1] == list[[i]]$Gender[4] & 
                              list[[i]]$Gender[1] =="Male", "Same-sex cluster (Male)", 
                            ifelse(list[[i]]$Gender[1] == list[[i]]$Gender[2] & 
                                     list[[i]]$Gender[1] == list[[i]]$Gender[3] & 
                                     list[[i]]$Gender[1] == list[[i]]$Gender[4] &
                                     list[[i]]$Gender[1] == "Female", "Same-sex cluster (Female)",
                                   "Heterosexual cluster"))
      
      cluster_status <- ifelse(list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[2] & 
                                 list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[3] & 
                                 list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[4] &
                                 list[[i]]$`HIV-1 infection`[1] =="Chronic", "Chronic cluster", 
                               ifelse(list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[2] & 
                                        list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[3] & 
                                        list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[4] &
                                        list[[i]]$`HIV-1 infection`[1] == "Incident", "Incident cluster",
                                      "Mixed cluster"))
      
      cluster_region <- ifelse(list[[i]]$Region[1] == list[[i]]$Region[2] & 
                                 list[[i]]$Region[1] == list[[i]]$Region[3] &
                                 list[[i]]$Region[1] == list[[i]]$Region[4], "Intra-region", 
                               "Inter-region")
      
      cluster_community <- ifelse(list[[i]]$Community[1] == list[[i]]$Community[2]  & 
                                    list[[i]]$Community[1] == list[[i]]$Community[3] &
                                    list[[i]]$Community[1] == list[[i]]$Community[4], "Intra-community", 
                                  "Inter-community")
      
      cluster_age <- ifelse(list[[i]]$`Age category`[1] == list[[i]]$`Age category`[2] & 
                              list[[i]]$`Age category`[1] == list[[i]]$`Age category`[3] & 
                              list[[i]]$`Age category`[1] == list[[i]]$`Age category`[4] &
                              list[[i]]$`Age category`[1] == "15 to 20", "15 to 20 cluster", 
                            ifelse(list[[i]]$`Age category`[1] == list[[i]]$`Age category`[2] & 
                                     list[[i]]$`Age category`[1] == list[[i]]$`Age category`[3] & 
                                     list[[i]]$`Age category`[1] == list[[i]]$`Age category`[4] &
                                     list[[i]]$`Age category`[1] == "21 to 49", "21 to 49 cluster",
                                   ifelse(list[[i]]$`Age category`[1] == list[[i]]$`Age category`[2] & 
                                            list[[i]]$`Age category`[1] == list[[i]]$`Age category`[3] & 
                                            list[[i]]$`Age category`[1] == list[[i]]$`Age category`[4] &
                                            list[[i]]$`Age category`[1] == "50 or above", "50 or above cluster",
                                          "Mixed cluster")))
      
      cluster_occupation <- ifelse(list[[i]]$Occupation[1] == list[[i]]$Occupation[2] & 
                                     list[[i]]$Occupation[1] == list[[i]]$Occupation[3] &
                                     list[[i]]$Occupation[1] == list[[i]]$Occupation[4], "Intra-occupation", 
                                   "Inter-occupation")
      
      cluster_group <- ifelse(list[[i]]$Group[1] == list[[i]]$Group[2] & 
                                list[[i]]$Group[1] == list[[i]]$Group[3] &
                                list[[i]]$Group[1] == list[[i]]$Group[4] &
                                list[[i]]$Group[1] == "Intervention", "Intervention cluster", 
                              ifelse(list[[i]]$Group[1] == list[[i]]$Group[2]  & 
                                       list[[i]]$Group[1] == list[[i]]$Group[3] & 
                                       list[[i]]$Group[1] == list[[i]]$Group[4] &
                                       list[[i]]$Group[1] == "Control", "Control cluster",
                                     "Mixed cluster"))
      
      cluster_subtype <- ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                  list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[3] & 
                                  list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[4] &
                                  list[[i]]$REGA_subtype[1] == "A1", "A1 cluster",
                                ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                         list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[3] & 
                                         list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[4] &
                                         list[[i]]$REGA_subtype[1] == "D", "D cluster",
                                       ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                                list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[3] &
                                                list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[4] &
                                                list[[i]]$REGA_subtype[1] == "C", "C cluster",
                                              ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                                       list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[3] & 
                                                       list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[4] &
                                                       list[[i]]$REGA_subtype[1] == "A1,C", "A1,C cluster",
                                                     ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                                              list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[3] & 
                                                              list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[4] &
                                                              list[[i]]$REGA_subtype[1] == "A1,D", "A1,D cluster",
                                                            "Other recombinant cluster")))))
      
      cluster[[i]] <- data.frame(cluster_sex, cluster_status, cluster_region, cluster_community, cluster_age, cluster_occupation, cluster_group, cluster_subtype)
      
    } else {if (nrow(list[[i]]) == 5){
      
      cluster_sex <- ifelse(list[[i]]$Gender[1] == list[[i]]$Gender[2] & 
                              list[[i]]$Gender[1] == list[[i]]$Gender[3] & 
                              list[[i]]$Gender[1] == list[[i]]$Gender[4] & 
                              list[[i]]$Gender[1] == list[[i]]$Gender[5] &
                              list[[i]]$Gender[1] =="Male", "Same-sex cluster (Male)", 
                            ifelse(list[[i]]$Gender[1] == list[[i]]$Gender[2] & 
                                     list[[i]]$Gender[1] == list[[i]]$Gender[3] & 
                                     list[[i]]$Gender[1] == list[[i]]$Gender[4] &
                                     list[[i]]$Gender[1] == list[[i]]$Gender[5] &
                                     list[[i]]$Gender[1] == "Female", "Same-sex cluster (Female)",
                                   "Heterosexual cluster"))
      
      cluster_status <- ifelse(list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[2] & 
                                 list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[3] & 
                                 list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[4] &
                                 list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[5] &
                                 list[[i]]$`HIV-1 infection`[1] =="Chronic", "Chronic cluster", 
                               ifelse(list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[2] & 
                                        list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[3] & 
                                        list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[4] &
                                        list[[i]]$`HIV-1 infection`[1] == list[[i]]$`HIV-1 infection`[5] &
                                        list[[i]]$`HIV-1 infection`[1] == "Incident", "Incident cluster",
                                      "Mixed cluster"))
      
      cluster_region <- ifelse(list[[i]]$Region[1] == list[[i]]$Region[2] & 
                                 list[[i]]$Region[1] == list[[i]]$Region[3] &
                                 list[[i]]$Region[1] == list[[i]]$Region[4] &
                                 list[[i]]$Region[1] == list[[i]]$Region[5], "Intra-region", 
                               "Inter-region")
      
      cluster_community <- ifelse(list[[i]]$Community[1] == list[[i]]$Community[2]  & 
                                    list[[i]]$Community[1] == list[[i]]$Community[3] &
                                    list[[i]]$Community[1] == list[[i]]$Community[4] &
                                    list[[i]]$Community[1] == list[[i]]$Community[5], "Intra-community", 
                                  "Inter-community")
      
      cluster_age <- ifelse(list[[i]]$`Age category`[1] == list[[i]]$`Age category`[2] & 
                              list[[i]]$`Age category`[1] == list[[i]]$`Age category`[3] & 
                              list[[i]]$`Age category`[1] == list[[i]]$`Age category`[4] &
                              list[[i]]$`Age category`[1] == list[[i]]$`Age category`[5] &
                              list[[i]]$`Age category`[1] == "15 to 20", "15 to 20 cluster", 
                            ifelse(list[[i]]$`Age category`[1] == list[[i]]$`Age category`[2] & 
                                     list[[i]]$`Age category`[1] == list[[i]]$`Age category`[3] & 
                                     list[[i]]$`Age category`[1] == list[[i]]$`Age category`[4] &
                                     list[[i]]$`Age category`[1] == list[[i]]$`Age category`[5] &
                                     list[[i]]$`Age category`[1] == "21 to 49", "21 to 49 cluster",
                                   ifelse(list[[i]]$`Age category`[1] == list[[i]]$`Age category`[2] & 
                                            list[[i]]$`Age category`[1] == list[[i]]$`Age category`[3] & 
                                            list[[i]]$`Age category`[1] == list[[i]]$`Age category`[4] &
                                            list[[i]]$`Age category`[1] == list[[i]]$`Age category`[5] &
                                            list[[i]]$`Age category`[1] == "50 or above", "50 or above cluster",
                                          "Mixed cluster")))
      
      cluster_occupation <- ifelse(list[[i]]$Occupation[1] == list[[i]]$Occupation[2] & 
                                     list[[i]]$Occupation[1] == list[[i]]$Occupation[3] &
                                     list[[i]]$Occupation[1] == list[[i]]$Occupation[4] &
                                     list[[i]]$Occupation[1] == list[[i]]$Occupation[5], "Intra-occupation", 
                                   "Inter-occupation")
      
      cluster_group <- ifelse(list[[i]]$Group[1] == list[[i]]$Group[2] & 
                                list[[i]]$Group[1] == list[[i]]$Group[3] &
                                list[[i]]$Group[1] == list[[i]]$Group[4] &
                                list[[i]]$Group[1] == list[[i]]$Group[5] &
                                list[[i]]$Group[1] == "Intervention", "Intervention cluster", 
                              ifelse(list[[i]]$Group[1] == list[[i]]$Group[2]  & 
                                       list[[i]]$Group[1] == list[[i]]$Group[3] & 
                                       list[[i]]$Group[1] == list[[i]]$Group[4] &
                                       list[[i]]$Group[1] == list[[i]]$Group[5] &
                                       list[[i]]$Group[1] == "Control", "Control cluster",
                                     "Mixed cluster"))
      
      cluster_subtype <- ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                  list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[3] & 
                                  list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[4] &
                                  list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[5] &
                                  list[[i]]$REGA_subtype[1] == "A1", "A1 cluster",
                                ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                         list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[3] & 
                                         list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[4] &
                                         list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[5] &
                                         list[[i]]$REGA_subtype[1] == "D", "D cluster",
                                       ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                                list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[3] &
                                                list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[4] &
                                                list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[5] &
                                                list[[i]]$REGA_subtype[1] == "C", "C cluster",
                                              ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                                       list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[3] & 
                                                       list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[4] &
                                                       list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[5] &
                                                       list[[i]]$REGA_subtype[1] == "A1,C", "A1,C cluster",
                                                     ifelse(list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[2] & 
                                                              list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[3] & 
                                                              list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[4] &
                                                              list[[i]]$REGA_subtype[1] == list[[i]]$REGA_subtype[5] &
                                                              list[[i]]$REGA_subtype[1] == "A1,D", "A1,D cluster",
                                                            "Other recombinant cluster")))))
      
      cluster[[i]] <- data.frame(cluster_sex, cluster_status, cluster_region, cluster_community, cluster_age, cluster_occupation, cluster_group, cluster_subtype)
    }}}
      "Cluster larger than n = 5"
    }
  } 
  print(cluster)
}

#### For 1.5% GD Clusters ####
# Run the functions ---- 
pol_1.5_clusters <- characterise.clusters(POL_1.5)
gag_1.5_clusters <- characterise.clusters(GAG_1.5)
gag_trim_1.5_clusters <- characterise.clusters(GAG_trim_1.5)

# Convert the list into a data frame ----

# POL
pol_1.5_summary <- data.frame(matrix(unlist(pol_1.5_clusters), nrow=length(pol_1.5_clusters), byrow=T))
colnames(pol_1.5_summary) <- c("cluster_sex", "cluster_status", "cluster_region", "cluster_community", "cluster_age", "cluster_occupation", "cluster_group", "cluster_subtype")

# GAG
gag_1.5_summary <- data.frame(matrix(unlist(gag_1.5_clusters), nrow=length(gag_1.5_clusters), byrow=T))
colnames(gag_1.5_summary) <- c("cluster_sex", "cluster_status", "cluster_region", "cluster_community", "cluster_age", "cluster_occupation", "cluster_group", "cluster_subtype")

# GG_trim
gag_trim_1.5_summary <- data.frame(matrix(unlist(gag_trim_1.5_clusters), nrow=length(gag_trim_1.5_clusters), byrow=T))
colnames(gag_trim_1.5_summary) <- c("cluster_sex", "cluster_status", "cluster_region", "cluster_community", "cluster_age", "cluster_occupation", "cluster_group", "cluster_subtype")

# Add the cluster numbers according to the 4.5% GD clusters, for comparison

# POL
numbers_in_pol_4.5 <- unique(POL_1.5_with_4.5$POL_4.5_Cluster)
numbers_in_pol_4.5 <- sort(numbers_in_pol_4.5)

pol_1.5_summary$`4.5_ID` <- numbers_in_pol_4.5
pol_1.5_summary <- pol_1.5_summary[, c(9, 1:8)]

# GAG
numbers_in_gag_4.5 <- unique(GAG_1.5_with_4.5$GAG_4.5_Cluster)
numbers_in_gag_4.5 <- sort(numbers_in_gag_4.5)
numbers_in_gag_4.5 <- replace(numbers_in_gag_4.5, c(15, 16), numbers_in_gag_4.5[c(16, 15)])
# we need to do the above replacement as the numbers aren't in order

gag_1.5_summary$`4.5_ID` <- numbers_in_gag_4.5
gag_1.5_summary <- gag_1.5_summary[, c(9, 1:8)]

# GAG_trim
numbers_in_gag_trim_4.5 <- unique(GAG_trim_1.5_with_4.5$GAG_trim_4.5_Cluster)
numbers_in_gag_trim_4.5 <- sort(numbers_in_gag_trim_4.5)
numbers_in_gag_trim_4.5 <- replace(numbers_in_gag_trim_4.5, c(10, 11), numbers_in_gag_trim_4.5[c(11, 10)])
# we need to do the above replacement as the numbers aren't in order

gag_trim_1.5_summary$`4.5_ID` <- numbers_in_gag_trim_4.5
gag_trim_1.5_summary <- gag_trim_1.5_summary[, c(9, 1:8)]


#### For 4.5% GD Clusters ####
# Run the functions ---- 
pol_4.5_clusters <- characterise.clusters(POL_4.5)
gag_4.5_clusters <- characterise.clusters(GAG_4.5)
gag_trim_4.5_clusters <- characterise.clusters(GAG_trim_4.5)

# Convert the list into a data frame ----

# POL
pol_4.5_summary <- data.frame(matrix(unlist(pol_4.5_clusters), nrow=length(pol_4.5_clusters), byrow=T))
colnames(pol_4.5_summary) <- c("cluster_sex", "cluster_status", "cluster_region", "cluster_community", "cluster_age", "cluster_occupation", "cluster_group", "cluster_subtype")

# GAG
gag_4.5_summary <- data.frame(matrix(unlist(gag_4.5_clusters), nrow=length(gag_4.5_clusters), byrow=T))
colnames(gag_4.5_summary) <- c("cluster_sex", "cluster_status", "cluster_region", "cluster_community", "cluster_age", "cluster_occupation", "cluster_group", "cluster_subtype")

# GG_trim
gag_trim_4.5_summary <- data.frame(matrix(unlist(gag_trim_4.5_clusters), nrow=length(gag_trim_4.5_clusters), byrow=T))
colnames(gag_trim_4.5_summary) <- c("cluster_sex", "cluster_status", "cluster_region", "cluster_community", "cluster_age", "cluster_occupation", "cluster_group", "cluster_subtype")


# Write out ----
write_csv(pol_1.5_summary, "POL_1.5_clusters.csv")
write_csv(gag_1.5_summary, "GAG_1.5_clusters.csv")
write_csv(gag_trim_1.5_summary, "GAG_trim_1.5_clusters.csv")

write_csv(pol_4.5_summary, "POL_4.5_clusters.csv")
write_csv(gag_4.5_summary, "GAG_4.5_clusters.csv")
write_csv(gag_trim_4.5_summary, "GAG_trim_4.5_clusters.csv")
