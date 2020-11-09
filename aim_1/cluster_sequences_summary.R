##%######################################################%##
#                                                          #
####         Summarising Sequences in Clusters          ####
#                       15/10/2020                         #
##%######################################################%##

setwd("~/Documents/PhD/Data/aim_1/SEARCH_clustering")

# Libraries ----

library(dplyr)
library(readr)
library(summarytools)

# Import data ----

# POL:
pol_1.5 <- read.table("Cluster_Picker_POL/1.5%/EPHpol.fasta_clusterPicks_list.txt", header = TRUE)
pol_1.5 <- dplyr::filter(pol_1.5, ClusterNumber != -1)
pol_4.5 <- read.table("Cluster_Picker_POL/4.5%/EPHpol.fasta_clusterPicks_list.txt", header = TRUE)
pol_4.5 <- dplyr::filter(pol_4.5, ClusterNumber != -1)
pol_common_clusters <- dplyr::semi_join(pol_4.5, pol_1.5, by = "SequenceName")

dplyr::count(pol_1.5, ClusterNumber)
dplyr::count(pol_4.5, ClusterNumber)
dplyr::count(pol_common_clusters, ClusterNumber)

pol_all_clusters <- dplyr::left_join(pol_4.5, pol_1.5, by = "SequenceName")
pol_all_clusters <- dplyr::rename(pol_all_clusters,
                                  "4.5_Cluster" = "ClusterNumber.x",
                                  "1.5_Cluster" = "ClusterNumber.y")

# GAG:

gag_1.5 <- read.table("Cluster_Picker_GAG/1.5%/EPHgag_549.fasta_clusterPicks_list.txt", header = TRUE)
gag_1.5 <- dplyr::filter(gag_1.5, ClusterNumber != -1)
gag_4.5 <- read.table("Cluster_Picker_GAG/4.5%/EPHgag_549.fasta_clusterPicks_list.txt", header = TRUE)
gag_4.5 <- dplyr::filter(gag_4.5, ClusterNumber != -1)
gag_common_clusters <- dplyr::semi_join(gag_4.5, gag_1.5, by = "SequenceName")

dplyr::count(gag_1.5, ClusterNumber)
dplyr::count(gag_4.5, ClusterNumber)
dplyr::count(gag_common_clusters, ClusterNumber)

gag_all_clusters <- dplyr::left_join(gag_4.5, gag_1.5, by = "SequenceName")
gag_all_clusters <- dplyr::rename(gag_all_clusters,
                                  "4.5_Cluster" = "ClusterNumber.x",
                                  "1.5_Cluster" = "ClusterNumber.y")

# GAG_trim:

gag_trim_1.5 <- read.table("Cluster_Picker_GAG_trim/1.5%/EPHgag_549_3trim.fasta_clusterPicks_list.txt", header = TRUE)
gag_trim_1.5 <- dplyr::filter(gag_trim_1.5, ClusterNumber != -1)
gag_trim_4.5 <- read.table("Cluster_Picker_GAG_trim/4.5%/EPHgag_549_3trim.fasta_clusterPicks_list.txt", header = TRUE)
gag_trim_4.5 <- dplyr::filter(gag_trim_4.5, ClusterNumber != -1)
gag_trim_common_clusters <- dplyr::semi_join(gag_trim_4.5, gag_trim_1.5, by = "SequenceName")

dplyr::count(gag_trim_1.5, ClusterNumber)
dplyr::count(gag_trim_4.5, ClusterNumber)
dplyr::count(gag_trim_common_clusters, ClusterNumber)

gag_trim_all_clusters <- dplyr::left_join(gag_trim_4.5, gag_trim_1.5, by = "SequenceName")
gag_trim_all_clusters <- dplyr::rename(gag_trim_all_clusters,
                                  "4.5_Cluster" = "ClusterNumber.x",
                                  "1.5_Cluster" = "ClusterNumber.y")



# Matching with Subtypes and EpiData ----

ALL <- read_csv("~/Documents/PhD/Data/aim_1/SEARCH_subtyping/A2_distinction/all_SEARCH.csv", 
                col_types = cols(digit_ID = col_character(), 
                                 gag_seq = col_factor(levels = c("0", "1")), 
                                 pol_seq = col_factor(levels = c("0", "1")), 
                                 has_seq = col_factor(levels = c("1")), 
                                 collection_date = col_date(format = "%Y-%m-%d"), 
                                 region_name = col_factor(levels = c("Western Uganda", "Eastern Uganda", "Kenya")), 
                                 incident = col_factor(levels = c("0", "1")), 
                                 male = col_factor(levels = c("0", "1")), 
                                 Gender = col_factor(levels = c("Female", "Male")),
                                 `HIV-1 infection` = col_factor(levels = c("Chronic", "Incident")), 
                                 Region = col_factor(levels = c("Western Uganda", "Eastern Uganda", "Kenya")), 
                                 `Sample collection date` = col_date(format = "%Y-%m-%d"), 
                                 `Age category` = col_factor(levels = c("15 to 20", "21 to 49", "50 or above")), 
                                 Occupation = col_factor(levels = c("Formal sector", "High-risk informal sector", "Low-risk informal sector", "Other", "No job or disabled", "21")), 
                                 Group = col_factor(levels = c("Intervention", "Control"))))

ALL$community_name <- as.factor(ALL$community_name)
ALL$Community <- as.factor(ALL$Community)
ALL$occupation <- as.factor(ALL$occupation)
ALL$Occupations <- as.factor(ALL$Occupations)
ALL$gag_REGA_Subtype <- as.factor(ALL$gag_REGA_Subtype)
ALL$gag_trim_REGA_Subtype <- as.factor(ALL$gag_trim_REGA_Subtype)
ALL$pol_REGA_Subtype <- as.factor(ALL$pol_REGA_Subtype)


# POL:

pol_all_clusters <- dplyr::rename(pol_all_clusters, "POL_Name" = "SequenceName")
pol_metadata_clusters <- dplyr::left_join(ALL, pol_all_clusters, by = "POL_Name")
pol_metadata_clusters <- dplyr::filter(pol_metadata_clusters, !is.na(`4.5_Cluster`))
pol_metadata_clusters <- pol_metadata_clusters[,c(1:22,26,30,34:36)]


# GAG:

gag_all_clusters <- dplyr::rename(gag_all_clusters, "GAG_Name" = "SequenceName")
gag_metadata_clusters <- dplyr::left_join(ALL, gag_all_clusters, by = "GAG_Name")
gag_metadata_clusters <- dplyr::filter(gag_metadata_clusters, !is.na(`4.5_Cluster`))
gag_metadata_clusters <- gag_metadata_clusters[,c(1:22,26,30,34, 37, 38)]


# GAG_trim:

gag_trim_all_clusters <- dplyr::rename(gag_trim_all_clusters, "GAG_trim_Name" = "SequenceName")
gag_trim_metadata_clusters <- dplyr::left_join(ALL, gag_trim_all_clusters, by = "GAG_trim_Name")
gag_trim_metadata_clusters <- dplyr::filter(gag_trim_metadata_clusters, !is.na(`4.5_Cluster`))
gag_trim_metadata_clusters <- gag_trim_metadata_clusters[,c(1:22,26,30,34, 39, 40)]



# Writing out ----
write_csv(pol_metadata_clusters, "POL_cluster_sequences.csv")
write_csv(gag_metadata_clusters, "GAG_cluster_sequences.csv")
write_csv(gag_trim_metadata_clusters, "GAG_trim_cluster_sequences.csv")

# Tidying up ----

pol_metadata_clusters <- pol_metadata_clusters[, c(1, 13:22, 25:27)]
gag_metadata_clusters <- gag_metadata_clusters[, c(1, 13:23, 26:27)]
gag_trim_metadata_clusters <- gag_trim_metadata_clusters[,c(1, 13:22, 24, 26:27)]

# Demographically characterising the participants in 4.5% clusters (in general) ----

# POL:
summary_pol_4.5_clusters <- dplyr::select(pol_metadata_clusters, -digit_ID, -`POL_1.5_Cluster`, -`POL_4.5_Cluster`) %>% 
  dfSummary(max.distinct.values = 32)
view(summary_pol_4.5_clusters, file = "summary_pol_4.5_clusters.html")

# GAG:
summary_gag_4.5_clusters <- dplyr::select(gag_metadata_clusters, -digit_ID, -`GAG_1.5_Cluster`, -`GAG_4.5_Cluster`) %>% 
  dfSummary(max.distinct.values = 32)
view(summary_gag_4.5_clusters, file = "summary_gag_4.5_clusters.html")

# GAG_trim:
summary_gag_trim_4.5_clusters <- dplyr::select(gag_trim_metadata_clusters, -digit_ID, -`GAG_trim_1.5_Cluster`, -`GAG_trim_4.5_Cluster`) %>% 
  dfSummary(max.distinct.values = 32)
view(summary_gag_trim_4.5_clusters, file = "summary_gag_trim_4.5_clusters.html")

# Demographically characterising the participants in 1.5% clusters (in general) ----

# POL:
summary_pol_1.5_clusters <- pol_metadata_clusters %>% 
  dplyr::filter(!is.na(`POL_1.5_Cluster`)) %>%
  dplyr::select(-digit_ID, -`POL_1.5_Cluster`, -`POL_4.5_Cluster`) %>% 
  dfSummary(max.distinct.values = 32)
view(summary_pol_1.5_clusters, file = "summary_pol_1.5_clusters.html")

# GAG:
summary_gag_1.5_clusters <- gag_metadata_clusters %>% 
  dplyr::filter(!is.na(`GAG_1.5_Cluster`)) %>%
  dplyr::select(-digit_ID, -`GAG_1.5_Cluster`, -`GAG_4.5_Cluster`) %>% 
  dfSummary(max.distinct.values = 32)
view(summary_gag_1.5_clusters, file = "summary_gag_1.5_clusters.html")

# GAG_trim:
summary_gag_trim_1.5_clusters <- gag_trim_metadata_clusters %>% 
  dplyr::filter(!is.na(`GAG_trim_1.5_Cluster`)) %>%
  dplyr::select(-digit_ID, -`GAG_trim_1.5_Cluster`, -`GAG_trim_4.5_Cluster`) %>% 
  dfSummary(max.distinct.values = 32)
view(summary_gag_trim_1.5_clusters, file = "summary_gag_trim_1.5_clusters.html")
