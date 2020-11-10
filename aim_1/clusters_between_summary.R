##%######################################################%##
#                                                          #
####       Find Clusters between 1.5 and 4.5% GD        ####
#                       10/11/2020                         #
##%######################################################%##

# So far I have generated .csv files for clusters with less than 1.5% GD and
# less than 4.5% GD, but no file for those between 1.5 and 4.5% GD. While this
# is fairly straightforward to do by hand, for reproducibility I will do this in R. 

setwd("~/Documents/PhD/Data/aim_1/SEARCH_clustering")

# Libraries ----

library(readr)
library(dplyr)

# Import data ----

pol_1.5 <- read_csv("POL_1.5_clusters.csv")
pol_4.5 <- read_csv("POL_4.5_clusters.csv")

gag_1.5 <- read_csv("GAG_1.5_clusters.csv")
gag_4.5 <- read_csv("GAG_4.5_clusters.csv")

gag_trim_1.5 <- read_csv("GAG_trim_1.5_clusters.csv")
gag_trim_4.5 <- read_csv("GAG_trim_4.5_clusters.csv")

# Function to find the clusters that are between 1.5 and 4.5% GD ----
# This function includes nested clusters in the between group, so nested clusters
# are reported twice, once as a smaller 1.5% cluster and once as a larger 4.5% cluster.

find.between.clusters <- function(x, y){
  
  # Make sure the IDs are in ascending order:
  x <- x[order(x$`4.5_ID`),]
  
  # Get the cluster IDs for clusters that are less than 1.5% GD:
  numbers_1.5 <- x$`4.5_ID`
  
  # Get the clusters that are not nested and greater than 1.5% GD:
  between_1.5_and_4.5 <- dplyr::filter(y,
                                       !(X1 %in% numbers_1.5))
  
  # Get the clusters that are less than 1.5% GD:
  less_than_1.5 <- dplyr::filter(y,
                                 X1 %in% numbers_1.5)
  
  # Find out if the clusters that are less than 1.5% GD are nested within one greater than 1.5% GD:
  nested <- ifelse(less_than_1.5$X1 == x$`4.5_ID` & less_than_1.5$cluster_size != x$cluster_size, "nested", "not-nested")
  
  # Append that information:
  less_than_1.5 <- cbind(less_than_1.5, nested)
  
  # Select the nested clusters:
  less_than_1.5 <- dplyr::filter(less_than_1.5, nested == "nested")
  less_than_1.5 <- dplyr::select(less_than_1.5, -nested)
  
  
  # Make a new dataframe for clusters between 1.5 and 4.5% GD, including nested clusters:
  final <- dplyr::union(between_1.5_and_4.5, less_than_1.5)
  
  # Rename X1:
  final <- dplyr::rename(final, "4.5_ID" = "X1")
  
  final # Remember to add this or the function won't return anything!
}

# Run the function ----
pol_between_clusters <- find.between.clusters(x = pol_1.5, y = pol_4.5)
gag_between_clusters <- find.between.clusters(x = gag_1.5, y = gag_4.5)
gag_trim_between_clusters <- find.between.clusters(x = gag_trim_1.5, y = gag_trim_4.5)

# Write out ----
write_csv(pol_between_clusters, "POL_between_clusters.csv")
write_csv(gag_between_clusters, "GAG_between_clusters.csv")
write_csv(gag_trim_between_clusters, "GAG_trim_between_clusters.csv")
