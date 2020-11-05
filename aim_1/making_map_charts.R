##%######################################################%##
#                                                          #
####               Making Pie Chart Maps                ####
#                       05/11/2020                         #
##%######################################################%##

# This script is to make figures where pie charts for subtype distribution
# are overlayed on a map. In this case, it is done using only POL subtypes,
# for making a figure for my EID Symposium poster.

setwd("~/Documents/PhD/Data/aim_1/SEARCH_subtyping/A2_distinction/map_charts")

# Libraries ----

library(maps)
library(mapdata)
library(ggmap)
library(ggplot2)
library(scatterpie)

# Aesthetic pre-requisites ----

# Create a blank theme to simplify the plots (if needed)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

# Add in the colour palette for the subtypes ----

subtype_colours5 <- c("#A3A380", "#D6CE93", "#EFEBCE", "#D8A48F", "#BB8588", "#CD8162")
names(subtype_colours5) <- c("A1", "A2", "D", "C", "G", "Recombinants")

# Make background maps ----

# Get the outlines of Uganda and Kenya
uganda_kenya <- map_data(map = "world", regions = c('kenya', 'uganda'))

# Plot the Uganda and Kenya map
map <- ggplot() + 
  geom_polygon(data = uganda_kenya, aes(x=long, y = lat, group = group), 
               fill = "#FFC65CB1", color = "black") + 
  coord_fixed(1)

# Plot the world map
world <- ggplot() + 
  geom_polygon(data = map_data(map = "world"), aes(x=long, y = lat, group = group), 
               fill = "#82828247", color = "black") + 
  coord_fixed(1.3) +
  blank_theme # this custom theme makes it minimalistic

# Highlight Uganda and Kenya on the world map
world_highlighted <- world + geom_polygon(data = uganda_kenya, aes(x=long, y = lat, group = group), 
                     fill = "#F0B00E", color = "black") + 
  coord_fixed(1.3) +
  blank_theme

# Export world map (we will later add it in PowerPoint to pinpoint the location of Uganda and Kenya)
ggsave("world_highlighted.pdf", plot = world_highlighted, width = 13, height = 10, units = c("in"))

# Assemble the data for the pie charts ----

# In this case, I am making three pie charts, for Western Uganda, Eastern Uganda 
# and Kenya, showing the distribution of HIV-1 subtypes accross the regions 
# (but only from POL sequences, n = 475).

# As it isn't that much data, I decided to do this manually. Each row of the data
# frame represents Western Uganda, Eastern Uganda and Kenya, respectively. 

long <- c(31, 34, 35)               # the longitude of where I want to place each pie chart.
lat <- c(-0.25, 1, -0.5)            # the latitude of where I want to place each chart.
region <- c("1", "2", "3")          # not entirely sure what this does
A1 <- c(48.6, 62.5, 65.1)           # the percentage of A1 subtype in each region (from previous calculations in "subtype_summaries.R")
A2 <- c(0, 0, 1.74)                 # the percentage of A2
D <- c(24.7, 25, 11.6)              # the percentage of D
C <- c(8.1, 0, 4.65)                # the percentage of C
G <- c(0.405, 0, 0)                 # the percentage of G
Recombinants <- c(18.2, 12.5, 16.9) # the percentage of recombinants
radius <- c(247, 56, 172)           # the size of the pie chart, proportional to number of sequences

# Assemble the data frame
scatterpie <- data.frame(long, lat, region, A1, A2, D, C, G, Recombinants, radius)

# Overlay the pie charts on the map ----

# There are different designs that I have come up with. One with the pie chart
# legend on Kenya, the other with the legend outside, at the top of the plot. 
# For now the latter is my favourite, but keeping both in here for reference. 

# Remember that if I want to remove the latitude and longitude I can do so like this:

      # blank_theme +
      # theme(axis.text.x=element_blank()) +
      # theme(axis.text.y=element_blank())

# Plot with legend in Kenya ----

legend_in <- map + geom_scatterpie(aes(x=long, y=lat, group=region, r = radius/250), data=scatterpie,
                      cols= c("A1", "A2", "D", "C", "G", "Recombinants"), pie_scale = 8.5, color = "black") +
  scale_fill_manual(values = subtype_colours5, 
                    breaks = c("A1", "A2", "D", "C", "G", "Recombinants"),
                    labels = c("A1", "A2", "D", "C", "G", "Recombinants"),
                    name = "HIV-1 subtypes") +
  theme(legend.position = c(0.75, 0.55),
        legend.background = element_rect(fill= "#FFFFFF7C", size=0.5, linetype="solid"),
        legend.key.size = unit(1.25, "cm"),
        legend.key.width = unit(1.25,"cm"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14)) +
  labs(x ="\nLongitude", y = "Latitude\n") +
  theme(axis.title = element_text(size = 15, face = "italic"),
        axis.text = element_text(size = 13)) +
  annotate("label", x = 31, y = 1, label = "Western Uganda", fill = "white") +
  annotate("label", x = 34, y = 1.5, label = "Eastern Uganda", fill = "white") +
  annotate("label", x = 35, y = 0.45, label = "Kenya", fill = "white")

ggsave("pol_distribution_legend_in.pdf", plot = legend_in, width = 13, height = 10, units = c("in"))


# Plot with legend at the top (favourite!) ----

legend_out <- map + geom_scatterpie(aes(x=long, y=lat, group=region, r = radius/250), data=scatterpie,
                      cols= c("A1", "A2", "D", "C", "G", "Recombinants"), pie_scale = 8.5, color = "black") +
  scale_fill_manual(values = subtype_colours5, 
                    breaks = c("A1", "A2", "D", "C", "G", "Recombinants"),
                    labels = c("A1", "A2", "D", "C", "G", "Recombinants"),
                    name = "HIV-1 subtypes",
                    guide = guide_legend(
                      direction = "horizontal",
                      title.position = "top",
                      label.position = "bottom",
                      nrow = 1,
                      label.hjust = 0.5,
                      label.vjust = 1,
                      label.theme = element_text(angle = 0))) +
  theme(legend.position = c(0.205, 0.935),
        legend.background = element_rect(fill= NA, size=0.5, linetype="solid"),
        legend.key.size = unit(1.25, "cm"),
        legend.key.width = unit(1.25,"cm"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14)) +
  labs(x ="\nLongitude", y = "Latitude\n") +
  theme(axis.title = element_text(size = 15, face = "italic"),
        axis.text = element_text(size = 13)) +
  annotate("label", x = 31, y = 1, label = "Western Uganda", fill = "white", size = 4.5) +
  annotate("label", x = 34, y = 1.5, label = "Eastern Uganda", fill = "white", size = 4.5) +
  annotate("label", x = 35, y = 0.45, label = "Kenya", fill = "white", size = 4.5)

ggsave("pol_distribution_legend_out.pdf", plot = legend_out, width = 13, height = 10, units = c("in"))


# Finishing the plots ----

# For now, I am finishing the plots in PowerPoint, where I am overlaying the world
# map and the distribution pie charts to show the location of the countries.