##%######################################################%##
#                                                          #
####            Joining GAG and POL datasets            ####
#                       27/11/2020                         #
##%######################################################%##

setwd("~/Documents/PhD/Data/aim_1/SEARCH_sliding_window")

# The role of this script is to join both the gag and pol datasets into one large
# dataset (with a lot of missing data) to carry out cluster analysis using a sliding
# window technique using all the sequences together.

# The concat() function used to concatenate the sequences together was modified from
# a script by T.F. Khang (https://figshare.com/articles/concatenation_R/3839538)

# Libraries ----

library(seqinr)     # To concatenate and edit the sequences
library(stringr)    # To extract the sequence names
library(phylotools) # To rename the sequences (rename.fasta)

# Import data ----

GAG_sequences <- seqinr::read.fasta("SEARCH_gag_549_3primetrim.fasta")
POL_sequences <- seqinr::read.fasta("SEARCH_pol_488.fasta")

# Generate new sequence names so they match between gag and pol ----

# GAG
names_GAG <- names(GAG_sequences)
for(s in 1:length(names_GAG)){
  names_GAG_new <- str_extract(names_GAG, "\\d{5,}")
}

names_GAG_new
length(unique(names_GAG_new))

names_GAG <- as.data.frame(names_GAG)
names_GAG$names_GAG_new <- names_GAG_new

# POL
names_POL <- names(POL_sequences)
for(s in 1:length(names_POL)){
  names_POL_new <- str_extract(names_POL, "\\d{5,}")
}

names_POL_new
length(unique(names_POL_new))

names_POL <- as.data.frame(names_POL)
names_POL$names_POL_new <- names_POL_new

# Make manual changes so that sequences match EpiData (for when needed) ----

# GAG Changes
names_GAG$names_GAG_new[names_GAG$names_GAG_new == "31393849"] <- "3139384"
names_GAG$names_GAG_new[names_GAG$names_GAG_new == "31357999"] <- "3135799"
names_GAG$names_GAG_new[names_GAG$names_GAG_new == "31350159"] <- "3135015"
names_GAG$names_GAG_new[names_GAG$names_GAG_new == "31393549"] <- "3139354"
names_GAG$names_GAG_new[names_GAG$names_GAG_new == "31384309"] <- "3138430"
names_GAG$names_GAG_new[names_GAG$names_GAG_new == "31402659"] <- "3140265"
names_GAG$names_GAG_new[names_GAG$names_GAG_new == "31358089"] <- "3135808"

# POL Changes
names_POL$names_POL_new[names_POL$names_POL_new == "31393849"] <- "3139384"
names_POL$names_POL_new[names_POL$names_POL_new == "31357999"] <- "3135799"
names_POL$names_POL_new[names_POL$names_POL_new == "31350159"] <- "3135015"
names_POL$names_POL_new[names_POL$names_POL_new == "31393549"] <- "3139354"
names_POL$names_POL_new[names_POL$names_POL_new == "31384309"] <- "3138430"
names_POL$names_POL_new[names_POL$names_POL_new == "31402659"] <- "3140265"

# Rename the sequence datasets ----

rename.fasta(infile = "SEARCH_gag_549_3primetrim.fasta", ref_table = names_GAG, outfile = "renamed_GAG.fasta")
GAG_renamed <- seqinr::read.fasta("renamed_GAG.fasta")

rename.fasta(infile = "SEARCH_pol_488.fasta", ref_table = names_POL, outfile = "renamed_POL.fasta")
POL_renamed <- seqinr::read.fasta("renamed_POL.fasta")

# Join GAG and POL datasets ----

sum(names(POL_renamed) %in% names(GAG_renamed)) 
# 292 sequence names match, so the joint dataset will be 745 sequences

# Functions to concatenate multiple aligned sequences in .fasta files, modified from
# T.F. Khang (2016). They will concatenate seq2 (POL) after seq1 (GAG) in a new .fasta file.

# concat() will not leave any gaps between GAG and POL sequences.
# concat.gap() will add a 50bp gap between GAG and POL for legibility. 

concat <- function(seq1, seq2){
  
  both <- sort(union (names(seq1), names(seq2)))
  seqjoin <- vector("list", length(both))
  names(seqjoin) <- both
  
  for(i in 1:length(both)){
    # Find the index of the matching sequence in the other .fasta file
    idmatch_seq1 <- which(names(seq1) == both[i])
    idmatch_seq2 <- which(names(seq2) == both[i])
    
    if(length(idmatch_seq1) > 0 & length(idmatch_seq2) > 0){
      # Join the two aligned sequences of the same ID
      seqjoin[[i]] <- c(seq1[[idmatch_seq1]], seq2[[idmatch_seq2]])
      
    }
    else if(length(idmatch_seq1) > 0 & length(idmatch_seq2) == 0){
      seqjoin[[i]] <- c(seq1[[idmatch_seq1]], rep("-",length(seq2[[1]])))
      
    }
    
    else if(length(idmatch_seq1) == 0 & length(idmatch_seq2) > 0){
      seqjoin[[i]] <- c(rep("-",length(seq1[[1]])), seq2[[idmatch_seq2]])
    }
    
  }
  
  # Coerce into a SeqFastadna object, so that write.fasta can be applied
  seqjoin <- lapply(seqjoin, as.SeqFastadna)
  
  return(seqjoin)
  
}
concat.gap <- function(seq1, seq2){
  
  both <- sort(union (names(seq1), names(seq2)))
  seqjoin <- vector("list", length(both))
  names(seqjoin) <- both
  
  for(i in 1:length(both)){
    # Find the index of the matching sequence in the other .fasta file
    idmatch_seq1 <- which(names(seq1) == both[i])
    idmatch_seq2 <- which(names(seq2) == both[i])
    
    if(length(idmatch_seq1) > 0 & length(idmatch_seq2) > 0){
      # Join the two aligned sequences of the same ID
      seqjoin[[i]] <- c(seq1[[idmatch_seq1]], rep("-", 50), seq2[[idmatch_seq2]])
      
    }
    else if(length(idmatch_seq1) > 0 & length(idmatch_seq2) == 0){
      seqjoin[[i]] <- c(seq1[[idmatch_seq1]], rep("-", 50), rep("-",length(seq2[[1]])))
      
    }
    
    else if(length(idmatch_seq1) == 0 & length(idmatch_seq2) > 0){
      seqjoin[[i]] <- c(rep("-",length(seq1[[1]])), rep("-", 50), seq2[[idmatch_seq2]])
    }
    
  }
  
  # Coerce into a SeqFastadna object, so that write.fasta can be applied
  seqjoin <- lapply(seqjoin, as.SeqFastadna)
  
  return(seqjoin)
  
}

GAG_POL <- concat(seq1 = GAG_renamed, seq2 = POL_renamed)
GAG_POL_gap <- concat.gap(seq1 = GAG_renamed, seq2 = POL_renamed)

write.fasta(GAG_POL, names(GAG_POL), file = "GAG_POL_no_gap.fasta")
write.fasta(GAG_POL_gap, names(GAG_POL_gap), file = "GAG_POL_50bp_gap.fasta")



