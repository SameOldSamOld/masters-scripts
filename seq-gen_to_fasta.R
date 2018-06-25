# Sam Old 1/09/17
# Read PHYLIP files from seq-gen
# Seperate data based on species name

# R version 3.4.2 (2017-09-28)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows >= 8 x64 (build 9200)
# Biostrings_2.44.2, phylotools_0.1.2

## Multiple PHYLIP files to be concatenated into one FASTA file

###################################
# Files must have certain format  #
# - filename: runXXX              #
# - species name: G01_XX          #
# - sequence: separated by spaces #
###################################

## Load libraries
library(Biostrings)
library(phylotools)

## Set way directory

setwd("C:/Users/samol/Desktop/4/")
files = list.files("./")
files = sort(files)

## This method uses much less memory but runs a lot slower 
# Benchmark: 1 minute

for( i in files ){
  
  # Read each PHLYIP file
  currentPhy = read.phy(i)
  for( j in 2:length(currentPhy) ){
    
    # Skip header line. 
    # Read info about each protein sequence
    ID = strsplit(currentPhy[j]," +")
    species = strsplit(ID[[1]][1], "_")[[1]][2]; sequence = ID[[1]][2];
    uniqueID = paste0(i, "_", ID[[1]][1])
    
    # Format data as FASTA data
    AA.formatted = BStringSet(sequence)
    names(AA.formatted) = uniqueID
    
    # Append/create data file
    if(file.exists(paste0(species,".fa"))){
      writeXStringSet(AA.formatted, filepath=paste0("./", species,".fa"), format="fasta", append=T)
    } else {
      writeXStringSet(AA.formatted, filepath=paste0("./", species,".fa"), format="fasta", append=F)
    }
  }
}

# THIS WRITES ENTIRE LIST AT ONCE - Very fast but big memory cost 
# Benchmark: 1 second
# #### IMPORTANT #### I broke it somehow. Lists not working as intended any more.
# sequenceList = list()
# ## Read PHYLIP files into a List
# for( i in files ){  
#   # Read in each files 1-by-1
#   currentPhy = read.phy(i)
#   # Nested for loop starts with a 2 to skip header line from phylip file
#   for( j in 2:length(currentPhy) ){
#     # Extract species, sequence, and unique ID for each gene in phylogeny tree
#     ID = strsplit(currentPhy[j]," +")
#     species = strsplit(ID[[1]][1], "_")[[1]][2]; 
#     sequence = ID[[1]][2];
#     uniqueID = paste0(i, "_", ID[[1]][1])
#     # Format data for FASTA file appending
#     AA.formatted = BStringSet(sequence)
#     names(AA.formatted) = uniqueID
#     # Save data in list
#     sequenceList[[species]][uniqueID] = AA.formatted
#   }
# }
# 
# ## Write a FASTA file from the list
# for( k in 1:length(sequenceList) ){
#   WriteXStringSet(sequenceList[[k]], filepath=paste0(k,".fa"))
# }

## EOF