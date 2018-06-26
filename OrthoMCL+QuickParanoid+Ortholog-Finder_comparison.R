# Sam Old
# Compare Analysis of OrthoMCL, QuickParanoid, and Ortholog-Finder in a biologically meaningful
#   manner. Simulated Phylogenies are input and scored as a whole. {Perfect, split-in: 2,3,4,5,6+}
# Input parsing loads output of each program and homogenises the data formats

## Input file locations
phy_string = "C:/Users/samol/OneDrive/MASTERS/Data/Simulated_PHY/2017.10.30/"
orthomcl_string = "C:/Users/samol/OneDrive/MASTERS/Data/orthomcl/2017.10.30/groups.txt"
ortholog_finder_string = "C:/Users/samol/OneDrive/MASTERS/Data/ortholog-finder/2017.10.30/orthodata.txt"
quick_paranoid_string = "C:/Users/samol/OneDrive/MASTERS/Data/QuickParanoid_18MAR19/result.txt"

#########################
## Input Parsing Phase ##
#########################

## Read OrthoMCL output
omcl_df = read.table(orthomcl_string, fill=T, row.names=1, na.strings = c(""))
omcl_df = data.frame(lapply(omcl_df, as.character), stringsAsFactors = FALSE)
# Remove pipe info
for (omcl_nrow in 1:nrow(omcl_df)) {
  for (omcl_ncol in 1:ncol(omcl_df)){
    if(!is.na(omcl_df[omcl_nrow, omcl_ncol])){
      omcl_df[omcl_nrow, omcl_ncol] = strsplit(omcl_df[omcl_nrow, omcl_ncol], "\\|")[[1]][2]
}}}


## Read Ortholog-Finder output and reformat to be similar to orthomcl
f_temp = read.table(ortholog_finder_string)
f_temp = data.frame(lapply(f_temp, as.character), stringsAsFactors = FALSE)
sorted_DF = data.frame(0)
# For each row of the orthodata.txt
for(a in 1:nrow(f_temp)) {
  # If a line matches txt <<< then extract group number and use to sort new col
  if(length(grep("<<<", f_temp[a,1]))==1) {
    group_number = as.numeric(strsplit(f_temp[a,1], "<<<")[[1]][2])
    new_column_counter = 1;
  }
  # fill the new dataframe with information in correct order
  sorted_DF[new_column_counter, group_number] = f_temp[a,1]
  new_column_counter = new_column_counter + 1
}
# Remove the sequence data and rename columns accordingly
colnames(sorted_DF) = sorted_DF[1,]
to_delete = seq(1, nrow(sorted_DF), 2)
of_df = sorted_DF[-to_delete,]; rownames(of_df) = 1:nrow(of_df); 
of_df = as.data.frame(t(of_df)); rm(f_temp); rm(sorted_DF);
of_df = data.frame(lapply(of_df, as.character), stringsAsFactors = FALSE)
# Remove piped info
for (of_nrow in 1:nrow(of_df)) {
  for (of_ncol in 1:ncol(of_df)){
    if(!is.na(of_df[of_nrow, of_ncol])){
      of_df[of_nrow, of_ncol] = strsplit(of_df[of_nrow, of_ncol], "\\|")[[1]][2]
}}}

## Read Quick Paranoid output
quick_paranoid_string = "C:/Users/samol/OneDrive/MASTERS/Data/QuickParanoid_18MAR19/result.txt"
g = read.table(quick_paranoid_string, fill=T)
# this allows non-full columns to be added to DF:
g = data.frame(lapply(g, as.character), stringsAsFactors = FALSE) 
colnames(g) = c("clusterID", "species", "gene", "is_seed_ortholog", "confidence score", "species_in_cluster", "tree_conflict")
g[,7] = do.call(paste0, g[c(7,8)]); g = g[,1:7];
qp_df = data.frame(0); last_id = 0; qp_incrementer = cluster_incrementer = 0
for (id in g$clusterID) {
  qp_incrementer = qp_incrementer + 1
  if(id == last_id) {
    cluster_incrementer = cluster_incrementer + 1
    qp_df[id, cluster_incrementer] = toString(g$gene[qp_incrementer])
  } else {
    last_id = id
    cluster_incrementer = 1
    qp_df[id, cluster_incrementer] = toString(g$gene[qp_incrementer])
  }
}; rm(g)

#################################
## Statistics Generating Phase ##
#################################

# Objective: One Gene family at a time, score how it performs within each of the three programs:
#   1. Ortholog-Finder    2. QuickParanoid    3. OrthoMCL

library(phylotools)
phylip_files = list.files(phy_string)
phy_scores = data.frame(matrix(data = 0, nrow = length(phylip_files), ncol = 4))
colnames(phy_scores) = c("omcl", "qp", "of", "max")

# For each run##
for (phy in phylip_files) {
  
  # Extract number+1 for indexing later 
  #### !!!!!!! Will fuck up later if runs don't start at run00 :( :( :(
  phy_n = as.numeric(unlist(strsplit(phy, "run"))[2])+1
  
  # Load run## phy file
  current_phy = read.phy(paste0(phy_string, phy))
  
  # What is the number of sequences in this phy file?
  expected_matches = as.numeric(strsplit(current_phy[1], " ")[[1]][2])
  phy_scores[phy_n, "max"] = expected_matches
  
  # Create temp storage files 
  matches = vector(mode="character", length = expected_matches)
  omcl_matches = qp_matches = of_matches = matches
  
  # Identify all matches for current phy file against omlc, of, & qp dataframes
  for (phy_name in 1:expected_matches) {
    matches[phy_name] = paste0(phy,"_",strsplit(current_phy, " ")[[phy_name+1]][1])
    omcl_matches[phy_name] = which(omcl_df == matches[phy_name], arr.ind = T)[1]
    qp_matches[phy_name] = which(qp_df == matches[phy_name], arr.ind = T)[1]
    of_matches[phy_name] = which(of_df == matches[phy_name], arr.ind = T)[1]
    # Count the number of groups each phy file is spread amongst
    if (sum(is.na(omcl_matches)) == 0) {
      phy_scores[phy_n, "omcl"] = length(unique(omcl_matches))
    } else {
      phy_scores[phy_n, "omcl"] = sum(is.na(omcl_matches)) - 1 + 
                                  length(unique(omcl_matches)) 
    }
    if (sum(is.na(qp_matches)) == 0) {
      phy_scores[phy_n, "qp"] = length(unique(qp_matches))
    } else {
      phy_scores[phy_n, "qp"] = sum(is.na(qp_matches)) - 1 + 
        length(unique(qp_matches)) 
    }
    if (sum(is.na(of_matches)) == 0) {
      phy_scores[phy_n, "of"] = length(unique(of_matches))
    } else {
      phy_scores[phy_n, "of"] = sum(is.na(of_matches)) - 1 + 
        length(unique(of_matches)) 
    }
  }
}


############## GRAPHICAL OUTPUT ##################

# Scores: 1 (PERFECT)
#         2 (split into two groups) (including NAs an individual groups)
#         3 (split into three groups)
#         ...
#         6 (split into six+ groups)
#         7 (Did not match anything)

groups_number = matrix(nrow=7, ncol =3, data=0); 
colnames(groups_number) = colnames(phy_scores[,1:3])
for ( coln in colnames(phy_scores[,1:3])) {
  for ( rown in 1:nrow(phy_scores) ) {
    match_score = phy_scores[rown,coln]
    if(match_score < 6) {
      groups_number[match_score,coln] = groups_number[match_score,coln] + 1
    } else if (match_score < phy_scores[rown, "max"]) {
      groups_number[6,coln] = groups_number[6,coln] + 1
    } else if(match_score == phy_scores[rown, "max"]) {
      groups_number[7,coln] = groups_number[7,coln] + 1
    }
  }
}

colnames(groups_number) = c("OrthoMCL", "QuickParanoid", "Ortholog-Finder")
barplot(as.matrix(groups_number), 
        col = c("#1fea00", "#27a102", "#1c6000", "#002f06", "gray20", "black", "red"),
        legend.text = c("1", "2", "3", "4", "5", "6+", "0"),
        ylab = "Number of simulated orthologs\ families identified")

