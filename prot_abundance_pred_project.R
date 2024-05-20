# 4/25/2024

# Mary Hathorn 
# Morgan Handojo 


# Our idea:

# Combine 3 random forests to try to predict protein abundance. The first two rf will be trained on the 2 most similar avg ribo cell lines
# to the new cell line because ribo has the highest correlation to protein abundance. The 3rd one be trained across all cell lines
# in order to still retain some of characteristics of that data.
#


# Overall Project Steps:

# 1. load in data, scale ribo and rna information using counts-per-million and log transform to make it more linear, 
#    as well as extract necessary features from mrna sequence data. From  mrna seq we extracted the following features:
#      - UTR5 length
#      - # of start codons in UTR5 (basically how many start codons before the actual one)
#      - UTR3 length
#      - CDS length
#      - GC content of the CDS


# 2. We made divided data frames based on cell line origin, as well as one for the mystery 7th cell line data.
#    These data frames have avg ribo and rna, as well as the first and third quadrants for ribo and rna, and the seq data.
#      - avg ribo and avg rna across cell-line specific experiments
#      - 25th percentile (q1) and 75th percentile (q3) for ribo and rna
#      - sequence data mentioned above
#
#    We also made a combined df that has all 6 provided cell lines.


# 3. Using PCA for dimensionality reduction, we use cosine similarity to calculate which cell lines (A549, HeLa, etc)
#    are the most similar to our new mystery cell line for ribo means. We then extract those names and train 2 models based
#    on that those 2 data frames.


# 4. Train the models.


# 5. Predict protein abundance using the myster data frame on each random forest. We will
#    take the average of these predictions to create our final predition, which will be
#    stored in the variable myst_final_pred.




# this code takes about 8 minutes to run


# install and load the necessary packages
install.packages("stringr")
library(stringr)

install.packages("randomForest")
library(randomForest)



# change path to wherever the file is in your computer

#rna_cl_train.rda
load("/your_path/rna_cl_train.rda")

#ribo_cl_train.rda
load("/your_path/ribo_cl_train.rda")

#prot_train.rda
load("/your_path/prot_train.rda")

#human_infor_train.csv
human_infor <- read.csv("your_path/human_infor_train.csv")
# Removing first column because it's useless
human_infor <- human_infor[, -1]

#mRNA_sequences.fa.gz
mrna <- gzfile("/your_path/mRNA_sequences.fa.gz", "r")

#RegionLengths.bed
mrna_info <- read.table("/your_path/RegionLengths.bed")


# mystery rna abundance goes here

load("")
myst_rna_filename <- load("")
myst_rna <- get(myst_rna_filename)

#mystery ribo abundance goes here

load("")
myst_ribo_filename <- load("")
myst_ribo <- get(myst_ribo_filename)




# NORMALIZING RNA AND RIBO DATA (COUNTS PER MILLION + LOG TRANSFORMATION)

counts_per_million <- function (count_matrix) {
  return (t(apply(count_matrix, 1, function(x){log(1000000*x/colSums(count_matrix)+1)})))
}

# scaling the A549, HeLa, etc
scaled_rna_cl <- counts_per_million(rna_cl_train)
scaled_ribo_cl <- counts_per_million(ribo_cl_train) 

# scaling the mystery 7th cell line data
scaled_test_rna <- counts_per_million(myst_rna)
scaled_test_ribo <- counts_per_million(myst_ribo)





# GATHERING DATA AND VARIABLES FROM MRNA SEQUENCES: 

# mrna seq has extra proteins that we don't actually want, so we want to isolate the 
# desired proteins, subset the mrna seq file, and then grab the protein, start index, end index,
# and class so we can make features for our trees
prot_of_interest <- rownames(prot_train)
mrna_info_subset <- mrna_info[mrna_info[, 1] %in% prot_of_interest, ]
colnames(mrna_info_subset) <- c("Protein", "Start_Index", "End_Index", "Classification")

# Making a DF from the mrna seq file
protein_names <- c()
mRNA_seqs <- c()
current_seq <- ""

# this while loop basically extracts information from each line 
while (TRUE) {
  line <- readLines(mrna, n = 1)
  if (length(line) == 0) {  
    break
  }
  
  if (substr(line, 1, 1) == ">") { 
    if (current_seq != "") {
      protein_names = c(protein_names, current_name)  
      mRNA_seqs = c(mRNA_seqs, current_seq)  
    }
    current_name = substring(line, 2)  
    current_seq = "" 
  } 
  else {
    current_seq = (line) 
  }
}

if (current_seq != "") {
  protein_names <- c(protein_names, current_name)
  mRNA_seqs <- c(mRNA_seqs, current_seq)
}

# put it all in a df
all_mrna <- data.frame(Protein = protein_names, Sequence = mRNA_seqs, stringsAsFactors = FALSE)

# selecting the proteins we are atually using xD
select_mrna <- all_mrna[all_mrna$Protein %in% prot_of_interest, ]

# adding subsequences to the MRNA data
# select_mrna holds the sequences while mrna_info_subset has like the classification, start and ends, etc
clean_mrna_data <- merge(select_mrna, mrna_info_subset, by = "Protein")

clean_mrna_data$Subsequence <- c("", nrow(clean_mrna_data))

extract_subsequence <- function(sequence, start, end) {
  start <- max(1, start)
  
  return(substr(sequence, start, end))
}

# this shoves the subsequence to the corresponding protein and classification
for (prot in seq_len(nrow(clean_mrna_data))) {
  clean_mrna_data$Subsequence[prot] <- extract_subsequence(clean_mrna_data$Sequence[prot],clean_mrna_data$Start_Index[prot],clean_mrna_data$End_Index[prot])
}

# removing extra PINX bc there's like more than one CDS, so we just selected the first one
clean_mrna_data <- clean_mrna_data[-(8356:8364),]

# subsets data to only have UTR5
UTR5_data <- subset(clean_mrna_data, Classification == "UTR5")
# gets length of the sequence
UTR5_data$UTR5_length <- nchar(UTR5_data$Subsequence)
#finds all start codons
UTR5_data$Start_Codons_UTR5 <- str_count(UTR5_data$Subsequence, fixed("ATG"))

# subsetting to only have UTR3
UTR3_data <- subset(clean_mrna_data, Classification == "UTR3")
# length of UTR3
UTR3_data$UTR3_length <- nchar(UTR3_data$Subsequence)

# subsetting to get CDS data
CDS_data <- subset(clean_mrna_data, Classification == "CDS")
# length of CDS
CDS_data$CDS_length <- nchar(CDS_data$Subsequence)
# GC content of CDS
CDS_data$GC_Content <- str_count(CDS_data$Subsequence, "[GC]") / nchar(CDS_data$Subsequence)

# combining all variables of into their respective larger DFs
UTR5_data <- UTR5_data[, c("Protein", "UTR5_length", "Start_Codons_UTR5")]
UTR3_data <- UTR3_data[, c("Protein", "UTR3_length")]
CDS_data <- CDS_data[, c("Protein", "CDS_length", "GC_Content")]

# merging everything together
mrna_char <- merge(UTR5_data, UTR3_data, by = "Protein", all = TRUE)
mrna_char <- merge(mrna_char, CDS_data, by = "Protein", all = TRUE)

# anything that is na is 0
# alternative is to impute na with the mean of the column
mrna_char[is.na(mrna_char)] <- 0


# change rownames to protein and remove the protein column
rownames(mrna_char) <- mrna_char$Protein
mrna_char$Protein <- NULL




# CREATING DATAFRAMES FOR EACH CELL LINE CONTAINING ALL NECESSARY DATA
create_cell_line_df <- function(cell_line_name, human_infor, scaled_ribo_cl, scaled_rna_cl, prot_train, mrna_char) {
  # get specific cell line
  cell_line_data <- subset(human_infor, cell_line == cell_line_name)
  
  # get column names from said line
  column_names <- cell_line_data[, 1]
  
  # grab the cell line specific experiments
  ribo_abundance <- scaled_ribo_cl[, column_names]
  
  rna_abundance <- scaled_rna_cl[, column_names, drop = FALSE]
  
  # need to convert to df otherwise single experiment cell lines will throw an error
  ribo_abundance <- as.data.frame(ribo_abundance)
  
  rna_abundance <- as.data.frame(rna_abundance)
  
  
  # calculate means
  if(ncol(ribo_abundance) > 1){
    ribo_mean <- rowMeans(ribo_abundance, na.rm = TRUE)
    rna_mean <- rowMeans(rna_abundance, na.rm = TRUE)
  } else{
    # if there is only one experiment
    # if you don't do this then K549 does some real funky stuff and doesn't work for
    # our combined data frame below LOL, for single lines, just set the only experiment as the "mean"
    ribo_mean <- setNames(as.data.frame(ribo_abundance[, 1, drop = FALSE]), "ribo_mean")
    rna_mean <- setNames(as.data.frame(rna_abundance[, 1, drop = FALSE]), "rna_mean")
  }
  
  
  # calculate row-wise first and third quartiles for ribo and rna abundances
  # for single-experiment cell lines these numbers are just the same as the abundance
  ribo_q1 <- apply(ribo_abundance, 1, quantile, probs = 0.25, na.rm = TRUE)
  ribo_q3 <- apply(ribo_abundance, 1, quantile, probs = 0.75, na.rm = TRUE)
  rna_q1 <- apply(rna_abundance, 1, quantile, probs = 0.25, na.rm = TRUE)
  rna_q3 <- apply(rna_abundance, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  # protein abundance data
  prot <- prot_train[cell_line_name]
  names(prot) <- "prot"
  
  # combine to df
  complete_df <- data.frame(
    ribo_mean,
    ribo_q1,
    ribo_q3,
    rna_mean,
    rna_q1,
    rna_q3,
    mrna_char,
    prot=prot
  )
  

  return(complete_df)
}


# make df for each cell line
A549_df <- create_cell_line_df("A549", human_infor, scaled_ribo_cl, scaled_rna_cl, prot_train, mrna_char)
HeLa_df <- create_cell_line_df("HeLa", human_infor, scaled_ribo_cl, scaled_rna_cl, prot_train, mrna_char)
HepG2_df <- create_cell_line_df("HepG2", human_infor, scaled_ribo_cl, scaled_rna_cl, prot_train, mrna_char)
K562_df <- create_cell_line_df("K562", human_infor, scaled_ribo_cl, scaled_rna_cl, prot_train, mrna_char)
MCF7_df <- create_cell_line_df("MCF7", human_infor, scaled_ribo_cl, scaled_rna_cl, prot_train, mrna_char)
U2OS_df <- create_cell_line_df("U2OS", human_infor, scaled_ribo_cl, scaled_rna_cl, prot_train, mrna_char)


# combined df with all cell lines
combined_df <- rbind(A549_df, HeLa_df, HepG2_df, K562_df, MCF7_df, U2OS_df)


# now let's make the dataframe for our myster protein
myst_protein_df_maker <- function(scaled_test_ribo, scaled_test_rna, mrna_char){
  # convert to df, mostly if its a single experiment cause otherwise it's a vector
  ribo_abundance <- as.data.frame(scaled_test_ribo)
  rna_abundance <- as.data.frame(scaled_test_rna)
  
  
  # calculate means
  if(ncol(ribo_abundance) > 1){
    ribo_mean <- rowMeans(ribo_abundance, na.rm = TRUE)
    rna_mean <- rowMeans(rna_abundance, na.rm = TRUE)
  } else{
    # in case the cell line for some reason only has 1 experiment
    # our combined data frame below LOL, for single lines, just set the only experiment as the "mean"
    # gotta use setNames otherwise it sets the column to "x_abundance" which is bad for the models later
    ribo_mean <- setNames(as.data.frame(ribo_abundance[, 1, drop = FALSE]), "ribo_mean")
    rna_mean <- setNames(as.data.frame(rna_abundance[, 1, drop = FALSE]), "rna_mean")
  }
  
  
  # calculate row-wise first and third quartiles for ribo and rna abundances
  # for single-experiment cell lines these numbers are just the same as the abundance
  ribo_q1 <- apply(ribo_abundance, 1, quantile, probs = 0.25, na.rm = TRUE)
  ribo_q3 <- apply(ribo_abundance, 1, quantile, probs = 0.75, na.rm = TRUE)
  rna_q1 <- apply(rna_abundance, 1, quantile, probs = 0.25, na.rm = TRUE)
  rna_q3 <- apply(rna_abundance, 1, quantile, probs = 0.75, na.rm = TRUE)
  
  done_df <- data.frame(
    ribo_mean,
    ribo_q1,
    ribo_q3,
    rna_mean,
    rna_q1,
    rna_q3,
    mrna_char
  )
  
  return(done_df)
}


# make df with myst protein to use for prediction
myst_df <- myst_protein_df_maker(scaled_test_ribo, scaled_test_rna, mrna_char)




# now trying to find the most similar cell lines to our mystery cell line

# using ribo to compare because it has the highest corr with protein abundance
other_points <- list(A549 = A549_df$ribo_mean, HeLa = HeLa_df$ribo_mean, HepG2 = HepG2_df$ribo_mean, K562 = K562_df$ribo_mean, MCF7 = MCF7_df$ribo_mean, U2OS = U2OS_df$ribo_mean)

# PCA doesn't work when SD is 0 so just remove 'em
remove_constant_columns <- function(data) {
  sds <- apply(data, 2, sd)
  non_constant_data <- data[, sds != 0, drop = FALSE]
  return(non_constant_data)
}

# Dimensionality reduction using PCA
perform_pca <- function(data) {
  # removes columns where sd is constant (no variation and causes calc issues)
  clean_data <- remove_constant_columns(data)
  
  pca <- prcomp(clean_data, center = TRUE, scale. = TRUE)
  
  # just use all the pca components, you can change this but when we ran it it only had 6 LOL
  num_pca_components = ncol(pca$x)
  
  return(pca$x[, 1:num_pca_components])
}

# cosine similarity used for distance, you could use Euclidean but apparently cosine is better for higher dim
cosine_similarity <- function(v1, v2) {
  # you can use the norm from the stats library but for some reason it wasn't working for me
  v1_norm <- sqrt(sum(v1^2))
  v2_norm <- sqrt(sum(v2^2))
  
  return(sum(v1 * v2) / (v1_norm * v2_norm))
}


# this ranks the similarities and spits out 
rank_similarities_highdim <- function(x, other_points) {
  # basically performs PCA on all the points including the myst, not sure if do.call is necessary though
  all_points <- do.call(cbind, c(list(x), other_points))
  reduced_points <- perform_pca(t(all_points))
  
  similarities <- sapply(2:nrow(reduced_points), function(i) {
    cosine_similarity(reduced_points[1, ], reduced_points[i, ])
  })
  
  # name vector so we can extract the actual names
  names(similarities) <- names(other_points)
  
  # descending order (higher is more similar)
  sorted_similarities <- sort(similarities, decreasing = TRUE)
  
  return(sorted_similarities)
}



similarities <- rank_similarities_highdim(myst_df$ribo_mean, other_points)

# we only want the top 2 similar cell lines - we found top two often had the best predictive power
two_similar_cell_lines = names(similarities[1:2])


# so we can access the dfs from the cell line names from the similarities we found
dataframes <- list(
  A549 = A549_df,
  HeLa = HeLa_df,
  HepG2 = HepG2_df, 
  K562 = K562_df, 
  MCF7 = MCF7_df,
  U2OS = U2OS_df
)


# this builds the cell line models
build_models <- function(cell_line_names) {

  # this was originally in a list bc we also wanted to store the df that was trained
  # for another idea, otherwise we would just store it in a vector to access the item
  models <- list()
  for (name in cell_line_names) {
    df <- dataframes[[name]]
      
    # train the model
    model <- randomForest(prot ~ ., df)
      
    # store the model and df in the list using the name as the key
    models[[name]] <- model
      
  }
    

  return(models)
  
}



most_similar_models <- build_models(cell_line_names = two_similar_cell_lines)


# training the large model
# this takes about 5 minutes to run, could be more or less for u LOLLL
big_model <- randomForest(prot ~ ., combined_df)

# access cell line models, storing in here mostly for better readability
cell_line_model_1 <- most_similar_models[[1]]
cell_line_model_2 <- most_similar_models[[2]]



# put everything in to a df
myst_data_predictions <- data.frame(
  pred1 = predict(big_model, myst_df),
  pred2 = predict(cell_line_model_1, myst_df),
  pred3 = predict(cell_line_model_2, myst_df)
)

# average all predictions
myst_final_pred <- rowMeans(myst_data_predictions)

final_message = "prediction is stored in myst_final_pred as a named numeric vector."

print(final_message)
print(myst_final_pred)

print(final_message)



