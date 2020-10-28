setwd('/Users/raphaelsinclair/Desktop/MSc Health Data Analytics - IC/HDA/SPH028 - Computational Epidemiology/Project/')

# Packages
library(dplyr)
library(lme4)
library(corrplot)
library(pheatmap)

# Load data
covars <- readRDS('Data/Covariates.rds')
proteins <- readRDS('Data/Proteins.rds')
transcripts <- readRDS('Data/Transcripts.rds')
metabolites <- readRDS('Data/Metabolites.rds')
methylation <- readRDS('Data/Methylation.rds')

# exposures <- readRDS('Data/Exposures.rds')

# Calculate percentage of missing data for each column and create list of column names
# with too much missing data (for removal, if necessary)
pro_na <- as.data.frame(apply(covars,2, function(x){sum(is.na(x))/length(x)})*100)
trans_na <- as.data.frame(apply(covars,2, function(x){sum(is.na(x))/length(x)})*100)
meta_na <- as.data.frame(apply(covars,2, function(x){sum(is.na(x))/length(x)})*100)
methyl_na <- as.data.frame(apply(covars,2, function(x){sum(is.na(x))/length(x)})*100)
colnames(pro_na) <- colnames(trans_na) <- colnames(meta_na) <- colnames(methyl_na) <- c('Percentage NA')

data_frames <- c('pro_na', 'trans_na', 'meta_na', 'methyl_na')

missing_cols <- NULL
missing_tol <- 10     # 10% threshold
for (df in data_frames) {
for (row_ in 1:nrow(paste0(df))) {
  if (paste0(df, '[', row_, ',]') >= missing_tol) {
    missing_cols <- c(missing_cols, row_)
    df <- paste0(df, '[-',missing_cols, ',]')
  }
}
}

# Compute correlation matrix (Note that correlation cannot be calculated for 
# unordered categorical variables)
temp <- subset(covars, select = -c(subjectidp, gender, season, session, city, date, id, plate, isolation, labeling, hybridization, chip, position))

temp <- cor(temp, use = "pairwise.complete.obs")
corrplot(temp, order = 'hclust')
pheatmap(temp)
