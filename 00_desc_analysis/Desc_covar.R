setwd('/Users/raphaelsinclair/Desktop/MSc Health Data Analytics - IC/HDA/SPH028 - Computational Epidemiology/Project/')

# Packages
library(dplyr)
library(lme4)
library(corrplot)
library(pheatmap)

# Load data
covars <- readRDS('Data/Covariates.rds')
# exposures <- readRDS('Data/Exposures.rds')

# Calculate percentage of missing data for each column and create list of column names
# with too much missing data (for removal, if necessary)
prop_na <- as.data.frame(apply(covars,2, function(x){sum(is.na(x))/length(x)})*100)
colnames(prop_na) <- c('Percentage NA')

missing_cols <- NULL
missing_tol <- 20
for (row_ in 1:nrow(prop_na)) {
  if (prop_na[row_,] >= missing_tol) {
    missing_cols <- c(missing_cols, rownames(prop_na)[row_])
  }
}

# Compute correlation matrix (Note that correlation cannot be calculated for 
# unordered categorical variables)
temp <- subset(covars, select = -c(subjectidp, gender, season, session, city, date, id, plate, isolation, labeling, hybridization, chip, position))

temp <- cor(temp, use = "pairwise.complete.obs")
corrplot(temp, order = 'hclust')
pheatmap(temp)
