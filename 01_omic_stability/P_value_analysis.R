setwd('/Users/raphaelsinclair/Desktop/MSc Health Data Analytics - IC/HDA/SPH028 - Computational Epidemiology/Project/')

# Packages
library(dplyr)
library(lme4)

# Load p-value data
stab_pro <- readRDS('Results/stab_pro.rds')
stab_methyl <- readRDS('Results/stab_methyl.rds')
stab_trans <- readRDS('Results/stab_trans.rds')
stab_meta <- readRDS('Results/stab_meta.rds')

##### APPLY MULTIPLE CORRECTIONS TO P-VALUES #####
# We must adjust the p-value to account for the false positives that will appear due to having so many
# hypothesis tests to carry out. It is not clear which method is best in this case, so we will
# investigate a few


# Bonferroni multiple corrections
pro_bonf <- p.adjust(stab_pro, method = "bonferroni")
methyl_bonf <- p.adjust(stab_methyl, method = "bonferroni")
trans_bonf <- p.adjust(stab_trans, method = "bonferroni")
meta_bonf <- p.adjust(stab_meta, method = "bonferroni")

# Benjamini-Hochbegrg multiple corrections
pro_bh <- p.adjust(stab_pro, method = "BH")
methyl_bh <- p.adjust(stab_methyl, method = "BH")
trans_bh <- p.adjust(stab_trans, method = "BH")
meta_bh <- p.adjust(stab_meta, method = "BH")

# False Discovery Rate multiple corrections
pro_fdr <- p.adjust(stab_pro, method = "fdr")
methyl_fdr <- p.adjust(stab_methyl, method = "fdr")
trans_fdr <- p.adjust(stab_trans, method = "fdr")
meta_fdr <- p.adjust(stab_meta, method = "fdr")

# Append rownames
names(pro_bonf) <- names(pro_bh) <- names(pro_fdr) <- names(stab_pro)
names(methyl_bonf) <- names(methyl_bh) <- names(methyl_fdr) <- names(stab_methyl)
names(trans_bonf) <- names(trans_bh) <- names(trans_fdr) <- names(stab_trans)
names(meta_bonf) <- names(meta_bh) <- names(meta_fdr) <- names(stab_meta)

