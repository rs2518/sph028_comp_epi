# comp_epi
MSc Health Data Analytics w/ Machine Learning (2018/19)

SPH028 - Computational Epidemiology project

Project 2: Assessing OMICs stability
------------------------------------

Aim: Evaluate if OMICs profile vary in time, and what are the potential drivers of these changes

1. Identify sets of OMICs signals that are stable in time and Identify sets of OMICs signals that vary in time

  -	Run linear mixed models modelling ID and technical covariates as random effects to calculate variances within ID (between time points) whilst accounting for technical variability for each biomarker (univariate analysis). Note that linear models must be re-run where we obtain singular fits among the technical sources of variability (zero variance between technical variants) WITHOUT the technical variants are random effects. Forgetting this will alter the intra-class variance in ID as we are adjusting for a random effect which is not.
  -	Now calculate ICC = (intra-class variation due to ID)/(Total variability). This gives a measure of the correlation between the variance within ID. Note that an ICC of 0 means that there is NO correlation between the variance within ID and so we have INSTABILITY in the OMICs profile for that particular protein. An ICC of 1 means a that there is a perfect correlation structure in the variance within ID and so the profile is STABLE
  -	Need to determine whether observed ICCs are due to chance. To do this, we need to mimic a hypothesis test where we compare our observed ICCs to a ‘null distribution’ of ICCs. To create a null distribution, we decompose the correlation structure within IDs by using permutations. Therefore, the ICCs obtained will have no correlation structure within IDs and so we should have INSTABILITY (ICCs close to 0) – hence our null hypothesis for our observed ICCs is that the profiles on the original dataset are UNSTABLE.
  -	NOTE: We carry out multiple iterations of the permutation to completely decompose the (potential) associations contained in the permuted structure of the data. Over the numerous iterations, we obtain the null distribution and can therefore look towards comparing the obtained ICCs against the null distribution
  -	The p-value for each ICC is given by its position within the distribution of null ICCs – i.e. the proportion of null ICCs which are higher than the observed ICC. Note that only one tail is required because variability is non-negative.
  -	Adjust cut-off/threshold using multiple correction (Bonferroni etc) and compare results from different thresholds
  -	Analyse correlation between biomarkers to see whether stable/unstable signals are correlated
  
2. Explore the rate at which OMICs levels change

  -	Compare distributions of ICCs. Plot all ICC distributions (which should be on the same scale 0-1) on the same graph for visual representation. Note that we need to consider the rate of change from A to B (C may not be particularly relevant but should be examined)
  -	Examine change in the OMIC profile from session A to session B (then B to C)

3. Interpret the involved (metabolic) pathways and possible mechanisms

  -	Pathway analysis
  -	Look into network analysis

4. Evaluate if these changes are related to exposures / exposure changes

  -	Merge exposure with covariates to be able to account for exposure and confounders. As exposure is now involved in the model equation, we must account for confounders of the exposures
  -	Select fixed effects whilst ensuring shared variance between fixed effects is minimal
  -	Run 2 models in parallel (1 with exposure, 1 without) and compute difference in ICCs. Create distribution of differences in ICCs by permuting exposure. 
  -	Run linear mixed model on EACH exposure (plot each exposure separately)
  -	Compare the stable/unstable OMIC signals from the previous analysis with the signals with significant differences in ICCs (associated with the exposure). NOTE: ICCs are directly comparable here because the decomposition of variability in ICCs depends only on the random effects. The only other source of heterogeneity comes from differing sample sizes (i.e. samples that did not have exposure data are omitted)
  -	What percentage of biomarkers are related (overlap) with all 3 exposures?
  -	Analyse correlation among exposures (at each time point too).

5. Investigate the link between stable and variable signals across OMIC platforms

  -	Show gradient of stability (in terms of proportion of stable/unstable signals)

Student CIDs: 01559055; 01621068

Lecturers: Dr Marc Chadeau, Barbara Bordinier
