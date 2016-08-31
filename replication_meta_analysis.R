

# This is the second stage : organizing the data and analysing it 
# The file I am using here has entirely imputed SDs when it was possible

# installing and loading packages
# Make sure packages are installed first
library(compute.es)
library(meta)
library(metafor)
library(MAd)
library(ggplot2)
library(lawstat)
library(tidyr)
library(metagear)


## Importing the data
path <- file.path("MA_SD.csv")
MA <- read.csv(path) 

# Removing a few not useful (anymore) columns
MA$X <- NULL
MA$d_var <- NULL
MA$corr <- NULL
MA$excel_calculated_d <- NULL

#Fixing variable types
MA$study_ID <- as.character(MA$study_ID)
MA$long_cite <- as.character(MA$long_cite)
MA$FYI...inclusion.notes <- as.character(MA$FYI...inclusion.notes)
MA$expt_num <- as.character(MA$expt_num)
MA$expt_condition <- as.character(MA$expt_condition)
MA$notes <- as.character(MA$notes)

MA$diff_means <- MA$x_1 - MA$x_2



# Setting correlations 
ri_var = 0.5
MA$inputed_corr <- ifelse(!(MA$seed_ID == "wordWSeg") & !(MA$seed_ID == "speechWSeg"), ri_var, 0.641)
MA$inputed_corr <- ifelse(MA$participant_design == "between", 0, MA$inputed_corr)



# Calculating effect sizes 

# For now, I will ignore LTS because I haven't found a way to do a meta-analysis without standard deviations 
# In the mean time, I have discovered that it is impossible to do a meta-analysis without them, and there is not enough information to try and infer them
# So I will do the meta-analysis without using LTS


MA_sd <- subset(MA, MA$ES_measure == "sd_mean_diff")
MA_sd_b <- escalc(measure = "SMD", m1i = MA_sd$x_1, m2i = MA_sd$x_2, sd1i = MA_sd$SD_1, sd2i = MA_sd$SD_2, n1i = MA_sd$n_1, n2i = MA_sd$n_2, subset = MA_sd$participant_design == "between", var.names = c("Effect_size", "ES_var"), data = MA_sd, append = TRUE)
MA_sd_w <- escalc(measure = "SMCC", m1i = MA_sd$x_1, m2i = MA_sd$x_2, sd1i = MA_sd$SD_1, sd2i = MA_sd$SD_2, ni = MA_sd$n_1, ri = MA_sd$inputed_corr, var.names = c("Effect_size", "ES_var"), subset = MA_sd$participant_design == "within_two", append = TRUE, data = MA_sd)



MA_binom <- subset(MA, MA$ES_measure == "binom")
MA_binom <- escalc(measure = "PR", xi = MA_binom$x_1, mi = MA_binom$x_2, var.names = c("Effect_size", "ES_var"), append = TRUE, data = MA_binom)

# Merging everything back
MA_ES <- rbind(MA_sd_b, MA_sd_w, MA_binom)
MA_ES <- MA_ES[order(MA_ES$seed_ID, MA_ES$Study_type_.OR.R., MA_ES$year), ]

# Now we need to agregate studies that had repeated measures (used the same group of infants for several measure). I am using the function in the MAd package, which follows textbook recommendations for the ES and ES_var, and I am getting the median for the other values
MA_ES_diff <- subset(MA_ES, MA_ES$same_infant == "")
MA_ES_same <- subset(MA_ES, MA_ES$same_infant != "")

aggr_es <- agg(same_infant, Effect_size, ES_var, data = MA_ES_same)
aggr_es <- aggr_es[order(aggr_es$id),]

f <- function(x){
  if(class(x)=="numeric"){
    return(median(x))
  } else {
    return(x[1])
  }
}

aggr_same <- aggregate(MA_ES_same, by = list(infant_id = MA_ES_same$same_infant), FUN = f)
aggr_same <- aggr_same[order(aggr_same$infant_id),]
stopifnot(aggr_same$infant_id == aggr_es$id)
aggr_same$Effect_size <- aggr_es$es
aggr_same$ES_var <- aggr_es$var
aggr_same$expt_condition <- "aggregated_measure"
aggr_same$group_name_2 <- "aggregated_measure"
aggr_same$infant_id <- NULL

# Putting back together will only unique infants groups 

MA1 <- rbind(aggr_same, MA_ES_diff)
MA1 <- MA1[order(MA1$seed_ID, MA1$Study_type_.OR.R., MA1$year), ]

summary(MA1$p_inc)
# Descriptive plot 
plot(MA1$Effect_size, MA1$ES_var, xlab = "Effect size", ylab = "Pooled SD")
length(MA1$ES_var < 0)

# Creating a data.frame for each meta-analysis
FBl <- subset(MA1, seed_ID == "FBl")
NDi <- subset(MA1, seed_ID == "NDis")
wWS <- subset(MA1, seed_ID == "wordWSeg")
sWS <- subset(MA1, seed_ID == "speechWSeg")
HTS <- subset(MA1, seed_ID == "HTScram")
sIS <- subset(MA1, seed_ID == "samInfS")
pIS <- subset(MA1, seed_ID == "popInfS")
VST <- subset(MA1, seed_ID == "VisStats")
# LTS <- subset(MA1, seed_ID == "LTScram")
IACT <- subset(MA1, seed_ID == "IntAct")
GRS <- subset(MA1, seed_ID == "GrSoc")


##### Calculating Summary Effect Sizes  ##### 

#1  False Belief
mm_FBl <- rma(yi = FBl$Effect_size, vi = FBl$ES_var, data = FBl)

ES_mm_FBl <- data.frame(name = "False Belief", Type = "Diff btw means", ES = mm_FBl$b, pval = round(mm_FBl$pval, 3))

#2 Number discrimination 
mm_NDi <- rma(yi = NDi$Effect_size, vi = NDi$ES_var, data = NDi)
ES_mm_NDi <- data.frame(name = "Number Discrimination", Type = "Diff btw means", ES = mm_NDi$b, pval = round(mm_NDi$pval,3))

#3 Word first Word Segmentation
mm_wWS <- rma(yi = wWS$Effect_size, vi = wWS$ES_var, data = wWS)
ES_mm_wWS <- data.frame(name = "w Word Seg", Type = "Diff btw means", ES = mm_wWS$b, pval = round(mm_wWS$pval,3))

#4 Speech first Word Segmentation
mm_sWS <- rma(yi = sWS$Effect_size, vi = sWS$ES_var, data = sWS)
ES_mm_sWS <- data.frame(name = "sp Word Seg", Type = "Diff btw means", ES = mm_sWS$b, pval = round(mm_sWS$pval,3))
mm_sWS
#5 Head Turn to Scrambled Faces
mm_HTS <- rma(yi = HTS$Effect_size, vi = HTS$ES_var, data = HTS)
ES_mm_HTS <- data.frame(name = "Head Turn", Type = "Diff btw means", ES = mm_HTS$b, pval = round(mm_HTS$pval,3))

#6 Sample first Inferential Statistics
mm_sIS <- rma(yi = sIS$Effect_size, vi = sIS$ES_var, data = sIS)
ES_mm_sIS <- data.frame(name = "s Inf Stats", Type = "Diff btw means", ES = mm_sIS$b, pval = round(mm_sIS$pval,3))

#7 Population first Inferential Statistics 
mm_pIS <- rma(yi = pIS$Effect_size, vi = pIS$ES_var, data = pIS)
ES_mm_pIS <- data.frame(name = "pop Inf Stats", Type = "Diff btw means", ES = mm_pIS$b, pval = round(mm_pIS$pval,3))

# 8 Visual statistics
mm_VST <- rma(yi = VST$Effect_size, vi = VST$ES_var, data = VST)
ES_mm_VST <- data.frame(name = "Visual Stats", Type = "Diff btw means", ES = mm_VST$b, pval = round(mm_VST$pval, 3))

# 9 Intential Action expected to be rational 
mm_IACT <- rma(yi = IACT$Effect_size, vi = IACT$ES_var, data = IACT)
ES_mm_IACT <- data.frame(name = "Int. Action", Type = "Diff btw means", ES = mm_IACT$b, pval = round(mm_IACT$pval, 3))

# 10 Preference for social individuals (grabing)
mm_GRS <- rma(yi = GRS$Effect_size, vi = GRS$ES_var, data = GRS)
ES_mm_GRS <- data.frame(name = "Preferring Social", Type = "Binomial", ES = mm_GRS$b, pval = round(mm_GRS$pval, 3))

ES_mm_total <- rbind(ES_mm_FBl, ES_mm_NDi, ES_mm_wWS, ES_mm_sWS, ES_mm_HTS, ES_mm_sIS, ES_mm_pIS, ES_mm_VST, ES_mm_IACT, ES_mm_GRS)
# ES_mm_total
# Uncomment to print table.

## Do experiments replicate the original results ? ####
# To do this, we will compare confidence intervals of the original studies with effect size of the replications. 
# To compute confidence intervals : either read it on a forest plot or from a summary effect size when there are several studies

# FBl

#1  False Belief
mm_FBl_OR <- rma(yi = FBl$Effect_size, vi = FBl$ES_var, data = FBl, subset = FBl$Study_type_.OR.R. == "OR")
mm_FBl_REP <- rma(yi = FBl$Effect_size, vi = FBl$ES_var, data = FBl, subset = FBl$Study_type_.OR.R. != "OR")
FBl_rep <- data.frame(estimate = c(coef(mm_FBl_OR), coef(mm_FBl_REP)), stderror = c(mm_FBl_OR$se, mm_FBl_REP$se),
                      meta = c("original","a-total"), tau2 = round(c(mm_FBl_OR$tau2, mm_FBl_REP$tau2),3))


FBl_rep_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = FBl_rep, digits =3)

REP_total <- data.frame(OR_ES = round(mm_FBl_OR$b, 3), rep_diff = round(FBl_rep_comp$b[2], 3), diff_pval = round(FBl_rep_comp$pval[2], 3))


#2 Number discrimination 
mm_NDi_OR <- rma(yi = NDi$Effect_size, vi = NDi$ES_var, data = NDi, subset = NDi$Study_type_.OR.R. == "OR")
mm_NDi_REP <- rma(yi = NDi$Effect_size, vi = NDi$ES_var, data = NDi, subset = NDi$Study_type_.OR.R. != "OR")
NDi_rep <- data.frame(estimate = c(coef(mm_NDi_OR), coef(mm_NDi_REP)), stderror = c(mm_NDi_OR$se, mm_NDi_REP$se),
                      meta = c("original","a-total"), tau2 = round(c(mm_NDi_OR$tau2, mm_NDi_REP$tau2),3))


NDi_rep_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = NDi_rep, digits =3)


REP_total <- rbind(REP_total, c(round(mm_NDi_OR$b, 3), round(NDi_rep_comp$b[2],3), round(NDi_rep_comp$pval[2],3)))

#3 Word first Word Segmentation
mm_wWS_OR <- rma(yi = wWS$Effect_size, vi = wWS$ES_var, data = wWS, subset = wWS$Study_type_.OR.R. == "OR")
mm_wWS_REP <- rma(yi = wWS$Effect_size, vi = wWS$ES_var, data = wWS, subset = wWS$Study_type_.OR.R. != "OR")
wWS_rep <- data.frame(estimate = c(coef(mm_wWS_OR), coef(mm_wWS_REP)), stderror = c(mm_wWS_OR$se, mm_wWS_REP$se),
                      meta = c("original","a-total"), tau2 = round(c(mm_wWS_OR$tau2, mm_wWS_REP$tau2),3))


wWS_rep_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = wWS_rep, digits =3)

REP_total <- rbind(REP_total, c(round(mm_wWS_OR$b, 3), round(wWS_rep_comp$b[2],3), round(wWS_rep_comp$pval[2],3)))

#4 Speech first Word Segmentation
mm_sWS_OR <- rma(yi = sWS$Effect_size, vi = sWS$ES_var, data = sWS, subset = sWS$Study_type_.OR.R. == "OR")

mm_sWS_REP = rma(yi = sWS$Effect_size, vi = sWS$ES_var, data = sWS, subset = sWS$Study_type_.OR.R. !="OR")

sWS_rep <- data.frame(estimate = c(coef(mm_sWS_OR), coef(mm_sWS_REP)), stderror = c(mm_sWS_OR$se, mm_sWS_REP$se),
                      meta = c("original","a-total"), tau2 = round(c(mm_sWS_OR$tau2, mm_sWS_REP$tau2),3))

sWS_rep_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = sWS_rep, digits =3)
sWS_rep_comp
REP_total <- rbind(REP_total, c(round(mm_sWS_OR$b, 3), round(sWS_rep_comp$b[2],3), round(sWS_rep_comp$pval[2],3)))

#5 Head Turn to Scrambled Faces
mm_HTS_OR <- rma(yi = HTS$Effect_size, vi = HTS$ES_var, data = HTS, subset = HTS$Study_type_.OR.R. == "OR")
mm_HTS_REP <- rma(yi = HTS$Effect_size, vi = HTS$ES_var, data = HTS, subset = HTS$Study_type_.OR.R. != "OR")
HTS_rep <- data.frame(estimate = c(coef(mm_HTS_OR), coef(mm_HTS_REP)), stderror = c(mm_HTS_OR$se, mm_HTS_REP$se),
                      meta = c("original","a-total"), tau2 = round(c(mm_HTS_OR$tau2, mm_HTS_REP$tau2),3))

HTS_rep_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = HTS_rep, digits =3)
HTS_rep_comp
REP_total <- rbind(REP_total, c(round(mm_HTS_OR$b, 3), round(HTS_rep_comp$b[2],3), round(HTS_rep_comp$pval[2],3)))

#6 Sample first Inferential Statistics
mm_sIS_OR <- rma(yi = sIS$Effect_size, vi = sIS$ES_var, data = sIS, subset = sIS$Study_type_.OR.R. == "OR")
mm_sIS_REP <- rma(yi = sIS$Effect_size, vi = sIS$ES_var, data = sIS, subset = sIS$Study_type_.OR.R. != "OR")
sIS_rep <- data.frame(estimate = c(coef(mm_sIS_OR), coef(mm_sIS_REP)), stderror = c(mm_sIS_OR$se, mm_sIS_REP$se),
                      meta = c("original","a-total"), tau2 = round(c(mm_sIS_OR$tau2, mm_sIS_REP$tau2),3))

sIS_rep_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = sIS_rep, digits =3)
sIS_rep_comp
REP_total <- rbind(REP_total, c(round(mm_sIS_OR$b, 3), round(sIS_rep_comp$b[2],3), round(sIS_rep_comp$pval[2],3)))

#7 Population first Inferential Statistics 
mm_pIS_OR <- rma(yi = pIS$Effect_size, vi = pIS$ES_var, data = pIS, subset = pIS$Study_type_.OR.R. == "OR")
mm_pIS_REP <- rma(yi = pIS$Effect_size, vi = pIS$ES_var, data = pIS, subset = pIS$Study_type_.OR.R. != "OR")
pIS_rep <- data.frame(estimate = c(coef(mm_pIS_OR), coef(mm_pIS_REP)), stderror = c(mm_pIS_OR$se, mm_pIS_REP$se),
                      meta = c("original","a-total"), tau2 = round(c(mm_pIS_OR$tau2, mm_pIS_REP$tau2),3))

pIS_rep_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = pIS_rep, digits =3)
pIS_rep_comp
REP_total <- rbind(REP_total, c(round(mm_pIS_OR$b, 3), round(pIS_rep_comp$b[2],3), round(pIS_rep_comp$pval[2],3)))

# 8 Visual statistics
mm_VST_OR <- rma(yi = VST$Effect_size, vi = VST$ES_var, data = VST, subset = VST$Study_type_.OR.R. == "OR")
mm_VST_REP <- rma(yi = VST$Effect_size, vi = VST$ES_var, data = VST, subset = VST$Study_type_.OR.R. !="OR")
VST_rep <- data.frame(estimate = c(coef(mm_VST_OR), coef(mm_VST_REP)), stderror = c(mm_VST_OR$se, mm_VST_REP$se),
                      meta = c("original","a-total"), tau2 = round(c(mm_VST_OR$tau2, mm_VST_REP$tau2),3))

VST_rep_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = VST_rep, digits =3)
VST_rep_comp
REP_total <- rbind(REP_total, c(round(mm_VST_OR$b, 3), round(VST_rep_comp$b[2],3), round(VST_rep_comp$pval[2],3)))

# 9 Intential Action expected to be rational 
mm_IACT_OR <- rma(yi = IACT$Effect_size, vi = IACT$ES_var, data = IACT, subset = IACT$Study_type_.OR.R. == "OR")
mm_IACT_REP <- rma(yi = IACT$Effect_size, vi = IACT$ES_var, data = IACT, subset =IACT$Study_type_.OR.R. !="OR")
IACT_rep <- data.frame(estimate = c(coef(mm_IACT_OR), coef(mm_IACT_REP)), stderror = c(mm_IACT_OR$se, mm_IACT_REP$se),
                       meta = c("original","a-total"), tau2 = round(c(mm_IACT_OR$tau2, mm_IACT_REP$tau2),3))

IACT_rep_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = IACT_rep, digits =3)
IACT_rep_comp
REP_total <- rbind(REP_total, c(round(mm_IACT_OR$b, 3), round(IACT_rep_comp$b[2],3), round(IACT_rep_comp$pval[2],3)))

# 10 Preference for social individuals (grabing)
mm_GRS_OR <- rma(yi = GRS$Effect_size, vi = GRS$ES_var, data = GRS, subset = GRS$Study_type_.OR.R. == "OR")
mm_GRS_REP <- rma(yi = GRS$Effect_size, vi = GRS$ES_var, data = GRS, subset = GRS$Study_type_.OR.R. != "OR")
mm_GRS_OR
mm_GRS_REP
GRS_rep <- data.frame(estimate = c(coef(mm_GRS_OR), coef(mm_GRS_REP)), stderror = c(mm_GRS_OR$se, mm_GRS_REP$se),
                      meta = c("original","a-total"), tau2 = round(c(mm_GRS_OR$tau2, mm_GRS_REP$tau2),3))

GRS_rep_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = GRS_rep, digits =3)
GRS_rep_comp
REP_total <- rbind(REP_total, c(round(mm_GRS_OR$b, 3), round(GRS_rep_comp$b[2],3), round(GRS_rep_comp$pval[2],3)))

REP_total
ES_mm_total <- cbind(ES_mm_total, REP_total)
ES_mm_total



#### Decrease in effect size over time ?  #### 



# Metaregression to test for the decrease : using the year as a moderator
decrES_metareg_total <- rma(Effect_size, ES_var, mods = ~ year + seed_ID, data = MA1)
decrES_metareg_total




####### Testing the hypothesis ##### 
# This time we will use two methods of analysis for each hypothesis : 
# the first one will compare accross the whole data-set 
# The second one will be used to ensure that the category with the more members doesn't favor some studies that might have higher/lower effect sizes, which might skew the results
# To do this, we will match the number of studies in the bigger category to the smaller one, using a procedure of random sampling 10 times so that there are as many experiments of each category sampled from a single seed study
# To test for bias, the assumption is that the factor causing bias leads to a higher effect size, so we will use a subgroup analysis
# To test for source of error, the assumption is that the factor causing error leads to a higher pooled variance. To compare them we will use the Brown-Forsythe test, which is used to check for equality of variance and is more robust that the Levene test

###### Are the methodological factors a source of bias ? ####

# Analysis 1 : Subgroup analysis
###### Mode of presentation of the stimuli ####

mm_MA1_m <- rma(yi = MA1$Effect_size, vi = MA1$ES_var, data = MA1, subset=stimulus_type_.M.A.== "M")
mm_MA1_a <- rma(yi = MA1$Effect_size, vi = MA1$ES_var, data = MA1, subset=stimulus_type_.M.A.== "A")
MA1_stim <- data.frame(estimate = c(coef(mm_MA1_m), coef(mm_MA1_a)), stderror = c(mm_MA1_m$se, mm_MA1_a$se),
                       meta = c("manual","automatic"), tau2 = round(c(mm_MA1_m$tau2, mm_MA1_a$tau2),3))

mm_MA1_m
mm_MA1_a

MA1_stim_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = MA1_stim, digits =3)
MA1_stim_comp

#### Coding method of the infants reaction #### We will make simple and simply compare live to non-live

mm_MA1_l <- rma(yi = MA1$Effect_size, vi = MA1$ES_var, data = MA1, subset=coding_type_.A.L.H.R. == "L")
mm_MA1_notl <- rma(yi = MA1$Effect_size, vi = MA1$ES_var, data = MA1, subset=coding_type_.A.L.H.R. != "L")
MA1_cod <- data.frame(estimate = c(coef(mm_MA1_l), coef(mm_MA1_notl)), stderror = c(mm_MA1_l$se, mm_MA1_notl$se),
                       meta = c("live","a_notlive"), tau2 = round(c(mm_MA1_l$tau2, mm_MA1_notl$tau2),3))

mm_MA1_l
mm_MA1_notl

MA1_cod_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = MA1_cod, digits =3)
MA1_cod_comp

### Ratio of included infants ####
MA1ex <- MA1[complete.cases(MA1$p_ratio), ]

mm_MA1ex_high <- rma(yi = MA1ex$Effect_size, vi = MA1ex$ES_var, data = MA1ex, subset=p_ratio >= median(MA1ex$p_ratio))
mm_MA1ex_low <- rma(yi = MA1ex$Effect_size, vi = MA1ex$ES_var, data = MA1ex, subset=p_ratio < median(MA1ex$p_ratio))
MA1ex_rat <- data.frame(estimate = c(coef(mm_MA1ex_high), coef(mm_MA1ex_low)), stderror = c(mm_MA1ex_high$se, mm_MA1ex_low$se),
                      meta = c("high","a_low"), tau2 = round(c(mm_MA1ex_high$tau2, mm_MA1ex_low$tau2),3))

mm_MA1ex_high
mm_MA1ex_low

MA1_rat_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = MA1ex_rat, digits =3)
MA1_rat_comp

### Analysis 2 : matching samples with randomisation

# First step : how many of each ? 
# Automatic 132 ; Manual : 53
# Not Live : 76 ; Live : 109
# Low  :117 ; High :  


#### Stimulus type ####
ES_and_var <- c("Effect_size", "ES_var")
randomizing <- c(1:100)
summary(MA1$seed_ID)
# FBl  ## A : 3, M : 8

FBl_A <- subset(FBl, subset = stimulus_type_.M.A. == "A" )
FBl_M <- subset(FBl, subset = stimulus_type_.M.A. == "M")
nrow(FBl_M) > nrow(FBl_A)
FBl_stim_samp1_ES <- data.frame(y = 1:nrow(FBl_A))
FBl_stim_samp1_VAR <- data.frame(y = 1:nrow(FBl_A))

for (i in randomizing) {
  samp1_FBl <- FBl_M[sample(1:nrow(FBl_M), nrow(FBl_A), replace=FALSE), ES_and_var]
  samp1_FBl_ES <- samp1_FBl[,"Effect_size"]
  samp1_FBl_VAR <- samp1_FBl[,"ES_var"]
  FBl_stim_samp1_ES <- cbind(FBl_stim_samp1_ES, samp1_FBl_ES, deparse.level = 1)
  FBl_stim_samp1_VAR <- cbind(FBl_stim_samp1_VAR, samp1_FBl_VAR, deparse.level = 1)
}

FBl_stim_samp1_ES$"Rand_ES" <- rowMeans(FBl_stim_samp1_ES)
FBl_stim_samp1_VAR$"Rand_VAR" <- rowMeans(FBl_stim_samp1_VAR)
FBl_rand_M <- cbind(FBl_stim_samp1_ES$Rand_ES, FBl_stim_samp1_VAR$Rand_VAR)
colnames(FBl_rand_M) <- ES_and_var
FBl_rand_M <- cbind(FBl_rand_M, stimulus_type_.M.A. = "M", seed_ID = "FBl")


# NDi ## A : 2 ; M : 6 
NDi_A <- subset(NDi, subset = stimulus_type_.M.A. == "A" )
NDi_M <- subset(NDi, subset = stimulus_type_.M.A. == "M")

NDi_stim_samp1_ES <- data.frame(y = 1:nrow(NDi_A))
NDi_stim_samp1_VAR <- data.frame(y = 1:nrow(NDi_A))

for (i in randomizing) {
  samp1_NDi <- NDi_M[sample(1:nrow(NDi_M), nrow(NDi_A), replace=FALSE), ES_and_var]
  samp1_NDi_ES <- samp1_NDi[,"Effect_size"]
  samp1_NDi_VAR <- samp1_NDi[,"ES_var"]
  NDi_stim_samp1_ES <- cbind(NDi_stim_samp1_ES, samp1_NDi_ES, deparse.level = 1)
  NDi_stim_samp1_VAR <- cbind(NDi_stim_samp1_VAR, samp1_NDi_VAR, deparse.level = 1)
}

NDi_stim_samp1_ES$"Rand_ES" <- rowMeans(NDi_stim_samp1_ES)
NDi_stim_samp1_VAR$"Rand_VAR" <- rowMeans(NDi_stim_samp1_VAR)
NDi_rand_M <- cbind(NDi_stim_samp1_ES$Rand_ES, NDi_stim_samp1_VAR$Rand_VAR)
colnames(NDi_rand_M) <- ES_and_var
NDi_rand_M <- cbind(NDi_rand_M, stimulus_type_.M.A. = "M", seed_ID = "NDis")


# "wordWSeg" ## A : 47 , M : 0, so we keep 0 of each
# "speechWSeg" ## A : 36, m : 0, so we keep 0 of each

# "HTScram"  ### A : 3, M : 10
HTS_A <- subset(HTS, subset = stimulus_type_.M.A. == "A" )
HTS_M <- subset(HTS, subset = stimulus_type_.M.A. == "M")

HTS_stim_samp1_ES <- data.frame(y = 1:nrow(HTS_A))
HTS_stim_samp1_VAR <- data.frame(y = 1:nrow(HTS_A))

for (i in randomizing) {
  samp1_HTS <- HTS_M[sample(1:nrow(HTS_M), nrow(HTS_A), replace=FALSE), ES_and_var]
  samp1_HTS_ES <- samp1_HTS[,"Effect_size"]
  samp1_HTS_VAR <- samp1_HTS[,"ES_var"]
  HTS_stim_samp1_ES <- cbind(HTS_stim_samp1_ES, samp1_HTS_ES, deparse.level = 1)
  HTS_stim_samp1_VAR <- cbind(HTS_stim_samp1_VAR, samp1_HTS_VAR, deparse.level = 1)
}

HTS_stim_samp1_ES$"Rand_ES" <- rowMeans(HTS_stim_samp1_ES)
HTS_stim_samp1_VAR$"Rand_VAR" <- rowMeans(HTS_stim_samp1_VAR)
HTS_rand_M <- cbind(HTS_stim_samp1_ES$Rand_ES, HTS_stim_samp1_VAR$Rand_VAR)
colnames(HTS_rand_M) <- ES_and_var
HTS_rand_M <- cbind(HTS_rand_M, stimulus_type_.M.A. = "M", seed_ID = "HTScram")



# "samInfS" ## A : 1, M : 2
sIS_A <- subset(sIS, subset = stimulus_type_.M.A. == "A" )
sIS_M <- subset(sIS, subset = stimulus_type_.M.A. == "M")

sIS_stim_samp1_ES <- data.frame(y = 1:nrow(sIS_A))
sIS_stim_samp1_VAR <- data.frame(y = 1:nrow(sIS_A))

for (i in randomizing) {
  samp1_sIS <- sIS_M[sample(1:nrow(sIS_M), nrow(sIS_A), replace=FALSE), ES_and_var]
  samp1_sIS_ES <- samp1_sIS[,"Effect_size"]
  samp1_sIS_VAR <- samp1_sIS[,"ES_var"]
  sIS_stim_samp1_ES <- cbind(sIS_stim_samp1_ES, samp1_sIS_ES, deparse.level = 1)
  sIS_stim_samp1_VAR <- cbind(sIS_stim_samp1_VAR, samp1_sIS_VAR, deparse.level = 1)
}

sIS_stim_samp1_ES$"Rand_ES" <- rowMeans(sIS_stim_samp1_ES)
sIS_stim_samp1_VAR$"Rand_VAR" <- rowMeans(sIS_stim_samp1_VAR)
sIS_rand_M <- cbind(sIS_stim_samp1_ES$Rand_ES, sIS_stim_samp1_VAR$Rand_VAR)
colnames(sIS_rand_M) <- ES_and_var
sIS_rand_M <- cbind(sIS_rand_M, stimulus_type_.M.A. = "M", seed_ID = "samInfS")
sIS_rand_M


# "popInfS" ## A : 6, M : 4 
pIS_A <- subset(pIS, subset = stimulus_type_.M.A. == "A" )
pIS_M <- subset(pIS, subset = stimulus_type_.M.A. == "M")

pIS_stim_samp1_ES <- data.frame(y = 1:nrow(pIS_M))
pIS_stim_samp1_VAR <- data.frame(y = 1:nrow(pIS_M))

for (i in randomizing) {
  samp1_pIS <- pIS_A[sample(1:nrow(pIS_A), nrow(pIS_M), replace=FALSE), ES_and_var]
  samp1_pIS_ES <- samp1_pIS[,"Effect_size"]
  samp1_pIS_VAR <- samp1_pIS[,"ES_var"]
  pIS_stim_samp1_ES <- cbind(pIS_stim_samp1_ES, samp1_pIS_ES, deparse.level = 1)
  pIS_stim_samp1_VAR <- cbind(pIS_stim_samp1_VAR, samp1_pIS_VAR, deparse.level = 1)
}
pIS_stim_samp1_ES$"Rand_ES" <- rowMeans(pIS_stim_samp1_ES)
pIS_stim_samp1_VAR$"Rand_VAR" <- rowMeans(pIS_stim_samp1_VAR)
pIS_rand_A <- cbind(pIS_stim_samp1_ES$Rand_ES, pIS_stim_samp1_VAR$Rand_VAR)
colnames(pIS_rand_A) <- ES_and_var
pIS_rand_A <- cbind(pIS_rand_A, stimulus_type_.M.A. = "A", seed_ID = "popInfS")

# "VisStats" ## A : 20, M : 0, so we keep 0 of each

# "IntAct"## A : 6 ; M : 1
IACT_A <- subset(IACT, subset = stimulus_type_.M.A. == "A" )
IACT_M <- subset(IACT, subset = stimulus_type_.M.A. == "M")

IACT_stim_samp1_ES <- data.frame(y = 1:nrow(IACT_M))
IACT_stim_samp1_VAR <- data.frame(y = 1:nrow(IACT_M))

for (i in randomizing) {
  samp1_IACT <- IACT_A[sample(1:nrow(IACT_A), nrow(IACT_M), replace=FALSE), ES_and_var]
  samp1_IACT_ES <- samp1_IACT[,"Effect_size"]
  samp1_IACT_VAR <- samp1_IACT[,"ES_var"]
  IACT_stim_samp1_ES <- cbind(IACT_stim_samp1_ES, samp1_IACT_ES, deparse.level = 1)
  IACT_stim_samp1_VAR <- cbind(IACT_stim_samp1_VAR, samp1_IACT_VAR, deparse.level = 1)
}
IACT_stim_samp1_ES$"Rand_ES" <- rowMeans(IACT_stim_samp1_ES)
IACT_stim_samp1_VAR$"Rand_VAR" <- rowMeans(IACT_stim_samp1_VAR)
IACT_rand_A <- cbind(IACT_stim_samp1_ES$Rand_ES, IACT_stim_samp1_VAR$Rand_VAR)
colnames(IACT_rand_A) <- ES_and_var
IACT_rand_A <- cbind(IACT_rand_A, stimulus_type_.M.A. = "A", seed_ID = "IntACT")
IACT_rand_A

# "GrSoc" ; A : 8, M : 22 
GRS_A <- subset(GRS, subset = stimulus_type_.M.A. == "A" )
GRS_M <- subset(GRS, subset = stimulus_type_.M.A. == "M")

GRS_stim_samp1_ES <- data.frame(y = 1:nrow(GRS_A))
GRS_stim_samp1_VAR <- data.frame(y = 1:nrow(GRS_A))

for (i in randomizing) {
  samp1_GRS <- GRS_M[sample(1:nrow(GRS_M), nrow(GRS_A), replace=FALSE), ES_and_var]
  samp1_GRS_ES <- samp1_GRS[,"Effect_size"]
  samp1_GRS_VAR <- samp1_GRS[,"ES_var"]
  GRS_stim_samp1_ES <- cbind(GRS_stim_samp1_ES, samp1_GRS_ES, deparse.level = 1)
  GRS_stim_samp1_VAR <- cbind(GRS_stim_samp1_VAR, samp1_GRS_VAR, deparse.level = 1)
}

GRS_stim_samp1_ES$"Rand_ES" <- rowMeans(GRS_stim_samp1_ES)
GRS_stim_samp1_VAR$"Rand_VAR" <- rowMeans(GRS_stim_samp1_VAR)
GRS_rand_M <- cbind(GRS_stim_samp1_ES$Rand_ES, GRS_stim_samp1_VAR$Rand_VAR)
colnames(GRS_rand_M) <- ES_and_var
GRS_rand_M <- cbind(GRS_rand_M, stimulus_type_.M.A. = "M", seed_ID = "GrSoc")

# Building Stim df with randomized studies
SEL_VEC <- c("seed_ID", "stimulus_type_.M.A.", "Effect_size", "ES_var")
rand_stim <- rbind(FBl_A[SEL_VEC], FBl_rand_M, NDi_A[SEL_VEC], NDi_rand_M, HTS_A[SEL_VEC], HTS_rand_M, sIS_A[SEL_VEC], sIS_rand_M, pIS_M[SEL_VEC], pIS_rand_A, IACT_M[SEL_VEC], IACT_rand_A, GRS_A[SEL_VEC], GRS_rand_M)
rand_stim
rand_stim$Effect_size <- as.numeric(rand_stim$Effect_size)
rand_stim$ES_var <- as.numeric(rand_stim$ES_var)
summary(rand_stim$seed_ID)
# Subgroup analysis for stimulus on this smaller randomized data.frame

mm_rand_stim_m <- rma(yi = rand_stim$Effect_size, vi = rand_stim$ES_var, data=rand_stim, subset=stimulus_type_.M.A.== "M")
mm_rand_stim_a <- rma(yi = rand_stim$Effect_size, vi = rand_stim$ES_var, data = rand_stim, subset=stimulus_type_.M.A.== "A")
rand_stim_stim <- data.frame(estimate = c(coef(mm_rand_stim_m), coef(mm_rand_stim_a)), stderror = c(mm_rand_stim_m$se, mm_rand_stim_a$se),
                       meta = c("manual","automatic"), tau2 = round(c(mm_rand_stim_m$tau2, mm_rand_stim_a$tau2),3))

mm_rand_stim_m
mm_rand_stim_a

rand_stim_stim_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = rand_stim_stim, digits =3)
rand_stim_stim_comp


### Coding type  #### 

# FBl  ## not L : 5, Live : 6

FBl_notlive <- subset(FBl, subset = coding_type_.A.L.H.R. != "L" )
FBl_live <- subset(FBl, subset = coding_type_.A.L.H.R. == "L")
nrow(FBl) == nrow(FBl_notlive) + nrow(FBl_live)


FBl_cod_samp1_ES <- data.frame(y = 1:nrow(FBl_notlive))
FBl_cod_samp1_VAR <- data.frame(y = 1:nrow(FBl_notlive))

for (i in randomizing) {
  samp1_FBl <- FBl_live[sample(1:nrow(FBl_live), nrow(FBl_notlive), replace=FALSE), ES_and_var]
  samp1_FBl_ES <- as.numeric(samp1_FBl[,"Effect_size"])
  samp1_FBl_VAR <- as.numeric(samp1_FBl[,"ES_var"])
  FBl_cod_samp1_ES <- cbind(FBl_cod_samp1_ES, samp1_FBl_ES, deparse.level = 1)
  FBl_cod_samp1_VAR <- cbind(FBl_cod_samp1_VAR, samp1_FBl_VAR, deparse.level = 1)
  
}
FBl_cod_samp1_ES$"Rand_ES" <- rowMeans(FBl_cod_samp1_ES)
FBl_cod_samp1_VAR$"Rand_VAR" <- rowMeans(FBl_cod_samp1_VAR)
FBl_rand_live <- cbind(FBl_cod_samp1_ES$Rand_ES, FBl_cod_samp1_VAR$Rand_VAR)
colnames(FBl_rand_live) <- ES_and_var
FBl_rand_live <- cbind(FBl_rand_live, coding_type_.A.L.H.R. = "L", seed_ID = "FBl")
FBl_notlive$coding_type_.A.L.H.R. <- "Not L"


# "NDis" ; Not Live : 5, Live : 3
NDi_notlive <- subset(NDi, subset = coding_type_.A.L.H.R. != "L" )
NDi_live <- subset(NDi, subset = coding_type_.A.L.H.R. == "L")
nrow(NDi) == nrow(NDi_notlive) + nrow(NDi_live)

NDi_cod_samp1_ES <- data.frame(y = 1:nrow(NDi_live))
NDi_cod_samp1_VAR <- data.frame(y = 1:nrow(NDi_live))

for (i in randomizing) {
  samp1_NDi <- NDi_notlive[sample(1:nrow(NDi_notlive), nrow(NDi_live), replace=FALSE), ES_and_var]
  samp1_NDi_ES <- samp1_NDi[,"Effect_size"]
  samp1_NDi_VAR <- samp1_NDi[,"ES_var"]
  NDi_cod_samp1_ES <- cbind(NDi_cod_samp1_ES, samp1_NDi_ES, deparse.level = 1)
  NDi_cod_samp1_VAR <- cbind(NDi_cod_samp1_VAR, samp1_NDi_VAR, deparse.level = 1)
  
}

NDi_cod_samp1_ES$"Rand_ES" <- rowMeans(NDi_cod_samp1_ES)
NDi_cod_samp1_VAR$"Rand_VAR" <- rowMeans(NDi_cod_samp1_VAR)
NDi_rand_notlive <- cbind(NDi_cod_samp1_ES$Rand_ES, NDi_cod_samp1_VAR$Rand_VAR)
colnames(NDi_rand_notlive) <- ES_and_var
NDi_rand_notlive <- cbind(NDi_rand_notlive, coding_type_.A.L.H.R. = "Not L", seed_ID = "NDis")
NDi_live$coding_type_.A.L.H.R. <- "L"

# "wordWSeg"  ## Not Live : 9 ; Live : 38
wWS_notlive <- subset(wWS, subset = coding_type_.A.L.H.R. != "L" )
wWS_live <- subset(wWS, subset = coding_type_.A.L.H.R. == "L")
nrow(wWS) == nrow(wWS_notlive) + nrow(wWS_live)

wWS_cod_samp1_ES <- data.frame(y = 1:nrow(wWS_notlive))
wWS_cod_samp1_VAR <- data.frame(y = 1:nrow(wWS_notlive))

for (i in randomizing) {
  samp1_wWS <- wWS_live[sample(1:nrow(wWS_live), nrow(wWS_notlive), replace=FALSE), ES_and_var]
  samp1_wWS_ES <- samp1_wWS[,"Effect_size"]
  samp1_wWS_VAR <- samp1_wWS[,"ES_var"]
  wWS_cod_samp1_ES <- cbind(wWS_cod_samp1_ES, samp1_wWS_ES, deparse.level = 1)
  wWS_cod_samp1_VAR <- cbind(wWS_cod_samp1_VAR, samp1_wWS_VAR, deparse.level = 1)
  
}
wWS_cod_samp1_ES$"Rand_ES" <- rowMeans(wWS_cod_samp1_ES)
wWS_cod_samp1_VAR$"Rand_VAR" <- rowMeans(wWS_cod_samp1_VAR)
wWS_rand_live <- cbind(wWS_cod_samp1_ES$Rand_ES, wWS_cod_samp1_VAR$Rand_VAR)
colnames(wWS_rand_live) <- ES_and_var
wWS_rand_live <- cbind(wWS_rand_live, coding_type_.A.L.H.R. = "L", seed_ID = "wordWSeg")
wWS_notlive$coding_type_.A.L.H.R. <- "Not L"

#"speechWSeg" ; Not Live : 11 ; Live : 25
sWS_notlive <- subset(sWS, subset = coding_type_.A.L.H.R. != "L" )
sWS_live <- subset(sWS, subset = coding_type_.A.L.H.R. == "L")
nrow(sWS) == nrow(sWS_notlive) + nrow(sWS_live)

sWS_cod_samp1_ES <- data.frame(y = 1:nrow(sWS_notlive))
sWS_cod_samp1_VAR <- data.frame(y = 1:nrow(sWS_notlive))

for (i in randomizing) {
  samp1_sWS <- sWS_live[sample(1:nrow(sWS_live), nrow(sWS_notlive), replace=FALSE), ES_and_var]
  samp1_sWS_ES <- samp1_sWS[,"Effect_size"]
  samp1_sWS_VAR <- samp1_sWS[,"ES_var"]
  sWS_cod_samp1_ES <- cbind(sWS_cod_samp1_ES, samp1_sWS_ES, deparse.level = 1)
  sWS_cod_samp1_VAR <- cbind(sWS_cod_samp1_VAR, samp1_sWS_VAR, deparse.level = 1)
  
}
sWS_cod_samp1_ES$"Rand_ES" <- rowMeans(sWS_cod_samp1_ES)
sWS_cod_samp1_VAR$"Rand_VAR" <- rowMeans(sWS_cod_samp1_VAR)
sWS_rand_live <- cbind(sWS_cod_samp1_ES$Rand_ES, sWS_cod_samp1_VAR$Rand_VAR)
colnames(sWS_rand_live) <- ES_and_var
sWS_rand_live <- cbind(sWS_rand_live, coding_type_.A.L.H.R. = "L", seed_ID = "speechWSeg")
sWS_notlive$coding_type_.A.L.H.R. <- "Not L"

# "HTScram" ; Not live : 11 ; Live : 2
HTS_notlive <- subset(HTS, subset = coding_type_.A.L.H.R. != "L" )
HTS_live <- subset(HTS, subset = coding_type_.A.L.H.R. == "L")
nrow(HTS) == nrow(HTS_notlive) + nrow(HTS_live)

HTS_cod_samp1_ES <- data.frame(y = 1:nrow(HTS_live))
HTS_cod_samp1_VAR <- data.frame(y = 1:nrow(HTS_live))

for (i in randomizing) {
  samp1_HTS <- HTS_notlive[sample(1:nrow(HTS_notlive), nrow(HTS_live), replace=FALSE), ES_and_var]
  samp1_HTS_ES <- samp1_HTS[,"Effect_size"]
  samp1_HTS_VAR <- samp1_HTS[,"ES_var"]
  HTS_cod_samp1_ES <- cbind(HTS_cod_samp1_ES, samp1_HTS_ES, deparse.level = 1)
  HTS_cod_samp1_VAR <- cbind(HTS_cod_samp1_VAR, samp1_HTS_VAR, deparse.level = 1)
  
}

HTS_cod_samp1_ES$"Rand_ES" <- rowMeans(HTS_cod_samp1_ES)
HTS_cod_samp1_VAR$"Rand_VAR" <- rowMeans(HTS_cod_samp1_VAR)
HTS_rand_notlive <- cbind(HTS_cod_samp1_ES$Rand_ES, HTS_cod_samp1_VAR$Rand_VAR)
colnames(HTS_rand_notlive) <- ES_and_var
HTS_rand_notlive <- cbind(HTS_rand_notlive, coding_type_.A.L.H.R. = "Not L", seed_ID = "HTScram")
HTS_live$coding_type_.A.L.H.R. <- "L"


# "samInfS" ; Not Live : 2 ; Live : 1
sIS_notlive <- subset(sIS, subset = coding_type_.A.L.H.R. != "L" )
sIS_live <- subset(sIS, subset = coding_type_.A.L.H.R. == "L")
nrow(sIS) == nrow(sIS_notlive) + nrow(sIS_live)

sIS_cod_samp1_ES <- data.frame(y = 1:nrow(sIS_live))
sIS_cod_samp1_VAR <- data.frame(y = 1:nrow(sIS_live))

for (i in randomizing) {
  samp1_sIS <- sIS_notlive[sample(1:nrow(sIS_notlive), nrow(sIS_live), replace=FALSE), ES_and_var]
  samp1_sIS_ES <- samp1_sIS[,"Effect_size"]
  samp1_sIS_VAR <- samp1_sIS[,"ES_var"]
  sIS_cod_samp1_ES <- cbind(sIS_cod_samp1_ES, samp1_sIS_ES, deparse.level = 1)
  sIS_cod_samp1_VAR <- cbind(sIS_cod_samp1_VAR, samp1_sIS_VAR, deparse.level = 1)
  
}

sIS_cod_samp1_ES$"Rand_ES" <- rowMeans(sIS_cod_samp1_ES)
sIS_cod_samp1_VAR$"Rand_VAR" <- rowMeans(sIS_cod_samp1_VAR)
sIS_rand_notlive <- cbind(sIS_cod_samp1_ES$Rand_ES, sIS_cod_samp1_VAR$Rand_VAR)
colnames(sIS_rand_notlive) <- ES_and_var
sIS_rand_notlive <- cbind(sIS_rand_notlive, coding_type_.A.L.H.R. = "Not L", seed_ID = "samInfS")
sIS_live$coding_type_.A.L.H.R. <- "L"


# "popInfS" ; Not Live : 2 ; Live : 8
pIS_notlive <- subset(pIS, subset = coding_type_.A.L.H.R. != "L" )
pIS_live <- subset(pIS, subset = coding_type_.A.L.H.R. == "L")
nrow(pIS) == nrow(pIS_notlive)+ nrow(pIS_live)

pIS_cod_samp1_ES <- data.frame(y = 1:nrow(pIS_notlive))
pIS_cod_samp1_VAR <- data.frame(y = 1:nrow(pIS_notlive))

for (i in randomizing) {
  samp1_pIS <- pIS_live[sample(1:nrow(pIS_live), nrow(pIS_notlive), replace=FALSE), ES_and_var]
  samp1_pIS_ES <- samp1_pIS[,"Effect_size"]
  samp1_pIS_VAR <- samp1_pIS[,"ES_var"]
  pIS_cod_samp1_ES <- cbind(pIS_cod_samp1_ES, samp1_pIS_ES, deparse.level = 1)
  pIS_cod_samp1_VAR <- cbind(pIS_cod_samp1_VAR, samp1_pIS_VAR, deparse.level = 1)
  
}
pIS_cod_samp1_ES$"Rand_ES" <- rowMeans(pIS_cod_samp1_ES)
pIS_cod_samp1_VAR$"Rand_VAR" <- rowMeans(pIS_cod_samp1_VAR)
pIS_rand_live <- cbind(pIS_cod_samp1_ES$Rand_ES, pIS_cod_samp1_VAR$Rand_VAR)
colnames(pIS_rand_live) <- ES_and_var
pIS_rand_live <- cbind(pIS_rand_live, coding_type_.A.L.H.R. = "L", seed_ID = "popInfS")
pIS_notlive$coding_type_.A.L.H.R. <- "Not L"

# "VisStats" ; Not Live : 4 ; Live : 16
VST_notlive <- subset(VST, subset = coding_type_.A.L.H.R. != "L" )
VST_live <- subset(VST, subset = coding_type_.A.L.H.R. == "L")
nrow(VST) == nrow(VST_notlive) +nrow(VST_live)

VST_cod_samp1_ES <- data.frame(y = 1:nrow(VST_notlive))
VST_cod_samp1_VAR <- data.frame(y = 1:nrow(VST_notlive))

for (i in randomizing) {
  samp1_VST <- VST_live[sample(1:nrow(VST_live), nrow(VST_notlive), replace=FALSE), ES_and_var]
  samp1_VST_ES <- samp1_VST[,"Effect_size"]
  samp1_VST_VAR <- samp1_VST[,"ES_var"]
  VST_cod_samp1_ES <- cbind(VST_cod_samp1_ES, samp1_VST_ES, deparse.level = 1)
  VST_cod_samp1_VAR <- cbind(VST_cod_samp1_VAR, samp1_VST_VAR, deparse.level = 1)
  
}
VST_cod_samp1_ES$"Rand_ES" <- rowMeans(VST_cod_samp1_ES)
VST_cod_samp1_VAR$"Rand_VAR" <- rowMeans(VST_cod_samp1_VAR)
VST_rand_live <- cbind(VST_cod_samp1_ES$Rand_ES, VST_cod_samp1_VAR$Rand_VAR)
colnames(VST_rand_live) <- ES_and_var
VST_rand_live <- cbind(VST_rand_live, coding_type_.A.L.H.R. = "L", seed_ID = "VisStats")
VST_notlive$coding_type_.A.L.H.R. <- "Not L"



# "IntAct" ; Not Live : 2 ; Live : 5
IACT_notlive <- subset(IACT, subset = coding_type_.A.L.H.R. != "L" )
IACT_live <- subset(IACT, subset = coding_type_.A.L.H.R. == "L")
nrow(IACT) == nrow(IACT_notlive) + nrow(IACT_live)

IACT_cod_samp1_ES <- data.frame(y = 1:nrow(IACT_notlive))
IACT_cod_samp1_VAR <- data.frame(y = 1:nrow(IACT_notlive))

for (i in randomizing) {
  samp1_IACT <- IACT_live[sample(1:nrow(IACT_live), nrow(IACT_notlive), replace=FALSE), ES_and_var]
  samp1_IACT_ES <- samp1_IACT[,"Effect_size"]
  samp1_IACT_VAR <- samp1_IACT[,"ES_var"]
  IACT_cod_samp1_ES <- cbind(IACT_cod_samp1_ES, samp1_IACT_ES, deparse.level = 1)
  IACT_cod_samp1_VAR <- cbind(IACT_cod_samp1_VAR, samp1_IACT_VAR, deparse.level = 1)
  
}
IACT_cod_samp1_ES$"Rand_ES" <- rowMeans(IACT_cod_samp1_ES)
IACT_cod_samp1_VAR$"Rand_VAR" <- rowMeans(IACT_cod_samp1_VAR)
IACT_rand_live <- cbind(IACT_cod_samp1_ES$Rand_ES, IACT_cod_samp1_VAR$Rand_VAR)
colnames(IACT_rand_live) <- ES_and_var
IACT_rand_live <- cbind(IACT_rand_live, coding_type_.A.L.H.R. = "L", seed_ID = "IntAct")
IACT_notlive$coding_type_.A.L.H.R. <- "Not L"

# "GrSoc" ; Not Live : 25 ; Live : 5
GRS_notlive <- subset(GRS, subset = coding_type_.A.L.H.R. != "L" )
GRS_live <- subset(GRS, subset = coding_type_.A.L.H.R. == "L")


GRS_cod_samp1_ES <- data.frame(y = 1:nrow(GRS_live))
GRS_cod_samp1_VAR <- data.frame(y = 1:nrow(GRS_live))

for (i in randomizing) {
  samp1_GRS <- GRS_notlive[sample(1:nrow(GRS_notlive), nrow(GRS_live), replace=FALSE), ES_and_var]
  samp1_GRS_ES <- samp1_GRS[,"Effect_size"]
  samp1_GRS_VAR <- samp1_GRS[,"ES_var"]
  GRS_cod_samp1_ES <- cbind(GRS_cod_samp1_ES, samp1_GRS_ES, deparse.level = 1)
  GRS_cod_samp1_VAR <- cbind(GRS_cod_samp1_VAR, samp1_GRS_VAR, deparse.level = 1)
  
}

GRS_cod_samp1_ES$"Rand_ES" <- rowMeans(GRS_cod_samp1_ES)
GRS_cod_samp1_VAR$"Rand_VAR" <- rowMeans(GRS_cod_samp1_VAR)
GRS_rand_notlive <- cbind(GRS_cod_samp1_ES$Rand_ES, GRS_cod_samp1_VAR$Rand_VAR)
colnames(GRS_rand_notlive) <- ES_and_var
GRS_rand_notlive <- cbind(GRS_rand_notlive, coding_type_.A.L.H.R. = "Not L", seed_ID = "GrSoc")
GRS_live$coding_type_.A.L.H.R. <- "L"

# Building Coding df with randomized studies

SEL_VEC <- c("seed_ID", "coding_type_.A.L.H.R.", "Effect_size", "ES_var")
rand_cod <- rbind(FBl_rand_live, FBl_notlive[SEL_VEC], NDi_rand_notlive, NDi_live[SEL_VEC], wWS_rand_live, wWS_notlive[SEL_VEC], sWS_rand_live, sWS_notlive[SEL_VEC], HTS_rand_notlive, HTS_live[SEL_VEC], sIS_rand_notlive, sIS_live[SEL_VEC], pIS_rand_live, pIS_notlive[SEL_VEC], VST_rand_live, VST_notlive[SEL_VEC], IACT_rand_live, IACT_notlive[SEL_VEC], GRS_rand_notlive, GRS_live[SEL_VEC], stringsAsFactors = FALSE)
rand_cod$Effect_size <- as.numeric(rand_cod$Effect_size)
rand_cod$ES_var <- as.numeric(rand_cod$ES_var)
rand_cod$coding_type_.A.L.H.R. <- factor(rand_cod$coding_type_.A.L.H.R.)
rand_cod$seed_ID <- factor(rand_cod$seed_ID)


# Subgroup analysis on randomized sample for Coding 
mm_rand_cod_l <- rma(yi = rand_cod$Effect_size, vi = rand_cod$ES_var, data = rand_cod, subset=coding_type_.A.L.H.R. == "L")
mm_rand_cod_notl <- rma(yi = rand_cod$Effect_size, vi = rand_cod$ES_var, data = rand_cod, subset=coding_type_.A.L.H.R. != "L")
rand_cod_cod <- data.frame(estimate = c(coef(mm_rand_cod_l), coef(mm_rand_cod_notl)), stderror = c(mm_rand_cod_l$se, mm_rand_cod_notl$se),
                      meta = c("live","a_notlive"), tau2 = round(c(mm_rand_cod_l$tau2, mm_rand_cod_notl$tau2),3))

mm_rand_cod_l
mm_rand_cod_notl

rand_cod_cod_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = rand_cod_cod, digits =3)
rand_cod_cod_comp



### Inclusion level ####
#"FBl"  # Low : 5 ; High : 6
FBlex <- FBl[complete.cases(FBl$p_ratio), ]
FBlex_low <- subset(FBlex, subset=p_ratio < median(MA1ex$p_ratio) )
FBlex_high <- subset(FBlex, subset=p_ratio >= median(MA1ex$p_ratio))


FBlex_samp1_ES <- data.frame(y = 1:nrow(FBlex_low))
FBlex_samp1_VAR <- data.frame(y = 1:nrow(FBlex_low))

for (i in randomizing) {
  samp1ex_FBl <- FBlex_high[sample(1:nrow(FBlex_high), nrow(FBlex_low), replace=FALSE), ES_and_var]
  samp1ex_FBl_ES <- as.numeric(samp1ex_FBl[,"Effect_size"])
  samp1ex_FBl_VAR <- as.numeric(samp1ex_FBl[,"ES_var"])
  FBlex_samp1_ES <- cbind(FBlex_samp1_ES, samp1ex_FBl_ES, deparse.level = 1)
  FBlex_samp1_VAR <- cbind(FBlex_samp1_VAR, samp1ex_FBl_VAR, deparse.level = 1)
  
}
FBlex_samp1_ES$"Rand_ES" <- rowMeans(FBlex_samp1_ES)
FBlex_samp1_VAR$"Rand_VAR" <- rowMeans(FBlex_samp1_VAR)
FBlex_rand_high <- cbind(FBlex_samp1_ES$Rand_ES, FBlex_samp1_VAR$Rand_VAR)
colnames(FBlex_rand_high) <- ES_and_var
FBlex_rand_high <- cbind(FBlex_rand_high, p_ratio = "high", seed_ID = "FBl")
FBlex_low$p_ratio <- "low"

#"NDis"  # Low : 4, High 4, so we simply take the full data set
NDiex <- NDi[complete.cases(NDi$p_ratio), ]
NDiex_low <- subset(NDiex, subset=p_ratio < median(MA1ex$p_ratio) )
NDiex_high <- subset(NDiex, subset=p_ratio >= median(MA1ex$p_ratio))
NDiex$p_ratio <- ifelse(NDiex$p_ratio >= median(MA1ex$p_ratio), "high", "low")


#"wordWSeg" # Low : 29, High : 13
wWSex <- wWS[complete.cases(wWS$p_ratio), ]
wWSex_low <- subset(wWSex, subset=p_ratio < median(MA1ex$p_ratio) )
wWSex_high <- subset(wWSex, subset=p_ratio >= median(MA1ex$p_ratio))


wWSex_samp1_ES <- data.frame(y = 1:nrow(wWSex_high))
wWSex_samp1_VAR <- data.frame(y = 1:nrow(wWSex_high))

for (i in randomizing) {
  samp1ex_wWS <- wWSex_low[sample(1:nrow(wWSex_low), nrow(wWSex_high), replace=FALSE), ES_and_var]
  samp1ex_wWS_ES <- as.numeric(samp1ex_wWS[,"Effect_size"])
  samp1ex_wWS_VAR <- as.numeric(samp1ex_wWS[,"ES_var"])
  wWSex_samp1_ES <- cbind(wWSex_samp1_ES, samp1ex_wWS_ES, deparse.level = 1)
  wWSex_samp1_VAR <- cbind(wWSex_samp1_VAR, samp1ex_wWS_VAR, deparse.level = 1)
  
}
wWSex_samp1_ES$"Rand_ES" <- rowMeans(wWSex_samp1_ES)
wWSex_samp1_VAR$"Rand_VAR" <- rowMeans(wWSex_samp1_VAR)
wWSex_rand_low <- cbind(wWSex_samp1_ES$Rand_ES, wWSex_samp1_VAR$Rand_VAR)
colnames(wWSex_rand_low) <- ES_and_var
wWSex_rand_low <- cbind(wWSex_rand_low, p_ratio = "low", seed_ID = "wordWSeg")
wWSex_high$p_ratio <- "high"


#"speechWSeg" ; Low : 13, High : 23
sWSex <- sWS[complete.cases(sWS$p_ratio), ]
sWSex_low <- subset(sWSex, subset=p_ratio < median(MA1ex$p_ratio) )
sWSex_high <- subset(sWSex, subset=p_ratio >= median(MA1ex$p_ratio))


sWSex_samp1_ES <- data.frame(y = 1:nrow(sWSex_low))
sWSex_samp1_VAR <- data.frame(y = 1:nrow(sWSex_low))

for (i in randomizing) {
  samp1ex_sWS <- sWSex_high[sample(1:nrow(sWSex_high), nrow(sWSex_low), replace=FALSE), ES_and_var]
  samp1ex_sWS_ES <- as.numeric(samp1ex_sWS[,"Effect_size"])
  samp1ex_sWS_VAR <- as.numeric(samp1ex_sWS[,"ES_var"])
  sWSex_samp1_ES <- cbind(sWSex_samp1_ES, samp1ex_sWS_ES, deparse.level = 1)
  sWSex_samp1_VAR <- cbind(sWSex_samp1_VAR, samp1ex_sWS_VAR, deparse.level = 1)
  
}
sWSex_samp1_ES$"Rand_ES" <- rowMeans(sWSex_samp1_ES)
sWSex_samp1_VAR$"Rand_VAR" <- rowMeans(sWSex_samp1_VAR)
sWSex_rand_high <- cbind(sWSex_samp1_ES$Rand_ES, sWSex_samp1_VAR$Rand_VAR)
colnames(sWSex_rand_high) <- ES_and_var
sWSex_rand_high <- cbind(sWSex_rand_high, p_ratio = "high", seed_ID = "speechWSeg")
sWSex_low$p_ratio <- "low"


#"HTScram", Low : 9, High : 3 
HTSex <- HTS[complete.cases(HTS$p_ratio), ]
HTSex_low <- subset(HTSex, subset=p_ratio < median(MA1ex$p_ratio) )
HTSex_high <- subset(HTSex, subset=p_ratio >= median(MA1ex$p_ratio))

HTSex_samp1_ES <- data.frame(y = 1:nrow(HTSex_high))
HTSex_samp1_VAR <- data.frame(y = 1:nrow(HTSex_high))

for (i in randomizing) {
  samp1ex_HTS <- HTSex_low[sample(1:nrow(HTSex_low), nrow(HTSex_high), replace=FALSE), ES_and_var]
  samp1ex_HTS_ES <- as.numeric(samp1ex_HTS[,"Effect_size"])
  samp1ex_HTS_VAR <- as.numeric(samp1ex_HTS[,"ES_var"])
  HTSex_samp1_ES <- cbind(HTSex_samp1_ES, samp1ex_HTS_ES, deparse.level = 1)
  HTSex_samp1_VAR <- cbind(HTSex_samp1_VAR, samp1ex_HTS_VAR, deparse.level = 1)
  
}
HTSex_samp1_ES$"Rand_ES" <- rowMeans(HTSex_samp1_ES)
HTSex_samp1_VAR$"Rand_VAR" <- rowMeans(HTSex_samp1_VAR)
HTSex_rand_low <- cbind(HTSex_samp1_ES$Rand_ES, HTSex_samp1_VAR$Rand_VAR)
colnames(HTSex_rand_low) <- ES_and_var
HTSex_rand_low <- cbind(HTSex_rand_low, p_ratio = "low", seed_ID = "HTScram")
HTSex_high$p_ratio <- "high"

#"samInfS" ; Low : 0, High 1, so we take 0
sISex <- sIS[complete.cases(sIS$p_ratio), ]
sISex_low <- subset(sISex, subset=p_ratio < median(MA1ex$p_ratio) )
sISex_high <- subset(sISex, subset=p_ratio >= median(MA1ex$p_ratio))


#"popInfS" ; Low : 3, High : 5
pISex <- pIS[complete.cases(pIS$p_ratio), ]
pISex_low <- subset(pISex, subset=p_ratio < median(MA1ex$p_ratio))
pISex_high <- subset(pISex, subset=p_ratio >= median(MA1ex$p_ratio))

pISex_samp1_ES <- data.frame(y = 1:nrow(pISex_low))
pISex_samp1_VAR <- data.frame(y = 1:nrow(pISex_low))

for (i in randomizing) {
  samp1ex_pIS <- pISex_high[sample(1:nrow(pISex_high), nrow(pISex_low), replace=FALSE), ES_and_var]
  samp1ex_pIS_ES <- as.numeric(samp1ex_pIS[,"Effect_size"])
  samp1ex_pIS_VAR <- as.numeric(samp1ex_pIS[,"ES_var"])
  pISex_samp1_ES <- cbind(pISex_samp1_ES, samp1ex_pIS_ES, deparse.level = 1)
  pISex_samp1_VAR <- cbind(pISex_samp1_VAR, samp1ex_pIS_VAR, deparse.level = 1)
  
}
pISex_samp1_ES$"Rand_ES" <- rowMeans(pISex_samp1_ES)
pISex_samp1_VAR$"Rand_VAR" <- rowMeans(pISex_samp1_VAR)
pISex_rand_high <- cbind(pISex_samp1_ES$Rand_ES, pISex_samp1_VAR$Rand_VAR)
colnames(pISex_rand_high) <- ES_and_var
pISex_rand_high <- cbind(pISex_rand_high, p_ratio = "high", seed_ID = "popInfs")
pISex_low$p_ratio <- "low"


#"VisStats" ; Low 4, High 16
VSTex <- VST[complete.cases(VST$p_ratio), ]
VSTex_low <- subset(VSTex, subset=p_ratio < median(MA1ex$p_ratio) )
VSTex_high <- subset(VSTex, subset=p_ratio >= median(MA1ex$p_ratio))

VSTex_samp1_ES <- data.frame(y = 1:nrow(VSTex_low))
VSTex_samp1_VAR <- data.frame(y = 1:nrow(VSTex_low))

for (i in randomizing) {
  samp1ex_VST <- VSTex_high[sample(1:nrow(VSTex_high), nrow(VSTex_low), replace=FALSE), ES_and_var]
  samp1ex_VST_ES <- as.numeric(samp1ex_VST[,"Effect_size"])
  samp1ex_VST_VAR <- as.numeric(samp1ex_VST[,"ES_var"])
  VSTex_samp1_ES <- cbind(VSTex_samp1_ES, samp1ex_VST_ES, deparse.level = 1)
  VSTex_samp1_VAR <- cbind(VSTex_samp1_VAR, samp1ex_VST_VAR, deparse.level = 1)
  
}
VSTex_samp1_ES$"Rand_ES" <- rowMeans(VSTex_samp1_ES)
VSTex_samp1_VAR$"Rand_VAR" <- rowMeans(VSTex_samp1_VAR)
VSTex_rand_high <- cbind(VSTex_samp1_ES$Rand_ES, VSTex_samp1_VAR$Rand_VAR)
colnames(VSTex_rand_high) <- ES_and_var
VSTex_rand_high <- cbind(VSTex_rand_high, p_ratio = "high", seed_ID = "VisStats")
VSTex_low$p_ratio <- "low"

#"IntAct" ; Low 7, High : 0
IACTex <- IACT[complete.cases(IACT$p_ratio), ]
IACTex_low <- subset(IACTex, subset=p_ratio < median(MA1ex$p_ratio) )
IACTex_high <- subset(IACTex, subset=p_ratio >= median(MA1ex$p_ratio))


#"GrSoc" ; Low : 13, High : 17
GRSex <- GRS[complete.cases(GRS$p_ratio), ]
GRSex_low <- subset(GRSex, subset=p_ratio < median(MA1ex$p_ratio) )
GRSex_high <- subset(GRSex, subset=p_ratio >= median(MA1ex$p_ratio))

GRSex_samp1_ES <- data.frame(y = 1:nrow(GRSex_low))
GRSex_samp1_VAR <- data.frame(y = 1:nrow(GRSex_low))

for (i in randomizing) {
  samp1ex_GRS <- GRSex_high[sample(1:nrow(GRSex_high), nrow(GRSex_low), replace=FALSE), ES_and_var]
  samp1ex_GRS_ES <- as.numeric(samp1ex_GRS[,"Effect_size"])
  samp1ex_GRS_VAR <- as.numeric(samp1ex_GRS[,"ES_var"])
  GRSex_samp1_ES <- cbind(GRSex_samp1_ES, samp1ex_GRS_ES, deparse.level = 1)
  GRSex_samp1_VAR <- cbind(GRSex_samp1_VAR, samp1ex_GRS_VAR, deparse.level = 1)
  
}
GRSex_samp1_ES$"Rand_ES" <- rowMeans(GRSex_samp1_ES)
GRSex_samp1_VAR$"Rand_VAR" <- rowMeans(GRSex_samp1_VAR)
GRSex_rand_high <- cbind(GRSex_samp1_ES$Rand_ES, GRSex_samp1_VAR$Rand_VAR)
colnames(GRSex_rand_high) <- ES_and_var
GRSex_rand_high <- cbind(GRSex_rand_high, p_ratio = "high", seed_ID = "GrSoc")
GRSex_low$p_ratio <- "low"


# Building Inclusion ratio df with randomized studies

SEL_VEC <- c("seed_ID", "p_ratio", "Effect_size", "ES_var")
rand_excl <- rbind(FBlex_rand_high, FBlex_low[SEL_VEC], NDiex[SEL_VEC], wWSex_rand_low, wWSex_high[SEL_VEC], sWSex_rand_high, sWSex_low[SEL_VEC], HTSex_rand_low, HTSex_high[SEL_VEC], pISex_rand_high, pISex_low[SEL_VEC], VSTex_rand_high, VSTex_low[SEL_VEC], GRSex_rand_high, GRSex_low[SEL_VEC], stringsAsFactors = FALSE)
rand_excl$Effect_size <- as.numeric(rand_excl$Effect_size)
rand_excl$ES_var <- as.numeric(rand_excl$ES_var)
rand_excl$p_ratio <- factor(rand_excl$p_ratio)
summary(rand_excl)


# Subgroup analysis on randomized sample for Ratio inclusion
mm_randex_high <- rma(yi = rand_excl$Effect_size, vi = rand_excl$ES_var, data = rand_excl, subset=p_ratio == "high")
mm_randex_low <- rma(yi = rand_excl$Effect_size, vi = rand_excl$ES_var, data = rand_excl, subset=p_ratio == "low")
randex_rat <- data.frame(estimate = c(coef(mm_randex_high), coef(mm_randex_low)), stderror = c(mm_randex_high$se, mm_randex_low$se),
                        meta = c("high","a_low"), tau2 = round(c(mm_randex_high$tau2, mm_randex_low$tau2),3))

mm_randex_high
mm_randex_low

rand_rat_comp <- rma(estimate, sei=stderror, mods = ~ meta, method = "FE", data = randex_rat, digits =3)
rand_rat_comp



#### Are the methodological factors a source of error ? #### 
# To compare variability, we compare the Pooled standard deviation of both groups : if it is significantly higher, then we can infer that probably the higher variability is lead by the methodological factor
### ONE : using the full data-set
# If the resulting p-value of the test is less than some significance level (typically 0.05),
# then the obtained differences in sample variances are unlikely to have occurred based on random sampling from a population.

#### Manual vs Automatic # 
stim_var_test <- t.test(ES_var ~ stimulus_type_.M.A., data=MA1)
MA1_m <- subset(MA1, subset= MA1$stimulus_type_.M.A. == "M")
sd_var_m <- sd(MA1_m$ES_var)
MA1_a <- subset(MA1, subset = MA1$stimulus_type_.M.A. == "A")
sd_var_a <- sd(MA1_a$ES_var)
stim_var_test
sd_var_m
sd_var_a

### Coding : Live vs not Live 
MA1_cod_var <- MA1
MA1_cod_var$coding_type_.A.L.H.R. <- ifelse(MA1_cod_var$coding_type_.A.L.H.R. == "L", "L", "Not L")
var_live <- subset(MA1_cod_var, subset = MA1_cod_var$coding_type_.A.L.H.R. == "L")
sd_var_live <- sd(var_live$ES_var)
var_notl <- subset(MA1_cod_var, subset = MA1_cod_var$coding_type_.A.L.H.R. == "Not L")
sd_var_notl <- sd(var_notl$ES_var)
cod_var_test <- t.test(ES_var ~ coding_type_.A.L.H.R., data=MA1_cod_var)
cod_var_test
sd_var_live
sd_var_notl


#### Inclusion ratio : High vs Low 
MA1ex_var <- MA1ex
MA1ex_var$p_ratio <- ifelse(MA1ex_var$p_ratio >= median(MA1ex_var$p_ratio), "high", "low" )
MA1ex_var$p_ratio <- factor(MA1ex_var$p_ratio)
rat_var_test <- t.test(ES_var ~ p_ratio, data=MA1ex_var)
rat_var_test
var_high <- subset(MA1ex_var, subset = MA1ex_var$p_ratio =="high")
var_low <- subset(MA1ex_var, subset = MA1ex_var$p_ratio == "low")
sd_var_high <- sd(var_high$ES_var)
sd_var_low <- sd(var_low$ES_var)
sd_var_high
sd_var_low


#### TWO : using the equilibrated (randomized sampling) datasets from earlier
## Man vs Auto 
stim_rand_var_test <- t.test(ES_var ~  stimulus_type_.M.A., data=rand_stim)
sd_var_rand_a <- sd(subset(rand_stim, subset= rand_stim$stimulus_type_.M.A. == "A")$ES_var)
sd_var_rand_m <- sd(subset(rand_stim, subset= rand_stim$stimulus_type_.M.A. == "M")$ES_var)
stim_rand_var_test
sd_var_rand_a
sd_var_rand_m

### Coding : Live vs Not Live  
rand_cod_var <- rand_cod
rand_cod_var$coding_type_.A.L.H.R. <- ifelse(rand_cod_var$coding_type_.A.L.H.R. == "L", "L", "Not L")
cod_rand_var_test <- t.test(ES_var ~ coding_type_.A.L.H.R., data=rand_cod_var)
var_rand_live <- subset(rand_cod_var, subset = rand_cod_var$coding_type_.A.L.H.R. == "L")
var_rand_notl <- subset(rand_cod_var, subset = rand_cod_var$coding_type_.A.L.H.R. == "Not L")
sd_var_rand_live <- sd(var_rand_live$ES_var)
sd_var_rand_notl <- sd(var_rand_notl$ES_var)
cod_rand_var_test
sd_var_rand_live
sd_var_rand_notl

## Inclusion ratio : High vs Low 
excl_rand_var_test <- t.test(ES_var ~ p_ratio, data=rand_excl)
excl_rand_var_test
var_rand_high <- subset(rand_excl, subset = rand_excl$p_ratio == "high")
var_rand_low <- subset(rand_excl, subset= rand_excl$p_ratio == "low")
sd_var_rand_high <- sd(var_rand_high$ES_var)
sd_var_rand_low <- sd(var_rand_low$ES_var)
sd_var_rand_high
sd_var_rand_low

## Building forest plots to make it all prettier

# FBl 
par(mar=c(4,4,1,2))
FBl_forest <- forest.rma(mm_FBl, slab = paste(FBl$study_ID, FBl$year, FBl$Study_type_.OR.R., sep = ", "), alim = c(-1, 5), refline = c(0.2), xlim=c(-10, 10), at =c(-1, 0, 0.2, 0.5, 2, 3.5, 5), ilab.xpos=c(-5,-3.5), ilab=cbind(FBl$x_1, FBl$x_2), ylim = c(0, 14), xlab = "Hedges' g : Effect Size")

op <- par(cex=.80, font=4)
text(c(-5,-3.5), 12.5, c("LT Unexp. /", "Exp."))
text(-8.5,                12.5, "Study",     pos=4)
text(10,                  12.5, "Relative Risk [95% CI]", pos=2)

par(op)
FBl_forest <- title("False Belief Experiment")

# NDi
par(mar=c(4,4,1,2))
NDi_forest <- forest.rma(mm_NDi, slab = paste(NDi$study_ID, NDi$year, NDi$Study_type_.OR.R., sep = ", "), refline = 0.2, at = c(-1, 0, 0.2, 0.5, 2), alim = c(-1, 2.5), xlim=c(-10, 10), ilab.xpos=c(-5,-3.5), ilab=cbind(NDi$x_1, NDi$x_2), ylim = c(0, 11), steps = 5, xlab = "Hedges' g : Effect Size")

op <- par(cex=.80, font=4)
text(c(-5,-3.5), 9.5, c("LT New. /", "Old #."))
text(-8.5,                9.5, "Study",     pos=4)
text(10,                  9.5, "Relative Risk [95% CI]", pos=2)

par(op)
NDi_forest <- title("Number Discrimination Experiment")

# wWS
wWS <- separate(wWS, study_ID, into = c("pres_type", "study_label"), sep = "/", remove = FALSE, convert = FALSE)
wWS <- separate(wWS, study_label, into = c("num", "label"), sep = " - ", remove = TRUE, convert = FALSE)
wWS$label <- "wWS"
wWS <- unite(wWS, "study_label", num:label, sep = " - ", remove = TRUE)
par(mar=c(4,4,1,2))
plot.new()
wWS_forest <- forest.rma(mm_wWS, slab = paste(wWS$study_label, wWS$year, wWS$Study_type_.OR.R., sep = ", "), refline = 0.2, at= c(-2, -1, 0, 1, 2), alim=c(-2.5, 2), xlim=c(-10, 10), ilab.xpos=c(-5,-3.5), ilab=cbind(wWS$x_1, wWS$x_2), ylim = c(-1, 50), steps = 5, xlab = "Hedges' g : Effect Size")
op <- par(cex=.80, font=4)
text(c(-5.3,-3.5), 48.5, c("Fam. /", "Unfam. word"))
text(-8.5,                48.5, "Study",     pos=4)
text(10,                  48.5, "Relative Risk [95% CI]", pos=2)

par(op)
wWS_forest <- title("w. Word Segmentation Experiment")

# sWS
sWS <- separate(sWS, study_ID, into = c("pres_type", "study_label"), sep = "/", remove = FALSE, convert = FALSE)
sWS <- separate(sWS, study_label, into = c("num", "label"), sep = " - ", remove = TRUE, convert = FALSE)
sWS$label <- "sWS"
sWS <- unite(sWS, "study_label", num:label, sep = " - ", remove = TRUE)
par(mar=c(4,4,1,2))
plot.new()
sWS_forest <- forest.rma(mm_sWS, slab = paste(sWS$study_label, sWS$year, sWS$Study_type_.OR.R., sep = ", "), alim = c(-1.5, 1.5), refline = 0.2, xlim=c(-10, 10), ilab.xpos=c(-5,-3.5), ilab=cbind(sWS$x_1, sWS$x_2), ylim = c(-1, 39), steps = 5, cex = 0.95, xlab = "Hedges' g : Effect Size")
op <- par(cex=.80, font=4)
text(c(-5.3,-3.5), 37.5, c("Fam. /", "Unfam. word"))
text(-8.5,                37.5, "Study",     pos=4)
text(10,                  37.5, "Relative Risk [95% CI]", pos=2)

par(op)
sWS_forest <- title("sp. Word Segmentation Experiment")

# sIS
sIS <- separate(sIS, study_ID, into = c("pres_type", "study_label"), sep = "/", remove = FALSE, convert = FALSE)
sIS <- separate(sIS, study_label, into = c("num", "label"), sep = " - ", remove = TRUE, convert = FALSE)
sIS$label <- "sIS"
sIS <- unite(sIS, "study_label", num:label, sep = " - ", remove = TRUE)
par(mar=c(4,4,1,2))
plot.new()
sIS_forest <- forest.rma(mm_sIS, slab = paste(sIS$study_label, sIS$year, sIS$Study_type_.OR.R., sep = ", "), refline = 0.2, xlim=c(-10, 10), ilab.xpos=c(-5,-3.5), alim = c(-0.5, 1.30), ilab=cbind(sIS$x_1, sIS$x_2), ylim = c(-1, 6), steps = 5, xlab = "Hedges' g : Effect Size")
op <- par(cex=.80, font=4)
text(c(-5.3,-3.5), 4.5, c("Unexp. /", "Expected"))
text(-8.5,                4.5, "Study",     pos=4)
text(10,                  4.5, "Relative Risk [95% CI]", pos=2)

par(op)
sIS_forest <- title("s. Inferential Statistics Experiment")



# pIS
pIS <- separate(pIS, study_ID, into = c("pres_type", "study_label"), sep = "/", remove = FALSE, convert = FALSE)
pIS <- separate(pIS, study_label, into = c("num", "label"), sep = " - ", remove = TRUE, convert = FALSE)
pIS$label <- "pIS"
pIS <- unite(pIS, "study_label", num:label, sep = " - ", remove = TRUE)
par(mar=c(4,4,1,2))
plot.new()
pIS_forest <- forest.rma(mm_pIS, slab = paste(pIS$study_label, pIS$year, pIS$Study_type_.OR.R., sep = ", "), alim =c(-1, 1.30), refline = 0.2, xlim=c(-10, 10), ilab.xpos=c(-5, -3.5), ilab=cbind(pIS$x_1, pIS$x_2), ylim = c(-1, 13), cex = 1, xlab = "Hedges' g : Effect Size")
op <- par(cex=.80, font=4)
text(c(-5.3,-3.5), 11.5, c("Unexp. /", "Expected"))
text(-8.5,                11.5, "Study",     pos=4)
text(10,                  11.5, "Relative Risk [95% CI]", pos=2)

par(op)
pIS_forest <- title("pop Inferential Statistics Experiment")


# HTS
HTS <- separate(HTS, study_ID, into = c("num", "label"), sep = "-", remove = TRUE, convert = FALSE)
HTS$label <- "HTS"
HTS <- unite(HTS, "study_label", num:label, sep = " - ", remove = TRUE)

par(mar=c(4,4,1,2))
plot.new()
HTS_forest <- forest.rma(mm_HTS, slab = paste(HTS$study_label, HTS$year, HTS$Study_type_.OR.R., sep = ", "), xlim=c(-10, 10), alim = c(-0.5, 4), refline = 0.2, ilab.xpos=c(-5,-3.5), ilab=cbind(HTS$x_1, HTS$x_2), ylim = c(-1, 16), cex = 1, xlab = "Hedges' g : Effect Size")
op <- par(cex=.80, font=4)
text(c(-5.3,-3.5), 14.5, c("Scramb. /", "Normal Face"))
text(-8.5,                14.5, "Study",     pos=4)
text(10,                  14.5, "Relative Risk [95% CI]", pos=2)

par(op)
HTS_forest <- title("Head Turn to Faces Experiment")

# VST
VST <- separate(VST, study_ID, into = c("num", "label"), sep = "-", remove = TRUE, convert = FALSE)
VST$label <- "VST"
VST <- unite(VST, "study_label", num:label, sep = "-", remove = TRUE)
par(mar=c(4,4,1,2))
plot.new()
VST_forest <- forest.rma(mm_VST, slab = paste(VST$study_label, VST$year, VST$Study_type_.OR.R., sep = ", "), xlim=c(-13, 13), refline = 0.2, ilab.xpos=c(-7.5,-5.5), ilab=cbind(round(VST$x_1, 2), round(VST$x_2,2)), ylim = c(-1, 23), cex = 1, xlab = "Hedges' g : Effect Size")
op <- par(cex=.80, font=4)
text(c(-7.8,-5.3), 21.5, c("Novel /", "Fam Seq."))
text(-11.5,                21.5, "Study",     pos=4)
text(12,                  21.5, "Relative Risk [95% CI]", pos=2)

par(op)
VST_forest <- title("Visual Statistics Experiment")
# IACT
par(mar=c(4,4,1,2))
plot.new()
IACT_forest <- forest.rma(mm_IACT, slab = paste(IACT$study_ID, IACT$year, IACT$Study_type_.OR.R., sep = ", "), xlim=c(-10, 10), refline = 0.2, ilab.xpos=c(-4.7,-3.5), ilab=cbind(IACT$x_1, IACT$x_2), ylim = c(-1, 10), cex = 1, xlab = "Hedges' g : Effect Size")
op <- par(cex=.80, font=4)
text(c(-4.7, -3.4), 8.5, c("Old /", "New Action"))
text(-8.5,                8.5, "Study",     pos=4)
text(10,                  8.5, "Relative Risk [95% CI]", pos=2)

par(op)
IACT_forest <- title("Intentional Action Experiment")

# GRS
par(mar=c(4,4,1,2))
plot.new()
GRS_forest <- forest.rma(mm_GRS, slab = paste(GRS$study_ID, GRS$year, GRS$Study_type_.OR.R., sep = ", "), xlim=c(-10, 10), alim = c(0, 1.5), refline = 0.2, ilab.xpos=c(-4.7,-3.5), ilab=cbind(GRS$x_1, GRS$x_2), ylim = c(-1, 33), cex = 1, steps = 5, xlab = "Hedges' g : Effect Size")
op <- par(cex=.80, font=4)
text(c(-4.8, -3.4), 31.5, c("Social /", "Non Social"))
text(-8.5,                31.5, "Study",     pos=4)
text(10,                  31.5, "Relative Risk [95% CI]", pos=2)

par(op)
GRS_forest <- title("Grabing (More) Social Actor Experiment")


