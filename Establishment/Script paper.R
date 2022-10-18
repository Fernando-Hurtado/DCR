#Index####
  #0 - Data and libraries
        #0.1 - Data
        #0.2 - Libraries
  #1 - Formal analyses and pairwise t tests
        #1.1 - Number of established propagules (propagules/mm2)
        #1.2 - Colonized surface (%)
        #1.3 - Viable biomass (%)
        #1.4 - Relative Growth Rate (g/month)
        #1.5 - Table S.3 Wilcoxon tests
        #1.6 - Table S.4 Pearson's Chi-squared tests
  #2 - Main figure (Figure 1)
        #2.1 - General preparation
        #2.2 - Number of established propagules (propagules/mm2)
        #2.3 - Colonized surface (%)
        #2.4 - Viable biomass (%)
        #2.5 - Relative Growth Rate (g/month)
        #2.6 - Plots
        #2.7 - Mosaic plot integration and export of the plot
  #3 - Correlations of indicators with traits 1 (Figure 2)
        #3.1 - General preparation
        #3.2 - Correlations
        #3.3 - Plots
        #3.4 - Mosaic plot integration and export of the plot
  #4 - Correlations of indicators with traits 2 (Figure S1)
        #4.1 - General preparation
        #4.2 - Same parts of code of figure 2
        #4.3 - Different part of code: final Scatterplots
        #4.4 - Mosaic plot: Export

############ 0 - Data and libraries --------------------------------------------------------------------
# 0.1 - Data ####
  # Data for the formal analysis would be called in each section. 
  # The former data from the experiment consist in 216 samples of mosses of 6 moss species and 3 propagule sizes.
  # Also data from the propagule characterization is used to study the correlation of different measured traits with the establishment indicators.
  # Propagule characterization dataset consist in 540 samples of these 6 species

# 0.2 - Libraries ####  
  # Libraries for analyses, models and plotting would be called in each section, with a comment of the specific use in the script and the version of the package used.
  # It is also specified in the main article the general usage of the packages in the Methods section:
  # Link to the publication: https://assets.researchsquare.com/files/rs-1358759/v1/6d414dbb-521d-4b6e-bbb0-c29072bb5461.pdf?c=1647365290

############ 1 - Formal analyses and pairwise t tests ---------------------------------------------------------
#### 1.1 - Number of established propagules (propagules/mm2) -----------------------------------------------------------------------
## General preparation -----------------------------------------------
# Libraries ---------------------------------------------------------------

library(car) # for Levene's test, Anova and boxCox (v.3.1-0)
library(robustbase) # for robust models (v. 0.95-0)
library(stats) # for shapiro.test, ks.test, bartlett.test, fligner.test, aov for doing the models and preparing the transformation, residuals for getting the overdispersion, cor for the R2 (v. 4.1.2)
library(rstatix) # for the pairwise test and Cohen's d (v. 0.7.0)
library(dplyr) # for managing the data and the pipelines (v. 2.2.1)

# Dataset and configuration -----------------------------------------------------------------

path.to.data <- "https://raw.githubusercontent.com/Fernando-Hurtado/DCR/main/Establishment/"
dat <- read.csv(paste0(path.to.data,"Dataset_Indicators.csv"), header=T, sep=";", dec= ",")

options(contrasts=c(factor="contr.sum", ordered="contr.poly")) ## We have factors in the models, so we change the contrast matrix as it do not corresponds with the theory in stadistics, by deafult R applies a different one

# I create the interaction column for some of the tests for the models

attach(dat)
dat$interaccion <- Sp.:Size
detach(dat)


# Statistical analyses ----------------------------------------------------
# Creation of the models 

eqt <- as.formula(Shoots_established ~ Sp.)
eqt2 <- as.formula(Shoots_established ~ Sp.+Size)
eqt3 <- as.formula(Shoots_established ~ Sp.*Size)


Shoots_established_aov <- aov(eqt, data= dat)
Shoots_established_aov2 <- aov(eqt2, data= dat)
Shoots_established_aov3 <- aov(eqt3, data= dat)

# We take a look at the residuals and the qqnorm

plot(Shoots_established_aov)
plot(Shoots_established_aov2)
plot(Shoots_established_aov3)


# We use Kolmogorov- Smirnov tests, as it works good for a big N

ks.test(residuals(Shoots_established_aov), y= pnorm)
ks.test(residuals(Shoots_established_aov2), y= pnorm)
ks.test(residuals(Shoots_established_aov3), y= pnorm)

# We look now at the homocedasticity, we use bartlett (and Levene for the interaction) as it accepts some deviation of normality

bartlett.test(Shoots_established~Sp., data= dat)
bartlett.test(Shoots_established~Size, data= dat)
leveneTest(Shoots_established~Sp.:Size, data= dat, center=median)

# We have to transform the variable. From this part we actually tested many transformations and finally we kept with Box-cox transformation

Shoots_established <- dat$Shoots_established + min(dat$Shoots_established[dat$Shoots_established!=0], na.rm = TRUE) # Box-cox transformations do not deal well with 0 values, so we prepare the data for the transformation
dat$Shoots_established <- Shoots_established

lambda.bc <- with(boxCox(Shoots_established_aov, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
lambda.bc2 <- with(boxCox(Shoots_established_aov2, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
lambda.bc3 <- with(boxCox(Shoots_established_aov3, lambda=seq(-2,2, 1/1000)), x[which.max(y)])

Shoots_established_bc <- ((Shoots_established^lambda.bc)-1)/(lambda.bc)
Shoots_established_bc2 <- ((Shoots_established^lambda.bc2)-1)/(lambda.bc2)
Shoots_established_bc3 <- ((Shoots_established^lambda.bc3)-1)/(lambda.bc3)


# we incorporate the variables, corrected for each model, to the dataset

dat$Shoots_established_bc <- Shoots_established_bc
dat$Shoots_established_bc2 <- Shoots_established_bc2
dat$Shoots_established_bc3 <- Shoots_established_bc3


# Creation of the models with the transformed variable

eqt <- as.formula(Shoots_established_bc ~ Sp.)
eqt2 <- as.formula(Shoots_established_bc2 ~ Sp.+Size)
eqt3 <- as.formula(Shoots_established_bc3 ~ Sp.*Size)

Shoots_established_bc_aov <- aov(eqt, data= dat)
Shoots_established_bc_aov2 <- aov(eqt2, data= dat)
Shoots_established_bc_aov3 <- aov(eqt3, data= dat)

# Creation of the robust models 

modelo.lmrob <- lmrob(eqt, data=dat, method="SMDM", setting = "KS2014")
modelo.lmrob2 <- lmrob(eqt2, data=dat, method="SMDM", setting = "KS2014")
modelo.lmrob3 <- lmrob(eqt3, data=dat, method="SMDM", setting = "KS2014")

# Analyses of the models, first a look at the residuals and qqnorm

plot(Shoots_established_bc_aov)
plot(modelo.lmrob)
plot(Shoots_established_bc_aov2)
plot(modelo.lmrob2)
plot(Shoots_established_bc_aov3)
plot(modelo.lmrob3)

# Even it is difficult to do the interpretation, it is an improvement in general in normality and more compensated leverages in the robust models.

# We go with the tests 

ks.test(residuals(Shoots_established_bc_aov), y= pnorm)
ks.test(residuals(Shoots_established_bc_aov2), y= pnorm)
ks.test(residuals(Shoots_established_bc_aov3), y= pnorm)
ks.test(residuals(modelo.lmrob), y= pnorm)
ks.test(residuals(modelo.lmrob2), y= pnorm)
ks.test(residuals(modelo.lmrob3), y= pnorm)

# Bartlett and Levene will do for the transformed variables

bartlett.test(Shoots_established_bc~Sp., data= dat)
bartlett.test(Shoots_established_bc2~Sp., data= dat)
bartlett.test(Shoots_established_bc3~Sp., data= dat)
bartlett.test(Shoots_established_bc2~Size, data= dat)
bartlett.test(Shoots_established_bc3~Size, data= dat)
leveneTest(Shoots_established_bc3~Sp.:Size, data= dat, center=median)

# We pay attention now to the dispersion of the models, due to overdispersion inflates tipe I error

phi <- sum((residuals(Shoots_established_bc_aov, type="pearson"))^2)/Shoots_established_bc_aov$df.residual
print(c("Pearson overdispersion =", round(phi, 3)), quote=FALSE)
phirob <- sum((residuals(modelo.lmrob, type="pearson"))^2)/modelo.lmrob$df.residual
print(c("Pearson overdispersion =", round(phirob, 3)), quote=FALSE)
phi2 <- sum((residuals(Shoots_established_bc_aov2, type="pearson"))^2)/Shoots_established_bc_aov2$df.residual
print(c("Pearson overdispersion =", round(phi2, 3)), quote=FALSE)
phirob2 <- sum((residuals(modelo.lmrob2, type="pearson"))^2)/modelo.lmrob2$df.residual
print(c("Pearson overdispersion =", round(phirob2, 3)), quote=FALSE)
phi3 <- sum((residuals(Shoots_established_bc_aov3, type="pearson"))^2)/Shoots_established_bc_aov3$df.residual
print(c("Pearson overdispersion =", round(phi3, 3)), quote=FALSE)
phirob3 <- sum((residuals(modelo.lmrob3, type="pearson"))^2)/modelo.lmrob3$df.residual
print(c("Pearson overdispersion =", round(phirob3, 3)), quote=FALSE)

# The full models do it better
# Finally in this first part of the statistical analyses, we look to the F and p values, corrected by sandwich estimators (hc4)
# We prefer to keep going with the robust models, concerning the overall results so far

Anova(modelo.lmrob, type=3, test="F", white.adjust="hc4")
Anova(modelo.lmrob2, type=3, test="F", white.adjust="hc4")
Anova(modelo.lmrob3, type=3, test="F", white.adjust="hc4")


# Pairwise t tests --------------------------------------------------------

## We select the full model transformed variable. We look now at the pairwise comparisons
# Pairwise comparisons
pwc <- dat %>% 
  pairwise_t_test(
    Shoots_established_bc3 ~ Sp., pool.sd = TRUE,
    p.adjust.method = "holm"
  ) 

pwc2 <- dat %>% 
  pairwise_t_test(
    Shoots_established_bc3 ~ Size, pool.sd = TRUE,
    p.adjust.method = "holm"
  ) 

pwc3 <- dat %>% 
  pairwise_t_test(
    Shoots_established_bc3 ~ interaccion, pool.sd = TRUE,
    p.adjust.method = "holm"
  ) 

# We calculate the means by groups

mean <- dat %>%
  group_by(Sp.) %>%
  summarise(
    mean = mean(Shoots_established_bc3)
  )
mean2 <- dat %>%
  group_by(Size) %>%
  summarise(
    mean = mean(Shoots_established_bc3)
  )
mean3 <- dat %>%
  group_by(interaccion) %>%
  summarise(
    mean = mean(Shoots_established_bc3)
  )

# Differences between means

mean_dif <- as.data.frame(outer(mean$mean, mean$mean, FUN="-"))
colnames(mean_dif) <- rownames(mean_dif) <- mean$Sp.

mean_dif2 <- as.data.frame(outer(mean2$mean, mean2$mean, FUN="-"))
colnames(mean_dif2) <- rownames(mean_dif2) <- mean2$Size

mean_dif3 <- as.data.frame(outer(mean3$mean, mean3$mean, FUN="-"))
colnames(mean_dif3) <- rownames(mean_dif3) <- mean3$interaccion

# We order rows by columns to make later the output of the results to coincide in order with the pwc

to_long <- function(x){
  
  df <- x[order(rownames(x)),order(colnames(x))]  
  df <- data.frame( t(combn(rownames(df),2)), dist=t(df)[lower.tri(df)])
  colnames(df)[1:3] <- c("trat1", "trat2", "mean_dif")
  return(df)
}

mean_dif_long <- to_long(mean_dif)

mean_dif_long2 <- to_long(mean_dif2)

mean_dif_long3 <- to_long(mean_dif3)

# we include the differences in the pwc

pwc$mean_dif <-mean_dif_long[,"mean_dif"]
pwc2$mean_dif <-mean_dif_long2[,"mean_dif"]
pwc3$mean_dif <-mean_dif_long3[,"mean_dif"]

# We finally incorporated the Cohen's D for the effect size and magnitudes

cohen_dif <- dat %>% cohens_d(Shoots_established_bc3 ~ Sp.,var.equal = FALSE)
cohen_dif <- as.data.frame(cohen_dif)
cohen_dif2 <- dat %>% cohens_d(Shoots_established_bc3 ~ Size,var.equal = FALSE)
cohen_dif2 <- as.data.frame(cohen_dif2)
cohen_dif3 <- dat %>% cohens_d(Shoots_established_bc3 ~ interaccion,var.equal = FALSE)
cohen_dif3 <- as.data.frame(cohen_dif3)

pwc$effsize <- cohen_dif[,"effsize"]
pwc$magnitude <- cohen_dif[,"magnitude"]
pwc2$effsize <- cohen_dif2[,"effsize"]
pwc2$magnitude <- cohen_dif2[,"magnitude"]
pwc3$effsize <- cohen_dif3[,"effsize"]
pwc3$magnitude <- cohen_dif3[,"magnitude"]

# We end this part exporting the results. These results will be used for creating the Compact Letter Displays of the Fig.1

pwc <- as.data.frame(pwc)
pwc2 <- as.data.frame(pwc2)
pwc3 <- as.data.frame(pwc3)

write.csv(pwc, file= "NESpwc.csv", sep= ",")
write.csv(pwc2, file= "NESpwc2.csv", sep= ",")
write.csv(pwc3, file= "NESpwc3.csv", sep= ",")


# Cleaning for next part of the code ####

rm(list = ls()) # to clean data before start

#### 1.2 - Colonized surface (%) -----------------------------------------------------------------------
## General preparation -----------------------------------------------
# Libraries ---------------------------------------------------------------

library(car) # for Levene's test, Anova and boxCox (v.3.1-0)
library(robustbase) # for robust models (v. 0.95-0)
library(stats) # for shapiro.test, ks.test, bartlett.test, fligner.test, aov for doing the models and preparing the transformation, residuals for getting the overdispersion, cor for the R2 (v.4.1.2)
library(rstatix) # for the pairwise test and Cohen's d (v. 0.7.0)
library(dplyr) # for managing the data and the pipelines (v. 2.2.1)

# Dataset and configuration -----------------------------------------------------------------

path.to.data <- "https://raw.githubusercontent.com/Fernando-Hurtado/DCR/main/Establishment/"
dat <- read.csv(paste0(path.to.data,"Dataset_Indicators.csv"), header=T, sep=";", dec= ",")

options(contrasts=c(factor="contr.sum", ordered="contr.poly")) ## We have factors in the models, so we change the contrast matrix as it do not corresponds with the theory in stadistics, by deafult R applies a different one

# I create the interaction column for some of the tests for the models

attach(dat)
dat$interaccion <- Sp.:Size
detach(dat)

# We create the variable colonized surface

dat$perc_col <- dat$S_Shootss/dat$S_substratum
dat$perc_col <- dat$perc_col*100

## Statistical analyses ----------------------------------------------------
# Creation of the models 

eqt <- as.formula(perc_col ~ Sp.)
eqt2 <- as.formula(perc_col ~ Sp.+Size)
eqt3 <- as.formula(perc_col ~ Sp.*Size)

perc_col_aov <- aov(eqt, data= dat)
perc_col_aov2 <- aov(eqt2, data= dat)
perc_col_aov3 <- aov(eqt3, data= dat)

# We take a look at the residuals and the qqnorm

plot(perc_col_aov)
plot(perc_col_aov2)
plot(perc_col_aov3)

# We use Kolmogorov- Smirnov tests, as it works good for a big N

ks.test(residuals(perc_col_aov), y= pnorm)
ks.test(residuals(perc_col_aov2), y= pnorm)
ks.test(residuals(perc_col_aov3), y= pnorm)

# We look now at the homocedasticity, we use bartlett (and Levene for the interaction) as it accepts some deviation of normality

bartlett.test(perc_col~Sp., data= dat)
bartlett.test(perc_col~Size, data= dat)
leveneTest(perc_col~Sp.:Size, data= dat, center=median)

# We have to transform the variable. From this part we actually tested many transformations and finally we kept with Box-cox transformation

perc_col <- dat$perc_col + min(dat$perc_col[dat$perc_col!=0], na.rm = TRUE) # Box-cox transformations do not deal well with 0 values, so we prepare the data for the transformation
dat$perc_col <- perc_col


lambda.bc <- with(boxCox(perc_col_aov, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
lambda.bc2 <- with(boxCox(perc_col_aov2, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
lambda.bc3 <- with(boxCox(perc_col_aov3, lambda=seq(-2,2, 1/1000)), x[which.max(y)])

perc_col_bc <- ((perc_col^lambda.bc )-1)/(lambda.bc )
perc_col_bc2 <- ((perc_col^lambda.bc2)-1)/(lambda.bc2)
perc_col_bc3 <- ((perc_col^lambda.bc3)-1)/(lambda.bc3)


# we incorporate the variables, corrected for each model, to the dataset

dat$perc_col_bc <- perc_col_bc
dat$perc_col_bc2 <- perc_col_bc2
dat$perc_col_bc3 <- perc_col_bc3


# Creation of the models with the transformed variable

eqt <- as.formula(perc_col_bc ~ Sp.)
eqt2 <- as.formula(perc_col_bc2 ~ Sp.+Size)
eqt3 <- as.formula(perc_col_bc3 ~ Sp.*Size)


perc_col_bc_aov <- aov(eqt, data= dat)
perc_col_bc_aov2 <- aov(eqt2, data= dat)
perc_col_bc_aov3 <- aov(eqt3, data= dat)

# Creation of the robust models 

modelo.lmrob <- lmrob(eqt, data=dat, method="SMDM", setting = "KS2014")
modelo.lmrob2 <- lmrob(eqt2, data=dat, method="SMDM", setting = "KS2014")
modelo.lmrob3 <- lmrob(eqt3, data=dat, method="SMDM", setting = "KS2014")

# Analyses of the models, first a look at the residuals and qqnorm

plot(perc_col_bc_aov)
plot(modelo.lmrob)
plot(perc_col_bc_aov2)
plot(modelo.lmrob2)
plot(perc_col_bc_aov3)
plot(modelo.lmrob3)

# Even it is difficult to do the interpretation, it is an improvement in general in normality and more compensated leverages in the robust models.

# We go with the tests 

ks.test(residuals(perc_col_bc_aov), y= pnorm)
ks.test(residuals(perc_col_bc_aov2), y= pnorm)
ks.test(residuals(perc_col_bc_aov3), y= pnorm)
ks.test(residuals(modelo.lmrob), y= pnorm)
ks.test(residuals(modelo.lmrob2), y= pnorm)
ks.test(residuals(modelo.lmrob3), y= pnorm)

# Bartlett and Levene will do for the transformed variables

bartlett.test(perc_col_bc~Sp., data= dat)
bartlett.test(perc_col_bc2~Sp., data= dat)
bartlett.test(perc_col_bc3~Sp., data= dat)
bartlett.test(perc_col_bc2~Size, data= dat)
bartlett.test(perc_col_bc3~Size, data= dat)
leveneTest(perc_col_bc3~Sp.:Size, data= dat, center=median)

# We pay attention now to the dispersion of the models, due to overdispersion inflates tipe I error

phi <- sum((residuals(perc_col_bc_aov, type="pearson"))^2)/perc_col_bc_aov$df.residual
print(c("Pearson overdispersion =", round(phi, 3)), quote=FALSE)
phirob <- sum((residuals(modelo.lmrob, type="pearson"))^2)/modelo.lmrob$df.residual
print(c("Pearson overdispersion =", round(phirob, 3)), quote=FALSE)
phi2 <- sum((residuals(perc_col_bc_aov2, type="pearson"))^2)/perc_col_bc_aov2$df.residual
print(c("Pearson overdispersion =", round(phi2, 3)), quote=FALSE)
phirob2 <- sum((residuals(modelo.lmrob2, type="pearson"))^2)/modelo.lmrob2$df.residual
print(c("Pearson overdispersion =", round(phirob2, 3)), quote=FALSE)
phi3 <- sum((residuals(perc_col_bc_aov3, type="pearson"))^2)/perc_col_bc_aov3$df.residual
print(c("Pearson overdispersion =", round(phi3, 3)), quote=FALSE)
phirob3 <- sum((residuals(modelo.lmrob3, type="pearson"))^2)/modelo.lmrob3$df.residual
print(c("Pearson overdispersion =", round(phirob3, 3)), quote=FALSE)

# The full models do it better
# Finally in this first part of the statistical analyses, we look to the F and p values, corrected by sandwich estimators (hc4)
# We prefer to keep going with the robust models, concerning the overall results so far

Anova(modelo.lmrob, type=3, test="F", white.adjust="hc4")
Anova(modelo.lmrob2, type=3, test="F", white.adjust="hc4")
Anova(modelo.lmrob3, type=3, test="F", white.adjust="hc4")

## Pairwise t tests --------------------------------------------------------

## We select the full model transformed variable. We look now at the pairwise comparisons
# Pairwise comparisons
pwc <- dat %>% 
  pairwise_t_test(
    perc_col_bc3 ~ Sp., pool.sd = TRUE,
    p.adjust.method = "holm"
  ) 

pwc2 <- dat %>% 
  pairwise_t_test(
    perc_col_bc3 ~ Size, pool.sd = TRUE,
    p.adjust.method = "holm"
  ) 

pwc3 <- dat %>% 
  pairwise_t_test(
    perc_col_bc3 ~ interaccion, pool.sd = TRUE,
    p.adjust.method = "holm"
  ) 

# We calculate the means by groups

mean <- dat %>%
  group_by(Sp.) %>%
  summarise(
    mean = mean(perc_col_bc3)
  )
mean2 <- dat %>%
  group_by(Size) %>%
  summarise(
    mean = mean(perc_col_bc3)
  )
mean3 <- dat %>%
  group_by(interaccion) %>%
  summarise(
    mean = mean(perc_col_bc3)
  )

# Differences between means

mean_dif <- as.data.frame(outer(mean$mean, mean$mean, FUN="-"))
colnames(mean_dif) <- rownames(mean_dif) <- mean$Sp.

mean_dif2 <- as.data.frame(outer(mean2$mean, mean2$mean, FUN="-"))
colnames(mean_dif2) <- rownames(mean_dif2) <- mean2$Size

mean_dif3 <- as.data.frame(outer(mean3$mean, mean3$mean, FUN="-"))
colnames(mean_dif3) <- rownames(mean_dif3) <- mean3$interaccion

# We order rows by columns to make later the output of the results to coincide in order with the pwc

to_long <- function(x){
  
  df <- x[order(rownames(x)),order(colnames(x))]  
  df <- data.frame( t(combn(rownames(df),2)), dist=t(df)[lower.tri(df)])
  colnames(df)[1:3] <- c("trat1", "trat2", "mean_dif")
  return(df)
}

mean_dif_long <- to_long(mean_dif)

mean_dif_long2 <- to_long(mean_dif2)

mean_dif_long3 <- to_long(mean_dif3)

# we include the differences in the pwc

pwc$mean_dif <-mean_dif_long[,"mean_dif"]
pwc2$mean_dif <-mean_dif_long2[,"mean_dif"]
pwc3$mean_dif <-mean_dif_long3[,"mean_dif"]

# We finally incorporated the Cohen's D for the effect size and magnitudes

cohen_dif <- dat %>% cohens_d(perc_col_bc3 ~ Sp.,var.equal = FALSE)
cohen_dif <- as.data.frame(cohen_dif)
cohen_dif2 <- dat %>% cohens_d(perc_col_bc3 ~ Size,var.equal = FALSE)
cohen_dif2 <- as.data.frame(cohen_dif2)
cohen_dif3 <- dat %>% cohens_d(perc_col_bc3 ~ interaccion,var.equal = FALSE)
cohen_dif3 <- as.data.frame(cohen_dif3)

pwc$effsize <- cohen_dif[,"effsize"]
pwc$magnitude <- cohen_dif[,"magnitude"]
pwc2$effsize <- cohen_dif2[,"effsize"]
pwc2$magnitude <- cohen_dif2[,"magnitude"]
pwc3$effsize <- cohen_dif3[,"effsize"]
pwc3$magnitude <- cohen_dif3[,"magnitude"]

# We end this part exporting the results. These results will be used for creating the Compact Letter Displays of the Fig.1

pwc <- as.data.frame(pwc)
pwc2 <- as.data.frame(pwc2)
pwc3 <- as.data.frame(pwc3)

write.csv(pwc, file= "CSpwc.csv", sep= ",")
write.csv(pwc2, file= "CSpwc2.csv", sep= ",")
write.csv(pwc3, file= "CSpwc3.csv", sep= ",")

# Cleaning for next part of the code ####

rm(list = ls()) # to clean data before start

#### 1.3 - Viable biomass (%) -----------------------------------------------------------------------
## General preparation -----------------------------------------------
# Libraries ---------------------------------------------------------------

library(car) # for Levene's test, Anova and boxCox (v. 3.1-0)
library(robustbase) # for robust models (v. 0.95-0)
library(stats) # for shapiro.test, ks.test, bartlett.test, fligner.test, aov for doing the models and preparing the transformation, residuals for getting the overdispersion, cor for the R2 (v. 4.1.2)
library(rstatix) # for the pairwise test and Cohen's d (v. 0.7.0)
library(dplyr) # for managing the data and the pipelines (v. 2.2.1)

# Dataset and configuration -----------------------------------------------------------------

path.to.data <- "https://raw.githubusercontent.com/Fernando-Hurtado/DCR/main/Establishment/"
dat <- read.csv(paste0(path.to.data,"Dataset_Indicators.csv"), header=T, sep=";", dec= ",")

options(contrasts=c(factor="contr.sum", ordered="contr.poly")) ## We have factors in the models, so we change the contrast matrix as it do not corresponds with the theory in stadistics, by deafult R applies a different one

# I create the interaction column for some of the tests for the models

attach(dat)
dat$interaccion <- Sp.:Size
detach(dat)

# We create the variable Viable biomass

dat$biomVper <- dat$per_BiomV

## Statistical analyses ----------------------------------------------------
# Creation of the models 

eqt <- as.formula(biomVper ~ Sp.)
eqt2 <- as.formula(biomVper ~ Sp.+Size)
eqt3 <- as.formula(biomVper ~ Sp.*Size)


biomVper_aov <- aov(eqt, data= dat)
biomVper_aov2 <- aov(eqt2, data= dat)
biomVper_aov3 <- aov(eqt3, data= dat)

# We take a look at the residuals and the qqnorm

plot(biomVper_aov)
plot(biomVper_aov2)
plot(biomVper_aov3)

# We use Kolmogorov- Smirnov tests, as it works good for a big N

ks.test(residuals(biomVper_aov), y= pnorm)
ks.test(residuals(biomVper_aov2), y= pnorm)
ks.test(residuals(biomVper_aov3), y= pnorm)

# We look now at the homocedasticity, we use bartlett (and Levene for the interaction) as it accepts some deviation of normality

bartlett.test(biomVper~Sp., data= dat)
bartlett.test(biomVper~Size, data= dat)
leveneTest(biomVper~Sp.:Size, data= dat, center=median)

# We have to transform the variable. From this part we actually tested many transformations and finally we kept with Box-cox transformation

biomVper <- dat$biomVper + min(dat$biomVper[dat$biomVper!=0], na.rm = TRUE) # Box-cox transformations do not deal well with 0 values, so we prepare the data for the transformation
dat$biomVper <- biomVper

lambda.bc <- with(boxCox(biomVper_aov, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
lambda.bc2 <- with(boxCox(biomVper_aov2, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
lambda.bc3 <- with(boxCox(biomVper_aov3, lambda=seq(-2,2, 1/1000)), x[which.max(y)])


biomVper_bc <- ((biomVper^lambda.bc)-1)/(lambda.bc)
biomVper_bc2 <- ((biomVper^lambda.bc2)-1)/(lambda.bc2)
biomVper_bc3 <- ((biomVper^lambda.bc3)-1)/(lambda.bc3)


# we incorporate the variables, corrected for each model, to the dataset

dat$biomVper_bc <- biomVper_bc
dat$biomVper_bc2 <- biomVper_bc2
dat$biomVper_bc3 <- biomVper_bc3


# Creation of the models with the transformed variable

eqt <- as.formula(biomVper_bc ~ Sp.)
eqt2 <- as.formula(biomVper_bc2 ~ Sp.+Size)
eqt3 <- as.formula(biomVper_bc3 ~ Sp.*Size)


biomVper_bc_aov <- aov(eqt, data= dat)
biomVper_bc_aov2 <- aov(eqt2, data= dat)
biomVper_bc_aov3 <- aov(eqt3, data= dat)

# Creation of the robust models 

modelo.lmrob <- lmrob(eqt, data=dat, method="SMDM", setting = "KS2014")
modelo.lmrob2 <- lmrob(eqt2, data=dat, method="SMDM", setting = "KS2014")
modelo.lmrob3 <- lmrob(eqt3, data=dat, method="SMDM", setting = "KS2014")

# Analyses of the models, first a look at the residuals and qqnorm

plot(biomVper_bc_aov)
plot(modelo.lmrob)
plot(biomVper_bc_aov2)
plot(modelo.lmrob2)
plot(biomVper_bc_aov3)
plot(modelo.lmrob3)

# Even it is difficult to do the interpretation, it is an improvement in general in normality and more compensated leverages in the robust models.

# We go with the tests 

ks.test(residuals(biomVper_bc_aov), y= pnorm)
ks.test(residuals(biomVper_bc_aov2), y= pnorm)
ks.test(residuals(biomVper_bc_aov3), y= pnorm)
ks.test(residuals(modelo.lmrob), y= pnorm)
ks.test(residuals(modelo.lmrob2), y= pnorm)
ks.test(residuals(modelo.lmrob3), y= pnorm)

# Bartlett and Levene will do for the transformed variables

bartlett.test(biomVper_bc~Sp., data= dat)
bartlett.test(biomVper_bc2~Sp., data= dat)
bartlett.test(biomVper_bc3~Sp., data= dat)
bartlett.test(biomVper_bc2~Size, data= dat)
bartlett.test(biomVper_bc3~Size, data= dat)
leveneTest(biomVper_bc3~Sp.:Size, data= dat, center=median)

# We pay attention now to the dispersion of the models, due to overdispersion inflates tipe I error

phi <- sum((residuals(biomVper_bc_aov, type="pearson"))^2)/biomVper_bc_aov$df.residual
print(c("Pearson overdispersion =", round(phi, 3)), quote=FALSE)
phirob <- sum((residuals(modelo.lmrob, type="pearson"))^2)/modelo.lmrob$df.residual
print(c("Pearson overdispersion =", round(phirob, 3)), quote=FALSE)
phi2 <- sum((residuals(biomVper_bc_aov2, type="pearson"))^2)/biomVper_bc_aov2$df.residual
print(c("Pearson overdispersion =", round(phi2, 3)), quote=FALSE)
phirob2 <- sum((residuals(modelo.lmrob2, type="pearson"))^2)/modelo.lmrob2$df.residual
print(c("Pearson overdispersion =", round(phirob2, 3)), quote=FALSE)
phi3 <- sum((residuals(biomVper_bc_aov3, type="pearson"))^2)/biomVper_bc_aov3$df.residual
print(c("Pearson overdispersion =", round(phi3, 3)), quote=FALSE)
phirob3 <- sum((residuals(modelo.lmrob3, type="pearson"))^2)/modelo.lmrob3$df.residual
print(c("Pearson overdispersion =", round(phirob3, 3)), quote=FALSE)

# The full models do it better
# Finally in this first part of the statistical analyses, we look to the F and p values, corrected by sandwich estimators (hc4)
# We prefer to keep going with the robust models, concerning the overall results so far

Anova(modelo.lmrob, type=3, test="F", white.adjust="hc4")
Anova(modelo.lmrob2, type=3, test="F", white.adjust="hc4")
Anova(modelo.lmrob3, type=3, test="F", white.adjust="hc4")


## Pairwise t tests --------------------------------------------------------

## We select the full model transformed variable. We look now at the pairwise comparisons
# Pairwise comparisons

pwc <- dat %>% 
  pairwise_t_test(
    biomVper_bc3 ~ Sp., pool.sd = TRUE,
    p.adjust.method = "holm"
  ) 

pwc2 <- dat %>% 
  pairwise_t_test(
    biomVper_bc3 ~ Size, pool.sd = TRUE,
    p.adjust.method = "holm"
  ) 

pwc3 <- dat %>% 
  pairwise_t_test(
    biomVper_bc3 ~ interaccion, pool.sd = TRUE,
    p.adjust.method = "holm"
  ) 

# We calculate the means by groups

mean <- dat %>%
  group_by(Sp.) %>%
  summarise(
    mean = mean(biomVper_bc3)
  )
mean2 <- dat %>%
  group_by(Size) %>%
  summarise(
    mean = mean(biomVper_bc3)
  )
mean3 <- dat %>%
  group_by(interaccion) %>%
  summarise(
    mean = mean(biomVper_bc3)
  )

# Differences between means

mean_dif <- as.data.frame(outer(mean$mean, mean$mean, FUN="-"))
colnames(mean_dif) <- rownames(mean_dif) <- mean$Sp.

mean_dif2 <- as.data.frame(outer(mean2$mean, mean2$mean, FUN="-"))
colnames(mean_dif2) <- rownames(mean_dif2) <- mean2$Size

mean_dif3 <- as.data.frame(outer(mean3$mean, mean3$mean, FUN="-"))
colnames(mean_dif3) <- rownames(mean_dif3) <- mean3$interaccion

# We order rows by columns to make later the output of the results to coincide in order with the pwc

to_long <- function(x){
  
  df <- x[order(rownames(x)),order(colnames(x))]  
  df <- data.frame( t(combn(rownames(df),2)), dist=t(df)[lower.tri(df)])
  colnames(df)[1:3] <- c("trat1", "trat2", "mean_dif")
  return(df)
}

mean_dif_long <- to_long(mean_dif)

mean_dif_long2 <- to_long(mean_dif2)

mean_dif_long3 <- to_long(mean_dif3)

# we include the differences in the pwc

pwc$mean_dif <-mean_dif_long[,"mean_dif"]
pwc2$mean_dif <-mean_dif_long2[,"mean_dif"]
pwc3$mean_dif <-mean_dif_long3[,"mean_dif"]

# We finally incorporated the Cohen's D for the effect size and magnitudes

cohen_dif <- dat %>% cohens_d(biomVper_bc3 ~ Sp.,var.equal = FALSE)
cohen_dif <- as.data.frame(cohen_dif)
cohen_dif2 <- dat %>% cohens_d(biomVper_bc3 ~ Size,var.equal = FALSE)
cohen_dif2 <- as.data.frame(cohen_dif2)
cohen_dif3 <- dat %>% cohens_d(biomVper_bc3 ~ interaccion,var.equal = FALSE)
cohen_dif3 <- as.data.frame(cohen_dif3)

pwc$effsize <- cohen_dif[,"effsize"]
pwc$magnitude <- cohen_dif[,"magnitude"]
pwc2$effsize <- cohen_dif2[,"effsize"]
pwc2$magnitude <- cohen_dif2[,"magnitude"]
pwc3$effsize <- cohen_dif3[,"effsize"]
pwc3$magnitude <- cohen_dif3[,"magnitude"]

# We end this part exporting the results. These results will be used for creating the Compact Letter Displays of the Fig.1

pwc <- as.data.frame(pwc)
pwc2 <- as.data.frame(pwc2)
pwc3 <- as.data.frame(pwc3)

write.csv(pwc, file= "VBpwc.csv", sep= ",")
write.csv(pwc2, file= "VBpwc2.csv", sep= ",")
write.csv(pwc3, file= "VBpwc3.csv", sep= ",")

# Cleaning for next part of the code ####

rm(list = ls()) # to clean data before start

#### 1.4 - Relative Growth Rate (g/month) ---------------------------------------------------------------------
## General preparation -----------------------------------------------
# Libraries ---------------------------------------------------------------

library(car) # for Levene's test, Anova and boxCox (v.3.1-0) 
library(robustbase) # for robust models (v. 0.95-0)
library(stats) # for shapiro.test, ks.test, bartlett.test, fligner.test, aov for doing the models and preparing the transformation, residuals for getting the overdispersion, cor for the R2 (v. 4.1.2)
library(rstatix) # for the pairwise test and Cohen's d (v. 0.7.0)
library(dplyr) # for managing the data and the pipelines (v. 2.2.1)

# Dataset and configuration -----------------------------------------------------------------

path.to.data <- "https://raw.githubusercontent.com/Fernando-Hurtado/DCR/main/Establishment/"
dat <- read.csv(paste0(path.to.data,"Dataset_Indicators.csv"), header=T, sep=";", dec= ",")

options(contrasts=c(factor="contr.sum", ordered="contr.poly")) ## We have factors in the models, so we change the contrast matrix as it do not corresponds with the theory in stadistics, by deafult R applies a different one

# I create the interaction column for some of the tests for the models

attach(dat)
dat$interaccion <- Sp.:Size
detach(dat)

# We create the variable RGR

dat$RGR <- (log(dat$Biom_total) - log(0.0040))/2

## Statistical analyses ----------------------------------------------------

# Creation of the models 

eqt <- as.formula(RGR ~ Sp.)
eqt2 <- as.formula(RGR ~ Sp.+Size)
eqt3 <- as.formula(RGR ~ Sp.*Size)

RGR_aov <- aov(eqt, data= dat)
RGR_aov2 <- aov(eqt2, data= dat)
RGR_aov3 <- aov(eqt3, data= dat)

# We take a look at the residuals and the qqnorm

plot(RGR_aov)
plot(RGR_aov2)
plot(RGR_aov3)

# We use Kolmogorov- Smirnov tests, as it works good for a big N

ks.test(residuals(RGR_aov), y= pnorm)
ks.test(residuals(RGR_aov2), y= pnorm)
ks.test(residuals(RGR_aov3), y= pnorm)

# We look now at the homocedasticity, we use bartlett (and Levene for the interaction) as it accepts some deviation of normality

bartlett.test(RGR~Sp., data= dat)
bartlett.test(RGR~Size, data= dat)
leveneTest(RGR~Sp.:Size, data= dat, center=median)

# RGR transformation was not needed, as we tested many transformations from this point and all of them with really few impact and in the case of the former variable, being acceptable in overall terms.

# We create the robust models

modelo.lmrob <- lmrob(eqt, data=dat, method="SMDM", setting = "KS2014")
modelo.lmrob2 <- lmrob(eqt2, data=dat, method="SMDM", setting = "KS2014")
modelo.lmrob3 <- lmrob(eqt3, data=dat, method="SMDM", setting = "KS2014")

# We pay attention now to the dispersion of the models, due to overdispersion inflates tipe I error

phi <- sum((residuals(RGR_aov, type="pearson"))^2)/RGR_aov$df.residual
print(c("Pearson overdispersion =", round(phi, 3)), quote=FALSE)
phirob <- sum((residuals(modelo.lmrob, type="pearson"))^2)/modelo.lmrob$df.residual
print(c("Pearson overdispersion =", round(phirob, 3)), quote=FALSE)
phi2 <- sum((residuals(RGR_aov2, type="pearson"))^2)/RGR_aov2$df.residual
print(c("Pearson overdispersion =", round(phi2, 3)), quote=FALSE)
phirob2 <- sum((residuals(modelo.lmrob2, type="pearson"))^2)/modelo.lmrob2$df.residual
print(c("Pearson overdispersion =", round(phirob2, 3)), quote=FALSE)
phi3 <- sum((residuals(RGR_aov3, type="pearson"))^2)/RGR_aov3$df.residual
print(c("Pearson overdispersion =", round(phi3, 3)), quote=FALSE)
phirob3 <- sum((residuals(modelo.lmrob3, type="pearson"))^2)/modelo.lmrob3$df.residual
print(c("Pearson overdispersion =", round(phirob3, 3)), quote=FALSE)

# The full models do it better
# Finally in this first part of the statistical analyses, we look to the F and p values, corrected by sandwich estimators (hc4)
# We prefer to keep going with the robust models, concerning the overall results so far

Anova(modelo.lmrob, type=3, test="F", white.adjust="hc4")
Anova(modelo.lmrob2, type=3, test="F", white.adjust="hc4")
Anova(modelo.lmrob3, type=3, test="F", white.adjust="hc4")


## Pairwise t tests --------------------------------------------------------

## We select the full model transformed variable. We look now at the pairwise comparisons
# Pairwise comparisons

pwc <- dat %>% 
  pairwise_t_test(
    RGR ~ Sp., pool.sd = TRUE,
    p.adjust.method = "holm"
  ) 

pwc2 <- dat %>% 
  pairwise_t_test(
    RGR ~ Size, pool.sd = TRUE,
    p.adjust.method = "holm"
  ) 

pwc3 <- dat %>% 
  pairwise_t_test(
    RGR ~ interaccion, pool.sd = TRUE,
    p.adjust.method = "holm"
  ) 

# We calculate the means by groups

mean <- dat %>%
  group_by(Sp.) %>%
  summarise(
    mean = mean(RGR)
  )
mean2 <- dat %>%
  group_by(Size) %>%
  summarise(
    mean = mean(RGR)
  )
mean3 <- dat %>%
  group_by(interaccion) %>%
  summarise(
    mean = mean(RGR)
  )

# Differences between means

mean_dif <- as.data.frame(outer(mean$mean, mean$mean, FUN="-"))
colnames(mean_dif) <- rownames(mean_dif) <- mean$Sp.

mean_dif2 <- as.data.frame(outer(mean2$mean, mean2$mean, FUN="-"))
colnames(mean_dif2) <- rownames(mean_dif2) <- mean2$Size

mean_dif3 <- as.data.frame(outer(mean3$mean, mean3$mean, FUN="-"))
colnames(mean_dif3) <- rownames(mean_dif3) <- mean3$interaccion

# We order rows by columns to make later the output of the results to coincide in order with the pwc

to_long <- function(x){
  
  df <- x[order(rownames(x)),order(colnames(x))]  
  df <- data.frame( t(combn(rownames(df),2)), dist=t(df)[lower.tri(df)])
  colnames(df)[1:3] <- c("trat1", "trat2", "mean_dif")
  return(df)
}

mean_dif_long <- to_long(mean_dif)

mean_dif_long2 <- to_long(mean_dif2)

mean_dif_long3 <- to_long(mean_dif3)

# we include the differences in the pwc

pwc$mean_dif <-mean_dif_long[,"mean_dif"]
pwc2$mean_dif <-mean_dif_long2[,"mean_dif"]
pwc3$mean_dif <-mean_dif_long3[,"mean_dif"]

# We finally incorporated the Cohen's D for the effect size and magnitudes

cohen_dif <- dat %>% cohens_d(RGR ~ Sp.,var.equal = FALSE)
cohen_dif <- as.data.frame(cohen_dif)
cohen_dif2 <- dat %>% cohens_d(RGR ~ Size,var.equal = FALSE)
cohen_dif2 <- as.data.frame(cohen_dif2)
cohen_dif3 <- dat %>% cohens_d(RGR ~ interaccion,var.equal = FALSE)
cohen_dif3 <- as.data.frame(cohen_dif3)

pwc$effsize <- cohen_dif[,"effsize"]
pwc$magnitude <- cohen_dif[,"magnitude"]
pwc2$effsize <- cohen_dif2[,"effsize"]
pwc2$magnitude <- cohen_dif2[,"magnitude"]
pwc3$effsize <- cohen_dif3[,"effsize"]
pwc3$magnitude <- cohen_dif3[,"magnitude"]

# We end this part exporting the results. These results will be used for creating the Compact Letter Displays of the Fig.1

pwc <- as.data.frame(pwc)
pwc2 <- as.data.frame(pwc2)
pwc3 <- as.data.frame(pwc3)

write.csv2(pwc, file= "RGRcpwc.csv", sep= ",")
write.csv2(pwc2, file= "RGRcpwc2.csv", sep= ",")
write.csv2(pwc3, file= "RGRcpwc3.csv", sep= ",")











# Cleaning for next part of the code ####

rm(list = ls()) # to clean data before start

#### 1.5 - Table S.3 Wilcoxon tests ---------------------------------------------------------------------
# Libraries ---------------------------------------------------------------

library(rstatix) # for the wilcox_test (v. 0.7.0)

# Dataset -----------------------------------------------------------------

path.to.data <- "https://raw.githubusercontent.com/Fernando-Hurtado/DCR/main/Establishment/"
dat <- read.csv(paste0(path.to.data,"Dataset_Traits.csv"), header=T, sep=";", dec= ",")

# Subsets -------------------------------------------------------------

dat$Size <- as.factor(dat$Size)

levels(dat$Size)[levels(dat$Size)=="Large"] <- "L"
levels(dat$Size)[levels(dat$Size)=="Medium"] <- "M"
levels(dat$Size)[levels(dat$Size)=="Small"] <- "S"

dat$Area_wet <- (dat$Area_wet)*100
dat$Area_dry <- (dat$Area_dry)*100
dat$Area_diff <- (dat$Area_wet - dat$Area_dry)/dat$Area_wet
dat$Interaction <- paste(dat$Sp.,dat$Size,sep=" ")

# Wilcoxon tests for the table S.4 ----------------------------------------------------------

cm.data.values <- wilcox_test(Area_diff~ Interaction,  data = dat, ref.group = ".all.")

# Cleaning for next part of the code ####

rm(list = ls()) # to clean data before start

#### 1.6 - Table S.4 Pearson's Chi-squared tests  ---------------------------------------------------------------------
## General preparation -----------------------------------------------
# Libraries ---------------------------------------------------------------

library(vcd) # for the use of residual base shading with mosaic (v. 1.4-10)

# Dataset -----------------------------------------------------------------

path.to.data <- "https://raw.githubusercontent.com/Fernando-Hurtado/DCR/main/Establishment/"
dat <- read.csv(paste0(path.to.data,"Dataset_Traits.csv"), header=T, sep=";", dec= ",")


# Subsets -------------------------------------------------------------

dat$Size <- as.factor(dat$Size)

levels(dat$Size)[levels(dat$Size)=="Large"] <- "L"
levels(dat$Size)[levels(dat$Size)=="Medium"] <- "M"
levels(dat$Size)[levels(dat$Size)=="Small"] <- "S"

vectorsp1 <- which(dat[,1]== "D. scoparium")
vectorsp2 <- which(dat[,1]== "H. aureum")
vectorsp3 <- which(dat[,1]== "H. cupressiforme")
vectorsp4 <- which(dat[,1]== "P. capillare")
vectorsp5 <- which(dat[,1]== "S. ruralis")
vectorsp6 <- which(dat[,1]== "T. squarrosa")

newdatasp1 <- dat[c(vectorsp1),]
newdatasp2 <- dat[c(vectorsp2),]
newdatasp3 <- dat[c(vectorsp3),]
newdatasp4 <- dat[c(vectorsp4),]
newdatasp5 <- dat[c(vectorsp5),]
newdatasp6 <- dat[c(vectorsp6),]

# chisq.tests -------------------------------------------------------------

Xsq <- chisq.test(newdatasp1$Size, newdatasp1$Type=="S")
Xsq2 <- chisq.test(newdatasp2$Size, newdatasp2$Type=="S")
Xsq3 <- chisq.test(newdatasp3$Size, newdatasp3$Type=="S")
Xsq4 <- chisq.test(newdatasp4$Size, newdatasp4$Type=="S")
Xsq5 <- chisq.test(newdatasp5$Size, newdatasp5$Type=="S")
Xsq6 <- chisq.test(newdatasp6$Size, newdatasp6$Type=="S")

# The p-values for the Table S.4 ------------------------------------------

Xsq$p.value 
Xsq2$p.value
Xsq3$p.value
Xsq4$p.value
Xsq5$p.value
Xsq6$p.value

# Pearson residuals by species and propagule size ----------------------------------------------------

# Comment for the use of Residual-based shading: The colors of the outputs represent significant differences we can gather after plotting
# They are used for the rest of the columns of the figure S.4

# Explanation for the interpretation of the colors of the mosaics:
      # The hue indicates the residuals sign: by default, blue for positive, and red for negative residuals. The
      # saturation of a residual is set according to its size: high saturation for large, and low saturation
      # for small residuals. Finally, the overall lightness is used to indicate the significance of a test
      # statistic: light colors for significant, and dark colors for non-significant results.
      # The heuristic for choosing the cut-off points 2 and 4 is that the Pearson residuals are approximately standard normal which
      # implies that the highlighted cells are those with residuals individually significant at approximately
      # the alpha ? = 0.05 and alpha ? = 0.0001 levels, respectively.

mosaic(~ Size + Type,
       direction = c("v", "h"),
       data = newdatasp1,
       shade = TRUE)

# "large" is bigger than expected(*) and "small" is smaller than expected(*)

mosaic(~ Size + Type,
       direction = c("v", "h"),
       data = newdatasp2,
       shade = TRUE)

# same as above

mosaic(~ Size + Type,
       direction = c("v", "h"),
       data = newdatasp3,
       shade = TRUE)

# same as above but "large" is much bigger than expected (**)

mosaic(~ Size + Type,
       direction = c("v", "h"),
       data = newdatasp4,
       shade = TRUE)

# no differences by groups

mosaic(~ Size + Type,
       direction = c("v", "h"),
       data = newdatasp5,
       shade = TRUE)

# "large" is bigger than expected (*)

mosaic(~ Size + Type,
       direction = c("v", "h"),
       data = newdatasp6,
       shade = TRUE,
)

# "large" is much bigger than expected (**) and "medium" and "small", less than expected (*)
# Cleaning for next part of the code ####

rm(list = ls()) # to clean data before start

############ 2 - Main figure (Figure 1) --------------------------------------------------------------------


### 2.1 - General preparation -----------------------------------------------
# Comments ####
  # The main figure is a mosaic plot integrated by 12 plots with two sets of compact letter displays (CLD, like tags) to annotate the significant differences across species within the same propagule size and across propagule sizes within the same species
  # Each plot is for each indicator(variable) and propagule size class 
# Libraries ####

library(grid) # for the annotations in the figures (v. 4.1.2)
library(ggpubr) # for the use of ggboxplot (v. 0.4.0)
library(cowplot) # for the mosaic plot (v. 1.1.1)
library(grDevices) # for using png for saving the images (v. 4.1.2)
library(car) # for Levene's test, Anova and boxCox (v.3.1-0)
library(stats) # for shapiro.test, ks.test, bartlett.test, fligner.test, aov for doing the models and preparing the transformation, residuals for getting the overdispersion, cor for the R2 (v. 4.1.2)

# Dataset ####

path.to.data <- "https://raw.githubusercontent.com/Fernando-Hurtado/DCR/main/Establishment/"
dat <- read.csv(paste0(path.to.data,"Dataset_Indicators.csv"), header=T, sep=";", dec= ",")

## Incorporation of the transformed indicators and RGR (non-transformed) -----------------------------

# Number of established propagules (propagules/mm2) ------------------------------------------------------

options(contrasts=c(factor="contr.sum", ordered="contr.poly")) ## We have factors in the models, so we change the contrast matrix as it do not corresponds with the theory in stadistics, by deafult R applies a different one
eqt3 <- as.formula(Shoots_established ~ Sp.*Size)
Shoots_established_aov3 <- aov(eqt3, data= dat)
Shoots_established <- dat$Shoots_established + min(dat$Shoots_established[dat$Shoots_established!=0], na.rm = TRUE) # Box-cox transformations do not deal well with 0 values, so we prepare the data for the transformation
dat$Shoots_established <- Shoots_established
lambda.bc3 <- with(boxCox(Shoots_established_aov3, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
Shoots_established_bc3 <- ((Shoots_established^lambda.bc3)-1)/(lambda.bc3)
dat$Shoots_established_bc3 <- Shoots_established_bc3

# Colonized surface -------------------------------------------------------

options(contrasts=c(factor="contr.sum", ordered="contr.poly")) ## We have factors in the models, so we change the contrast matrix as it do not corresponds with the theory in stadistics, by deafult R applies a different one
dat$perc_col <- dat$S_Shootss/dat$S_substratum
dat$perc_col <- dat$perc_col*100
eqt3 <- as.formula(perc_col ~ Sp.*Size)
perc_col_aov3 <- aov(eqt3, data= dat)
perc_col <- dat$perc_col + min(dat$perc_col[dat$perc_col!=0], na.rm = TRUE) # Box-cox transformations do not deal well with 0 values, so we prepare the data for the transformation
dat$perc_col <- perc_col
lambda.bc3 <- with(boxCox(perc_col_aov3, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
perc_col_bc3 <- ((perc_col^lambda.bc3)-1)/(lambda.bc3)
dat$perc_col_bc3 <- perc_col_bc3

# Viable biomass ----------------------------------------------------------

options(contrasts=c(factor="contr.sum", ordered="contr.poly")) ## We have factors in the models, so we change the contrast matrix as it do not corresponds with the theory in stadistics, by deafult R applies a different one
dat$biomVper <- dat$per_BiomV
eqt3 <- as.formula(biomVper ~ Sp.*Size)
biomVper_aov3 <- aov(eqt3, data= dat)
biomVper <- dat$biomVper + min(dat$biomVper[dat$biomVper!=0], na.rm = TRUE) # Box-cox transformations do not deal well with 0 values, so we prepare the data for the transformation
dat$biomVper <- biomVper
lambda.bc3 <- with(boxCox(biomVper_aov3, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
biomVper_bc3 <- ((biomVper^lambda.bc3)-1)/(lambda.bc3)
dat$biomVper_bc3 <- biomVper_bc3

# RGR ---------------------------------------------------------------------

# We create the variable RGR. It does not need transformation

dat$RGR <- (log(dat$Biom_total) - log(0.0040))/2

# Functions ####

## get_plot_limits from QsRutils would be quite handy here

get_plot_limits <- function(plot) {
  gb = ggplot2::ggplot_build(plot)
  xmin = gb$layout$panel_params[[1]]$x.range[1]
  xmax = gb$layout$panel_params[[1]]$x.range[2]
  ymin = gb$layout$panel_params[[1]]$y.range[1]
  ymax = gb$layout$panel_params[[1]]$y.range[2]
  list(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
}

# Subsets ####

# We define the subsets for the figures (factor levels of the treatment)

vectorsize1 <- which(dat[,10]=="Large")
vectorsize2 <- which(dat[,10]=="Medium")
vectorsize3 <- which(dat[,10]=="Small")

# change in names to two letters for better representation later in the figures

dat$Sp.<-as.factor(dat$Sp.)

levels(dat$Sp.)[levels(dat$Sp.)=="D. scoparium"] <- "Ds"
levels(dat$Sp.)[levels(dat$Sp.)=="H. aureum"] <- "Ha"
levels(dat$Sp.)[levels(dat$Sp.)=="H. cupressiforme"] <- "Hc"
levels(dat$Sp.)[levels(dat$Sp.)=="P. capillare"] <- "Pc"
levels(dat$Sp.)[levels(dat$Sp.)=="S. ruralis"] <- "Sr"
levels(dat$Sp.)[levels(dat$Sp.)=="T. squarrosa"] <- "Ts"


### 2.2 - Number of established propagules (propagules/mm2) ####

# We define the boxplot first to calculate later the plot limits

pltsp1 <- ggboxplot(dat[c(vectorsize1),],x= "Sp.", y= "Shoots_established_bc3") +     geom_boxplot()
pltsp2 <- ggboxplot(dat[c(vectorsize2),],x= "Sp.", y= "Shoots_established_bc3") +     geom_boxplot()
pltsp3 <- ggboxplot(dat[c(vectorsize3),],x= "Sp.", y= "Shoots_established_bc3") +     geom_boxplot()

# We use the coordinates that the boxplot give us for the upper whiskers to use it as reference for the compact letter display

box.rsltsp1 <- with(dat[c(vectorsize1),], graphics::boxplot(Shoots_established_bc3 ~ Sp., plot = FALSE))
box.rsltsp2 <- with(dat[c(vectorsize2),], graphics::boxplot(Shoots_established_bc3 ~ Sp., plot = FALSE))
box.rsltsp3 <- with(dat[c(vectorsize3),], graphics::boxplot(Shoots_established_bc3 ~ Sp., plot = FALSE))

# Compact Letters Display (CLD) ####

# CLD per Sizes

xsp1 <- c(1:6)
ysp1 <- box.rsltsp1$stats[5, ]
cbdsp1 <- c("a", "b","b","b","b","c")
ltr_dfsp1 <- data.frame(xsp1, ysp1, cbdsp1)


lmtssp1 <- get_plot_limits(pltsp1)
y.rangesp1 <- lmtssp1$ymax - lmtssp1$ymin
y.nudgesp1 <- 0.05 * y.rangesp1

xsp2 <- c(1:6)
ysp2 <- box.rsltsp2$stats[5, ]
cbdsp2 <- c("a","b","c","b","c","d")
ltr_dfsp2 <- data.frame(xsp2, ysp2, cbdsp2)

lmtssp2 <- get_plot_limits(pltsp2)
y.rangesp2 <- lmtssp2$ymax - lmtssp2$ymin
y.nudgesp2 <- 0.05 * y.rangesp2

xsp3 <- c(1:6)
ysp3 <- box.rsltsp3$stats[5, ]
cbdsp3 <- c("a", "b","ac","b","c","d")
ltr_dfsp3 <- data.frame(xsp3, ysp3, cbdsp3)

lmtssp3 <- get_plot_limits(pltsp3)
y.rangesp3 <- lmtssp3$ymax - lmtssp3$ymin
y.nudgesp3 <- 0.05 * y.rangesp3

# CLD alpha beta gamma ####

# We use unicode for greek letters to avoid conflicts and problems when revisiting the script

xABCsp1 <- c(0.6,1.6,2.6,3.6,4.6,5.6)
yABCsp1 <- c(7.7,7.7,7.7,7.7,7.7,7.7)
cABCsp1 <- c("\u03B1","\u03B1","\u03B1","\u03B1","\u03B1","\u03B1")
ltr_dfABCsp1 <- data.frame(xABCsp1, yABCsp1, cABCsp1)

xABCsp2 <- c(0.6,1.6,2.6,3.6,4.6,5.6)
yABCsp2 <- c(7.7,7.7,7.7,7.7,7.7,7.7)
cABCsp2 <- c("\u03B1\u03B2","\u03B2","\u03B2","\u03B1","\u03B2","\u03B2")
ltr_dfABCsp2 <- data.frame(xABCsp2, yABCsp2, cABCsp2)

xABCsp3 <- c(0.6,1.6,2.6,3.6,4.6,5.6)
yABCsp3 <- c(7.7,7.7,7.7,7.7,7.7,7.7)
cABCsp3 <- c("\u03B3","\u03B1\u03B2","\u03B3","\u03B2","\u03B3","\u03B2")
ltr_dfABCsp3 <- data.frame(xABCsp3, yABCsp3, cABCsp3)

### 2.3 - Colonized surface (%) ####

# We define the boxplot first to calculate later the plot limits

pltsp4 <- ggboxplot(dat[c(vectorsize1),],x= "Sp.", y= "perc_col_bc3") +     geom_boxplot()
pltsp5 <- ggboxplot(dat[c(vectorsize2),],x= "Sp.", y= "perc_col_bc3") +     geom_boxplot()
pltsp6 <- ggboxplot(dat[c(vectorsize3),],x= "Sp.", y= "perc_col_bc3") +     geom_boxplot()

## We use the coordinates that the boxplot give us for the upper whiskers to use it as reference for the compact letter display

box.rsltsp4 <- with(dat[c(vectorsize1),], graphics::boxplot(perc_col_bc3 ~ Sp., plot = FALSE))
box.rsltsp5 <- with(dat[c(vectorsize2),], graphics::boxplot(perc_col_bc3 ~ Sp., plot = FALSE))
box.rsltsp6 <- with(dat[c(vectorsize3),], graphics::boxplot(perc_col_bc3 ~ Sp., plot = FALSE))

# CLD ####

# CLD per Sizes

xsp4 <- c(1:6)
ysp4 <- box.rsltsp4$stats[5, ]
cbdsp4 <- c("a", "b","b","b","bc","c")
ltr_dfsp4 <- data.frame(xsp4, ysp4, cbdsp4)

lmtssp4 <- get_plot_limits(pltsp4)
y.rangesp4 <- lmtssp4$ymax - lmtssp4$ymin
y.nudgesp4 <- 0.05 * y.rangesp4

xsp5 <- c(1:6)
ysp5 <- box.rsltsp5$stats[5, ]
cbdsp5 <- c("a","b","c","b","c","b")
ltr_dfsp5 <- data.frame(xsp5, ysp5, cbdsp5)

lmtssp5 <- get_plot_limits(pltsp5)
y.rangesp5 <- lmtssp5$ymax - lmtssp5$ymin
y.nudgesp5 <- 0.05 * y.rangesp5

xsp6 <- c(1:6)
ysp6 <- box.rsltsp6$stats[5, ]
cbdsp6 <- c("a", "b","a","b","a","c")
ltr_dfsp6 <- data.frame(xsp6, ysp6, cbdsp6)

lmtssp6 <- get_plot_limits(pltsp6)
y.rangesp6 <- lmtssp6$ymax - lmtssp6$ymin
y.nudgesp6 <- 0.05 * y.rangesp6

# CLD alpha beta gamma ####

# We use unicode for greek letters to avoid conflicts and problems when revisiting the script

xABCsp4 <- c(0.6,1.6,2.6,3.6,4.6,5.6)
yABCsp4 <- c(4.7,4.7,4.7,4.7,4.7,4.7)
cABCsp4 <- c("\u03B1", "\u03B1","\u03B1","\u03B1","\u03B1","\u03B1")
ltr_dfABCsp4 <- data.frame(xABCsp4, yABCsp4, cABCsp4)


xABCsp5 <- c(0.6,1.6,2.6,3.6,4.6,5.6)
yABCsp5 <- c(4.7,4.7,4.7,4.7,4.7,4.7)
cABCsp5 <- c("\u03B2","\u03B1","\u03B2","\u03B1","\u03B2","\u03B1")
ltr_dfABCsp5 <- data.frame(xABCsp5, yABCsp5, cABCsp5)

xABCsp6 <- c(0.6,1.6,2.6,3.6,4.6,5.6)
yABCsp6 <- c(4.7,4.7,4.7,4.7,4.7,4.7)
cABCsp6 <- c("\u03B2","\u03B1","\u03B3","\u03B2","\u03B3","\u03B1")
ltr_dfABCsp6 <- data.frame(xABCsp6, yABCsp6, cABCsp6)

### 2.4 - Viable biomass (%) ####

# We define the boxplot first to calculate later the plot limits

pltsp7 <- ggboxplot(dat[c(vectorsize1),],x= "Sp.", y= "biomVper_bc3") +     geom_boxplot()
pltsp8 <- ggboxplot(dat[c(vectorsize2),],x= "Sp.", y= "biomVper_bc3") +     geom_boxplot()
pltsp9 <- ggboxplot(dat[c(vectorsize3),],x= "Sp.", y= "biomVper_bc3") +     geom_boxplot()


## We use the coordinates that the boxplot give us for the upper whiskers to use it as reference for the compact letter display

box.rsltsp7 <- with(dat[c(vectorsize1),], graphics::boxplot(biomVper_bc3 ~ Sp., plot = FALSE))
box.rsltsp8 <- with(dat[c(vectorsize2),], graphics::boxplot(biomVper_bc3 ~ Sp., plot = FALSE))
box.rsltsp9 <- with(dat[c(vectorsize3),], graphics::boxplot(biomVper_bc3 ~ Sp., plot = FALSE))

# CLD ####

# CLD per sizes

xsp7 <- c(1:6)
ysp7 <- box.rsltsp7$stats[5, ]
cbdsp7 <- c("a", "bc","b","c","bc","c")
ltr_dfsp7 <- data.frame(xsp7, ysp7, cbdsp7)

lmtssp7 <- get_plot_limits(pltsp7)
y.rangesp7 <- lmtssp7$ymax - lmtssp7$ymin
y.nudgesp7 <- 0.15 * y.rangesp7

xsp8 <- c(1:6)
ysp8 <- box.rsltsp8$stats[5, ]
cbdsp8 <- c("a","b","c","bd","ce","de")
ltr_dfsp8 <- data.frame(xsp8, ysp8, cbdsp8)

lmtssp8 <- get_plot_limits(pltsp8)
y.rangesp8 <- lmtssp8$ymax - lmtssp8$ymin
y.nudgesp8 <- 0.15 * y.rangesp8

xsp9 <- c(1:6)
ysp9 <- box.rsltsp9$stats[5, ]
cbdsp9 <- c("a", "bc","a","b","a","c")
ltr_dfsp9 <- data.frame(xsp9, ysp9, cbdsp9)

lmtssp9 <- get_plot_limits(pltsp9)
y.rangesp9 <- lmtssp9$ymax - lmtssp9$ymin
y.nudgesp9 <- 0.15 * y.rangesp9


# CLD alpha beta gamma ####

# We use unicode for greek letters to avoid conflicts and problems when revisiting the script

xABCsp7 <- c(0.6,1.6,2.6,3.6,4.6,5.6)
yABCsp7 <- c(3.2,3.2,3.2,3.2,3.2,3.2)
cABCsp7 <- c("\u03B1", "\u03B1","\u03B1","\u03B1","\u03B1","\u03B1")
ltr_dfABCsp7 <- data.frame(xABCsp7, yABCsp7, cABCsp7)

xABCsp8 <- c(0.6,1.6,2.6,3.6,4.6,5.6)
yABCsp8 <- c(3.2,3.2,3.2,3.2,3.2,3.2)
cABCsp8 <- c("\u03B1","\u03B1","\u03B2","\u03B1","\u03B2","\u03B2")
ltr_dfABCsp8 <- data.frame(xABCsp8, yABCsp8, cABCsp8)

xABCsp9 <- c(0.6,1.6,2.6,3.6,4.6,5.6)
yABCsp9 <- c(3.2,3.2,3.2,3.2,3.2,3.2)
cABCsp9 <- c("\u03B1","\u03B2","\u03B3","\u03B2","\u03B3","\u03B2")
ltr_dfABCsp9 <- data.frame(xABCsp9, yABCsp9, cABCsp9)

### 2.5 - Relative Growth Rate (g/month) ####

# We define the boxplot first to calculate later the plot limits

pltsp10 <- ggboxplot(dat[c(vectorsize1),],x= "Sp.", y= "RGR") +     geom_boxplot()
pltsp11 <- ggboxplot(dat[c(vectorsize2),],x= "Sp.", y= "RGR") +     geom_boxplot()
pltsp12 <- ggboxplot(dat[c(vectorsize3),],x= "Sp.", y= "RGR") +     geom_boxplot()


## We use the coordinates that the boxplot give us for the upper whiskers to use it as reference for the compact letter display

box.rsltsp10 <- with(dat[c(vectorsize1),], graphics::boxplot(RGR ~ Sp., plot = FALSE))
box.rsltsp11 <- with(dat[c(vectorsize2),], graphics::boxplot(RGR ~ Sp., plot = FALSE))
box.rsltsp12 <- with(dat[c(vectorsize3),], graphics::boxplot(RGR ~ Sp., plot = FALSE))

# CLD ####

# CLD per sizes

xsp10 <- c(1:6)
ysp10 <- box.rsltsp10$stats[5, ]
cbdsp10 <- c("a","a","a","ab","ab","b")
ltr_dfsp10 <- data.frame(xsp10, ysp10, cbdsp10)

lmtssp10 <- get_plot_limits(pltsp10)
y.rangesp10 <- lmtssp10$ymax - lmtssp10$ymin
y.nudgesp10 <- 0.15 * y.rangesp10

xsp11 <- c(1:6)
ysp11 <- box.rsltsp11$stats[5, ]
cbdsp11 <- c("a","a","a","ab","ab","b")
ltr_dfsp11 <- data.frame(xsp11, ysp11, cbdsp11)

lmtssp11 <- get_plot_limits(pltsp11)
y.rangesp11 <- lmtssp11$ymax - lmtssp11$ymin
y.nudgesp11 <- 0.15 * y.rangesp11

xsp12 <- c(1:6)
ysp12 <- box.rsltsp12$stats[5, ]
cbdsp12 <- c("a","abc","a","ab","ad","abd")
ltr_dfsp12 <- data.frame(xsp12, ysp12, cbdsp12)

lmtssp12 <- get_plot_limits(pltsp12)
y.rangesp12 <- lmtssp12$ymax - lmtssp12$ymin
y.nudgesp12 <- 0.15 * y.rangesp12


# CLD alpha beta gamma ####

# We use unicode for greek letters to avoid conflicts and problems when revisiting the script

xABCsp10 <- c(0.6,1.6,2.6,3.6,4.6,5.6)
yABCsp10 <- c(0.3,0.3,0.3,0.3,0.3,0.3)
cABCsp10 <- c("\u03B1", "\u03B1","\u03B1","\u03B1","\u03B1","\u03B1")
ltr_dfABCsp10 <- data.frame(xABCsp10, yABCsp10, cABCsp10)

xABCsp11 <- c(0.6,1.6,2.6,3.6,4.6,5.6)
yABCsp11 <- c(0.25,0.25,0.25,0.25,0.25,0.25)
cABCsp11 <- c("\u03B1","\u03B1","\u03B2","\u03B1","\u03B2","\u03B2")
ltr_dfABCsp11 <- data.frame(xABCsp11, yABCsp11, cABCsp11)

xABCsp12 <- c(0.6,1.6,2.6,3.6,4.6,5.6)
yABCsp12 <- c(0.25,0.25,0.25,0.25,0.25,0.25)
cABCsp12 <- c("\u03B1","\u03B2","\u03B3","\u03B2","\u03B3","\u03B2")
ltr_dfABCsp12 <- data.frame(xABCsp12, yABCsp12, cABCsp12)


# CLD under whiskers ####

ydbj1 <- box.rsltsp1$stats[1, ]
ydbj2 <- box.rsltsp2$stats[1, ]
ydbj3 <- box.rsltsp3$stats[1, ] # c(0.000000, 2.246720, 0.000000, 1.405889, 0.000000, 3.8)
ydbj4 <- box.rsltsp4$stats[1, ]
ydbj5 <- box.rsltsp5$stats[1, ]
ydbj6 <- box.rsltsp6$stats[1, ]
ydbj7 <- box.rsltsp7$stats[1, ]
ydbj8 <- box.rsltsp8$stats[1, ]
ydbj9 <- box.rsltsp9$stats[1, ]
ydbj10 <- box.rsltsp10$stats[1, ]
ydbj11 <- box.rsltsp11$stats[1, ]
ydbj12 <- box.rsltsp12$stats[1, ]
y.nudgedbj1 <-  -0.05 * y.rangesp1
y.nudgedbj2 <-  -0.05 * y.rangesp2
y.nudgedbj3 <-  -0.05 * y.rangesp3
y.nudgedbj4 <-  -0.05 * y.rangesp4
y.nudgedbj5 <-  -0.05 * y.rangesp5
y.nudgedbj6 <-  -0.05 * y.rangesp6
y.nudgedbj7 <-  -0.05 * y.rangesp7
y.nudgedbj8 <-  -0.05 * y.rangesp8
y.nudgedbj9 <-  -0.05 * y.rangesp9
y.nudgedbj10 <-  -0.05 * y.rangesp10
y.nudgedbj11 <-  -0.05 * y.rangesp11
y.nudgedbj12 <-  -0.05 * y.rangesp12

### 2.6 - Plots #######

par(mfrow=c(1,1)) # to avoid problems later in the mosaic plot

# Number of established propagules ####

pltsp1 <- ggboxplot(dat[c(vectorsize1),],x= "Sp.", y= "Shoots_established_bc3", xlab=FALSE, fill = "Sp.", color = "black", palette = c('#777776FF','#762a83','#FFEA46FF','#A0DA39FF','#FC8D64FF','#08519c')) +
  rremove("legend") + 
  geom_text(data= ltr_dfABCsp1, aes(x=xABCsp1,y=yABCsp1, label=cABCsp1),colour = "black", size= 6) +
  geom_text(data = ltr_dfsp1, aes(x=xsp1, y= ydbj1, label=cbdsp1), colour = "black", nudge_y = y.nudgedbj1) + 
  ggtitle("Large") +
  ylab(" ") + 
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"),
        axis.text.x = element_text(size = 12, face = "bold.italic"), 
        plot.margin = margin(1, 0, 1, 0, "cm")) + 
  annotate("rect", fill= "black", colour= "grey",xmin = 0.45, xmax = 1.45, ymin = -0.75, ymax = 8.1, alpha= .06) +
  annotate("rect", fill= "black", colour= "grey",xmin = 2.45, xmax = 3.45, ymin = -0.75, ymax = 8.1, alpha= .06) +
  annotate("rect",fill= "black", colour= "grey", xmin = 4.45, xmax = 5.45, ymin = -0.75, ymax = 8.1, alpha= .06)

pltsp2 <- ggboxplot(dat[c(vectorsize2),],x= "Sp.", y= "Shoots_established_bc3", fill = "Sp.", color = "black", palette = c('#777776FF','#762a83','#FFEA46FF','#A0DA39FF','#FC8D64FF','#08519c')) +
  rremove("legend") + 
  geom_text(data= ltr_dfABCsp2, aes(x=xABCsp2,y=yABCsp2, label=cABCsp2), colour = "black", size= 5) + 
  geom_text(data = ltr_dfsp2, aes(x=xsp2, y= ydbj2, label=cbdsp2), colour = "black", nudge_y = y.nudgedbj2 ) +
  ggtitle("Medium") +
  ylab(" ") + 
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"),
        axis.text.x = element_text(size = 12, face = "bold.italic"), 
        plot.margin = margin(1, 0, 1, 0, "cm")) + 
  annotate("rect",fill= "black", colour= "grey", xmin = 0.45, xmax = 1.45, ymin = -0.75, ymax = 8.1, alpha= .06) +
  annotate("rect",fill= "black", colour= "grey", xmin = 2.45, xmax = 3.45, ymin = -0.75, ymax = 8.1, alpha= .06) +
  annotate("rect",fill= "black", colour= "grey", xmin = 4.45, xmax = 5.45, ymin = -0.75, ymax = 8.1, alpha= .06)

pltsp3 <- ggboxplot(dat[c(vectorsize3),],x= "Sp.", y= "Shoots_established_bc3", fill = "Sp.", color = "black", palette = c('#777776FF','#762a83','#FFEA46FF','#A0DA39FF','#FC8D64FF','#08519c')) + 
  rremove("legend") +
  geom_text(data= ltr_dfABCsp3, aes(x=xABCsp3,y=yABCsp3, label=cABCsp3), colour = "black", size= 5) + 
  geom_text(data = ltr_dfsp3, aes(x=xsp3, y= ydbj3, label=cbdsp3), colour = "black", nudge_y = y.nudgedbj3 ) + 
  ggtitle("Small") +
  ylab(" ") + 
  xlab("") + 
  theme(plot.title = element_text(hjust = 0.5, face= "bold"), 
        axis.text.x = element_text(size = 12, face = "bold.italic"), plot.margin = margin(1, 0, 1, 0, "cm")) +
  annotate("rect",fill= "black", colour= "grey", xmin = 0.45, xmax = 1.45, ymin = -0.75, ymax = 8.1, alpha= .06) +
  annotate("rect",fill= "black", colour= "grey", xmin = 2.45, xmax = 3.45, ymin = -0.75, ymax = 8.1, alpha= .06) +
  annotate("rect",fill= "black", colour= "grey", xmin = 4.45, xmax = 5.45, ymin = -0.75, ymax = 8.1, alpha= .06)

# Colonized surface ####

pltsp4 <- ggboxplot(dat[c(vectorsize1),],x= "Sp.", y= "perc_col_bc3", xlab=FALSE, fill = "Sp.", color = "black", palette = c('#777776FF','#762a83','#FFEA46FF','#A0DA39FF','#FC8D64FF','#08519c')) +
  rremove("legend") + 
  geom_text(data= ltr_dfABCsp4, aes(x=xABCsp4,y=yABCsp4, label=cABCsp4),colour = "black", size= 6) + 
  geom_text(data = ltr_dfsp4, aes(x=xsp4, y= ydbj4, label=cbdsp4), colour = "black", nudge_y = y.nudgedbj4) + 
  ggtitle(" ") +
  ylab(" ") + 
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"),
        axis.text.x = element_text(size = 12, face = "bold.italic"), 
        plot.margin = margin(1, 0, 1, 0, "cm")) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 0.45, xmax = 1.45, ymin = -2.3, ymax = 5, alpha= .06) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 2.45, xmax = 3.45, ymin = -2.3, ymax = 5, alpha= .06) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 4.45, xmax = 5.45, ymin = -2.3, ymax = 5, alpha= .06)


pltsp5 <- ggboxplot(dat[c(vectorsize2),],x= "Sp.", y= "perc_col_bc3", fill = "Sp.", color = "black", palette = c('#777776FF','#762a83','#FFEA46FF','#A0DA39FF','#FC8D64FF','#08519c')) +
  rremove("legend") + 
  geom_text(data= ltr_dfABCsp5, aes(x=xABCsp5,y=yABCsp5, label=cABCsp5),colour = "black", size= 6) + 
  geom_text(data = ltr_dfsp5, aes(x=xsp5, y= ydbj5, label=cbdsp5), colour = "black", nudge_y = y.nudgedbj5 ) +
  ggtitle(" ") +
  ylab(" ") + 
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"),
        axis.text.x = element_text(size = 12, face = "bold.italic"), 
        plot.margin = margin(1, 0, 1, 0, "cm")) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 0.45, xmax = 1.45, ymin = -2.3, ymax = 5, alpha= .06) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 2.45, xmax = 3.45, ymin = -2.3, ymax = 5, alpha= .06) +  
  annotate("rect", fill= "black", colour= "grey", xmin = 4.45, xmax = 5.45, ymin = -2.3, ymax = 5, alpha= .06)

pltsp6 <- ggboxplot(dat[c(vectorsize3),],x= "Sp.", y= "perc_col_bc3", fill = "Sp.", color = "black", palette = c('#777776FF','#762a83','#FFEA46FF','#A0DA39FF','#FC8D64FF','#08519c')) + 
  rremove("legend") +
  geom_text(data= ltr_dfABCsp6, aes(x=xABCsp6,y=yABCsp6, label=cABCsp6),colour = "black", size= 6) + 
  geom_text(data = ltr_dfsp6, aes(x=xsp6, y= ydbj6, label=cbdsp6), colour = "black", nudge_y = y.nudgedbj6 ) +
  ggtitle(" ") +
  ylab(" ") + 
  xlab("") + 
  theme(plot.title = element_text(hjust = 0.5, face= "bold"), 
        axis.text.x = element_text(size = 12, face = "bold.italic"), plot.margin = margin(1, 0, 1, 0, "cm")) +
  annotate("rect", fill= "black", colour= "grey", xmin = 0.45, xmax = 1.45, ymin = -2.3, ymax = 5, alpha= .06) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 2.45, xmax = 3.45, ymin = -2.3, ymax = 5, alpha= .06) +  
  annotate("rect", fill= "black", colour= "grey", xmin = 4.45, xmax = 5.45, ymin = -2.3, ymax = 5, alpha= .06)

# Viable biomass ####

pltsp7 <- ggboxplot(dat[c(vectorsize1),],x= "Sp.", y= "biomVper_bc3", xlab=FALSE, fill = "Sp.", color = "black", palette = c('#777776FF','#762a83','#FFEA46FF','#A0DA39FF','#FC8D64FF','#08519c')) +
  rremove("legend") + 
  geom_text(data= ltr_dfABCsp7, aes(x=xABCsp7,y=yABCsp7, label=cABCsp7),colour = "black", size= 6) + 
  geom_text(data = ltr_dfsp7, aes(x=xsp7, y= ydbj7, label=cbdsp7),colour = "black", nudge_y = y.nudgedbj7) +  
  ggtitle(" ") +
  ylab(" ") + 
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"),
        axis.text.x = element_text(size = 12, face = "bold.italic"), 
        plot.margin = margin(1, 0, 1, 0, "cm")) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 0.45, xmax = 1.45, ymin = 0.5, ymax = 3.4, alpha= .06) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 2.45, xmax = 3.45, ymin = 0.5, ymax = 3.4, alpha= .06) +  
  annotate("rect", fill= "black", colour= "grey", xmin = 4.45, xmax = 5.45, ymin = 0.5, ymax = 3.4, alpha= .06)

pltsp8 <- ggboxplot(dat[c(vectorsize2),],x= "Sp.", y= "biomVper_bc3", fill = "Sp.", color = "black", palette = c('#777776FF','#762a83','#FFEA46FF','#A0DA39FF','#FC8D64FF','#08519c')) +
  rremove("legend") + 
  geom_text(data= ltr_dfABCsp8, aes(x=xABCsp8,y=yABCsp8, label=cABCsp8),colour = "black", size= 6) + 
  geom_text(data = ltr_dfsp8, aes(x=xsp8, y= ydbj8, label=cbdsp8),colour = "black", nudge_y = y.nudgedbj8 ) + 
  ggtitle(" ") +
  ylab(" ") + 
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"),
        axis.text.x = element_text(size = 12, face = "bold.italic"), 
        plot.margin = margin(1, 0, 1, 0, "cm")) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 0.45, xmax = 1.45, ymin = 0.5, ymax = 3.4, alpha= .06) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 2.45, xmax = 3.45, ymin = 0.5, ymax = 3.4, alpha= .06) +  
  annotate("rect", fill= "black", colour= "grey", xmin = 4.45, xmax = 5.45, ymin = 0.5, ymax = 3.4, alpha= .06)

pltsp9 <- ggboxplot(dat[c(vectorsize3),],x= "Sp.", y= "biomVper_bc3", fill = "Sp.", color = "black", palette = c('#777776FF','#762a83','#FFEA46FF','#A0DA39FF','#FC8D64FF','#08519c')) + 
  rremove("legend") +
  geom_text(data= ltr_dfABCsp9, aes(x=xABCsp9,y=yABCsp9, label=cABCsp9),colour = "black", size= 6) + 
  geom_text(data = ltr_dfsp9, aes(x=xsp9, y= ydbj9, label=cbdsp9),colour = "black", nudge_y = y.nudgedbj9 ) + 
  ggtitle(" ") +
  ylab(" ") + 
  xlab("") + 
  theme(plot.title = element_text(hjust = 0.5, face= "bold"), 
        axis.text.x = element_text(size = 12, face = "bold.italic"), plot.margin = margin(1, 0, 1, 0, "cm")) +
  annotate("rect", fill= "black", colour= "grey", xmin = 0.45, xmax = 1.45, ymin = 0.5, ymax = 3.4, alpha= .06) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 2.45, xmax = 3.45, ymin = 0.5, ymax = 3.4, alpha= .06) +  
  annotate("rect", fill= "black", colour= "grey", xmin = 4.45, xmax = 5.45, ymin = 0.5, ymax = 3.4, alpha= .06)

# Relative Growth Rate ####

pltsp10 <- ggboxplot(dat[c(vectorsize1),],x= "Sp.", y= "RGR", xlab=FALSE, fill = "Sp.", color = "black", palette = c('#777776FF','#762a83','#FFEA46FF','#A0DA39FF','#FC8D64FF','#08519c')) +
  rremove("legend") + 
  geom_text(data= ltr_dfABCsp10, aes(x=xABCsp10,y=yABCsp10, label=cABCsp10),colour = "black", size= 6) + 
  geom_text(data = ltr_dfsp10, aes(x=xsp10, y= ydbj10, label=cbdsp10),colour = "black", nudge_y = y.nudgedbj10) +  
  ggtitle(" ") +
  ylab(" ") + 
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"),
        axis.text.x = element_text(size = 12, face = "bold.italic"), 
        plot.margin = margin(1, 0, 1, 0, "cm")) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 0.45, xmax = 1.45, ymin = -0.3, ymax = 0.35, alpha= .06) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 2.45, xmax = 3.45, ymin = -0.3, ymax = 0.35, alpha= .06) +  
  annotate("rect", fill= "black", colour= "grey", xmin = 4.45, xmax = 5.45, ymin = -0.3, ymax = 0.35, alpha= .06)

pltsp11 <- ggboxplot(dat[c(vectorsize2),],x= "Sp.", y= "RGR", fill = "Sp.", color = "black", palette = c('#777776FF','#762a83','#FFEA46FF','#A0DA39FF','#FC8D64FF','#08519c')) +
  rremove("legend") + 
  geom_text(data= ltr_dfABCsp11, aes(x=xABCsp11,y=yABCsp11, label=cABCsp11),colour = "black", size= 6) + 
  geom_text(data = ltr_dfsp11, aes(x=xsp11, y= ydbj11, label=cbdsp11),colour = "black", nudge_y = y.nudgedbj11) + 
  ggtitle(" ") +
  ylab(" ") + 
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5, face= "bold"),
        axis.text.x = element_text(size = 12, face = "bold.italic"), 
        plot.margin = margin(1, 0, 1, 0, "cm")) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 0.45, xmax = 1.45, ymin = -0.3, ymax = 0.3, alpha= .06) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 2.45, xmax = 3.45, ymin = -0.3, ymax = 0.3, alpha= .06) +  
  annotate("rect", fill= "black", colour= "grey", xmin = 4.45, xmax = 5.45, ymin = -0.3, ymax = 0.3, alpha= .06)

pltsp12 <- ggboxplot(dat[c(vectorsize3),],x= "Sp.", y= "RGR", fill = "Sp.", color = "black", palette = c('#777776FF','#762a83','#FFEA46FF','#A0DA39FF','#FC8D64FF','#08519c')) + 
  rremove("legend") +
  geom_text(data= ltr_dfABCsp12, aes(x=xABCsp12,y=yABCsp12, label=cABCsp12),colour = "black", size= 6) + 
  geom_text(data = ltr_dfsp12, aes(x=xsp12, y= ydbj12, label=cbdsp12),colour = "black", nudge_y = y.nudgedbj12) + 
  ggtitle(" ") +
  ylab(" ") + 
  xlab("") + 
  theme(plot.title = element_text(hjust = 0.5, face= "bold"), 
        axis.text.x = element_text(size = 12, face = "bold.italic"), 
        plot.margin = margin(1, 0, 1, 0, "cm")) +
  annotate("rect", fill= "black", colour= "grey", xmin = 0.45, xmax = 1.45, ymin = -0.3, ymax = 0.3, alpha= .06) + 
  annotate("rect", fill= "black", colour= "grey", xmin = 2.45, xmax = 3.45, ymin = -0.3, ymax = 0.3, alpha= .06) +  
  annotate("rect", fill= "black", colour= "grey", xmin = 4.45, xmax = 5.45, ymin = -0.3, ymax = 0.3, alpha= .06)

### 2.7 - Mosaic plot integration and export of the plot ####

# Correction of margins for a proper visualization of the complex mosaic plot. Margins adjusted for letting space to the common Y axis

pltsp1 <- pltsp1 + theme(plot.margin = unit(c(0, 0, 0, 0.25), "cm"))
pltsp2 <- pltsp2 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
pltsp3 <- pltsp3 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
pltsp4 <- pltsp4 + theme(plot.margin = unit(c(0, 0, 0, 0.25), "cm"))
pltsp5 <- pltsp5 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
pltsp6 <- pltsp6 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
pltsp7 <- pltsp7 + theme(plot.margin = unit(c(0, 0, 0, 0.25), "cm"))
pltsp8 <- pltsp8 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
pltsp9 <- pltsp9 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
pltsp10 <- pltsp10 + theme(plot.margin = unit(c(0, 0, 0, 0.25), "cm"))
pltsp11 <- pltsp11 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
pltsp12 <- pltsp12 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

# Export

png(file="myfig1.png",width=1468,height=1117, units = "px", res = 100) # it will create the png file on the current directory. It helps to improve resolution and reproducibility when saving the file.
plot_grid(pltsp1, pltsp2, pltsp3,
          pltsp4, pltsp5,pltsp6,
          pltsp7, pltsp8, pltsp9,
          pltsp10, pltsp11, pltsp12, ncol=3, nrow = 4, align = "hv")
grid.text(label = "Number of established shoots", x= 0.01 , y=0.87, rot = 90, gp= gpar(fontsize=14, fontface = "bold"))
grid.text(label = "Colonized surface", x= 0.01 , y=0.625, rot = 90, gp= gpar(fontsize=14, fontface = "bold"))
grid.text(label = "Viable biomass", x= 0.01 , y=0.385, rot = 90, gp= gpar(fontsize=14, fontface = "bold"))
grid.text(label = "RGR (g/month)", x= 0.01 , y=0.125, rot = 90, gp= gpar(fontsize=14, fontface = "bold"))

dev.off() # once device is set off, the plot with the measurements specified is exported correctly to the current directory
 
# Cleaning for next part of the code ####

rm(list = ls()) # to clean data before start

############ 3 - Correlations of indicators with traits (Figure 2) -------------------------------------------------------------------




### 3.1 - General preparation -----------------------------------------------
# Comments ----------------------------------------------------------------

# This figure incorporates many difficulties. One of them is working with different datasets.
# The first part of this section is related to merge the dataframes and using them for testing and plotting correlations 
# Second, we have to export the subplots of the correlations and the network plot separately and call them as images and integrate them in a mosaic plot,
# as the result corplot is a class list and trying to change it to any kind of 
# plot class affects the quality of the mosaic plot.


# Libraries ####

library("ggpubr") # for plotting with the modification of network_plot base on ggplot (v. 0.4.0)
library(corrr) # for the network plot using correlate (v. 0.4.4)
library(dplyr) # for managing the data and the pipelines (v. 2.2.1)
library(psych) # for corr.test that includes adjust function (v. 2.2.5)
library(corrplot) # for the corrplot (v. 0.92)
library(gridExtra) # for the mosaic plot with marrangeGrob (v. 2.3)
library(grid) # for the annotations in the figures (v. 4.1.2)
library(grDevices) # for using png for saving the images (v. 4.1.2)
library(car) # for Levene's test, Anova and boxCox (v.3.1-0)
library(stats) # for shapiro.test, ks.test, bartlett.test, fligner.test, aov for doing the models and preparing the transformation, residuals for getting the overdispersion, cor for the R2 (v. 4.1.2)

# Dataset 1 ####

# Size and shape traits measurements

path.to.data <- "https://raw.githubusercontent.com/Fernando-Hurtado/DCR/main/Establishment/"
dat <- read.csv(paste0(path.to.data,"Dataset_Traits.csv"), header=T, sep=";", dec= ",")

# Dataset 2 ####

# Indicators of establishment success measurements

path.to.data2 <- "https://raw.githubusercontent.com/Fernando-Hurtado/DCR/main/Establishment/"
dat2 <- read.csv(paste0(path.to.data2,"Dataset_Indicators.csv"), header=T, sep=";", dec= ",")

## Incorporation of the transformed indicators and RGR (non-transformed) -----------------------------

# Number of established propagules (propagules/mm2) ------------------------------------------------------

options(contrasts=c(factor="contr.sum", ordered="contr.poly")) ## We have factors in the models, so we change the contrast matrix as it do not corresponds with the theory in stadistics, by deafult R applies a different one
eqt3 <- as.formula(Shoots_established ~ Sp.*Size)
Shoots_established_aov3 <- aov(eqt3, data= dat2)
Shoots_established <- dat2$Shoots_established + min(dat2$Shoots_established[dat2$Shoots_established!=0], na.rm = TRUE) # Box-cox transformations do not deal well with 0 values, so we prepare the data for the transformation
dat2$Shoots_established <- Shoots_established
lambda.bc3 <- with(boxCox(Shoots_established_aov3, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
Shoots_established_bc3 <- ((Shoots_established^lambda.bc3)-1)/(lambda.bc3)
dat2$Shoots_established_bc3 <- Shoots_established_bc3

# Colonized surface -------------------------------------------------------

options(contrasts=c(factor="contr.sum", ordered="contr.poly")) ## We have factors in the models, so we change the contrast matrix as it do not corresponds with the theory in stadistics, by deafult R applies a different one
dat2$perc_col <- dat2$S_Shootss/dat2$S_substratum
dat2$perc_col <- dat2$perc_col*100
eqt3 <- as.formula(perc_col ~ Sp.*Size)
perc_col_aov3 <- aov(eqt3, data= dat2)
perc_col <- dat2$perc_col + min(dat2$perc_col[dat2$perc_col!=0], na.rm = TRUE) # Box-cox transformations do not deal well with 0 values, so we prepare the data for the transformation
dat2$perc_col <- perc_col
lambda.bc3 <- with(boxCox(perc_col_aov3, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
perc_col_bc3 <- ((perc_col^lambda.bc3)-1)/(lambda.bc3)
dat2$perc_col_bc3 <- perc_col_bc3

# Viable biomass ----------------------------------------------------------

options(contrasts=c(factor="contr.sum", ordered="contr.poly")) ## We have factors in the models, so we change the contrast matrix as it do not corresponds with the theory in stadistics, by deafult R applies a different one
dat2$biomVper <- dat2$per_BiomV
eqt3 <- as.formula(biomVper ~ Sp.*Size)
biomVper_aov3 <- aov(eqt3, data= dat2)
biomVper <- dat2$biomVper + min(dat2$biomVper[dat2$biomVper!=0], na.rm = TRUE) # Box-cox transformations do not deal well with 0 values, so we prepare the data for the transformation
dat2$biomVper <- biomVper
lambda.bc3 <- with(boxCox(biomVper_aov3, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
biomVper_bc3 <- ((biomVper^lambda.bc3)-1)/(lambda.bc3)
dat2$biomVper_bc3 <- biomVper_bc3

# RGR ---------------------------------------------------------------------

# We create the variable RGR. It does not need transformation

dat2$RGR <- (log(dat2$Biom_total) - log(0.0040))/2

# Create the function for the merge ####

# This function is really useful, we use "Sample" as a common field for the merge for managing the datasets

MyMerge <- function(x, y){
  df <- merge(x, y, by= "Sample", all=TRUE, sort=FALSE)
  return(df)
}

# Modification of network_plot function -----------------------------------

# Modification of network_plot function to get dark points and customizable range for legend. Basically I added from source the legend range and I set dark points to default in-built instead of white for a better contrast

network_plot.cor_df_black<- function(rdf,
                                     min_cor = .30,
                                     legend = c("full", "range", "none"),
                                     colours = c("indianred2", "white", "skyblue1"),
                                     repel = TRUE,
                                     curved = TRUE,
                                     legend_range=c(-1,1),
                                     colors) {
  legend <- rlang::arg_match(legend)
  
  if (min_cor < 0 || min_cor > 1) {
    rlang::abort("min_cor must be a value ranging from zero to one.")
  }
  
  if (!missing(colors)) {
    colours <- colors
  }
  
  rdf <- as_matrix(rdf, diagonal = 1)
  distance <- 1 - abs(rdf)
  
  points <- if (ncol(rdf) == 1) {
    # 1 var: a single central point
    matrix(c(0, 0), ncol = 2, dimnames = list(colnames(rdf)))
  } else if (ncol(rdf) == 2) {
    # 2 vars: 2 opposing points
    matrix(c(0, -0.1, 0, 0.1), ncol = 2, dimnames = list(colnames(rdf)))
  } else {
    # More than 2 vars: multidimensional scaling to obtain x and y coordinates for points.
    suppressWarnings(stats::cmdscale(distance, k = 2))
  }
  
  if (ncol(points) < 2) {
    cont_flag <- FALSE
    shift_matrix <- matrix(1,
                           nrow = nrow(rdf),
                           ncol = ncol(rdf)
    )
    diag(shift_matrix) <- 0
    
    for (shift in 10^(-6:-1)) {
      shifted_distance <- distance + shift * shift_matrix
      points <- suppressWarnings(stats::cmdscale(shifted_distance))
      
      if (ncol(points) > 1) {
        cont_flag <- TRUE
        break
      }
    }
    
    if (!cont_flag) rlang::abort("Can't generate network plot.\nAttempts to generate 2-d coordinates failed.")
    
    rlang::warn("Plot coordinates derived from correlation matrix have dimension < 2.\nPairwise distances have been adjusted to facilitate plotting.")
  }
  
  
  points <- data.frame(points)
  colnames(points) <- c("x", "y")
  points$id <- rownames(points)
  
  # Create a proximity matrix of the paths to be plotted.
  proximity <- abs(rdf)
  proximity[upper.tri(proximity)] <- NA
  diag(proximity) <- NA
  proximity[proximity < min_cor] <- NA
  
  # Produce a data frame of data needed for plotting the paths.
  n_paths <- sum(!is.na(proximity))
  paths <- data.frame(matrix(nrow = n_paths, ncol = 6))
  colnames(paths) <- c("x", "y", "xend", "yend", "proximity", "sign")
  path <- 1
  for (row in 1:nrow(proximity)) {
    for (col in 1:ncol(proximity)) {
      path_proximity <- proximity[row, col]
      if (!is.na(path_proximity)) {
        path_sign <- sign(rdf[row, col])
        x <- points$x[row]
        y <- points$y[row]
        xend <- points$x[col]
        yend <- points$y[col]
        paths[path, ] <- c(x, y, xend, yend, path_proximity, path_sign)
        path <- path + 1
      }
    }
  }
  
  if(legend %in% c("full", "none")){
    legend_range = legend_range
  }
  else if(legend == "range"){
    legend_range = c(min(rdf[row(rdf)!=col(rdf)]),
                     max(rdf[row(rdf)!=col(rdf)]))
  }
  plot_ <- list(
    # For plotting paths
    if (curved) {
      geom_curve(
        data = paths,
        aes(
          x = x, y = y, xend = xend, yend = yend,
          alpha = proximity, size = proximity,
          colour = proximity * sign
        )
      )
    },
    if (!curved) {
      geom_segment(
        data = paths,
        aes(
          x = x, y = y, xend = xend, yend = yend,
          alpha = proximity, size = proximity,
          colour = proximity * sign
        )
      )
    },
    scale_alpha(limits = c(0, 1)),
    scale_size(limits = c(0, 1)),
    scale_colour_gradientn(limits = legend_range, colors = colours),
    # Plot the points
    geom_point(
      data = points,
      aes(x, y),
      size = 4, shape = 19, colour = "black"
    ),
    # Plot variable labels
    if (repel) {
      ggrepel::geom_text_repel(
        data = points,
        aes(x, y, label = id),
        fontface = "bold", size = 5,
        segment.size = 0.0,
        segment.color = "white"
      )
    },
    if (!repel) {
      geom_text(
        data = points,
        aes(x, y, label = id),
        fontface = "bold", size = 5
      )
    },
    # expand the axes to add space for curves
    expand_limits(
      x = c(
        min(points$x) - .1,
        max(points$x) + .1
      ),
      y = c(
        min(points$y) - .1,
        max(points$y) + .1
      )
    ),
    # Theme and legends
    theme_void(),
    guides(size = "none", alpha = "none"),
    if (legend != "none") labs(colour = NULL),
    if (legend == "none") theme(legend.position = "none")
  )
  
  ggplot() + plot_
}




# Creation of Area diff (wet vs. dry area) and adjust to mm2 ####

dat$Area_wet <- (dat$Area_wet)*100
dat$Area_dry <- (dat$Area_dry)*100
dat$Area_diff <- (dat$Area_wet - dat$Area_dry)/dat$Area_wet

# Creation of Feret Diameter diff (wet vs. dry length) and adjust to mm ####

dat$Feret_wet <- (dat$Feret_wet)*10
dat$Feret_dry <- (dat$Feret_dry)*10
dat$Feret_diff <- (dat$Feret_wet - dat$Feret_dry)/dat$Feret_wet

# Creation of new datasets with mean for correlations ####

# Shoots %

Shoots_mean <- as.data.frame(dat %>%
                               filter(Type == "S") %>%
                               count(Sample,Type == "S"))
Shoots_mean[,3]<-(Shoots_mean[,3]/30)*100
Shoots_mean<-Shoots_mean[,-2]
names(Shoots_mean)[names(Shoots_mean) == "n"] <- "Shoots_mean"

# Apparent viability

Index_viability <- as.data.frame(dat %>%
                                   group_by(Sample)%>%
                                   summarise(Index_viability=mean(Index.of.viability,na.rm=TRUE),
                                             Index_viability_Sd=sd(Index.of.viability,na.rm=TRUE)))

# Area (mm?)

Area_prop <- as.data.frame(dat %>%
                             group_by(Sample)%>%
                             summarise(Area_prop=mean(Area_wet,na.rm=TRUE),
                                       Area_prop_Sd=sd(Area_wet,na.rm=TRUE)))

# Length (mm)

Maximum_caliper <- as.data.frame(dat %>%
                                   group_by(Sample)%>%
                                   summarise(Maximum_caliper=mean(Feret_wet,na.rm=TRUE),
                                             Maximum_caliper_Sd=sd(Feret_wet,na.rm=TRUE)))

# Circularity

Circ_prop <- as.data.frame(dat %>%
                             group_by(Sample)%>%
                             summarise(Circ_prop=mean(Circ._wet,na.rm=TRUE),
                                       Circ_prop_Sd=sd(Circ._wet,na.rm=TRUE)))

# Wet vs. dry area

Area_diff_prop <- as.data.frame(dat %>%
                                  group_by(Sample)%>%
                                  summarise(Area_diff_prop=mean(Area_diff,na.rm=TRUE),
                                            Area_diff_prop_Sd=sd(Area_diff,na.rm=TRUE)))

# Wet vs. dry length

FeretD_diff <- as.data.frame(dat %>%
                               group_by(Sample)%>%
                               summarise(FeretD_diff=mean(Feret_diff,na.rm=TRUE),
                                         FeretD_diff_Sd=sd(Feret_diff,na.rm=TRUE)))

# Number of established propagules

Prop_established <- as.data.frame(dat2 %>%
                                    group_by(Sample)%>%
                                    summarise(Prop_established=mean(Shoots_established_bc3,na.rm=TRUE),
                                              Prop_established_Sd=sd(Shoots_established_bc3,na.rm=TRUE)))

# Colonized surface

Area_exp <- as.data.frame(dat2 %>%
                            group_by(Sample)%>%
                            summarise(Area_exp=mean(perc_col_bc3,na.rm=TRUE),
                                      Area_exp_Sd=sd(perc_col_bc3,na.rm=TRUE)))

# Viable biomass

Perc_Bviable <- as.data.frame(dat2 %>%
                                group_by(Sample)%>%
                                summarise(Perc_Bviable=mean(biomVper_bc3,na.rm=TRUE),
                                          Perc_Bviable_Sd=sd(biomVper_bc3,na.rm=TRUE)))

# RGR (g/month)

Growth_biomassdiff <- as.data.frame(dat2 %>%
                                      group_by(Sample)%>%
                                      summarise(Growth_biomassdiff=mean(RGR,na.rm=TRUE),
                                                Growth_biomassdiff_Sd=sd(RGR,na.rm=TRUE)))

# We reduce here using the function created for merging

dat3 <- Reduce(MyMerge, list(Shoots_mean, Index_viability,Area_prop,Maximum_caliper,Circ_prop,Area_diff_prop,FeretD_diff,
                             Prop_established,Area_exp,Perc_Bviable,Growth_biomassdiff))
dat3 <- as.data.frame(dat3) # we convert it in a data frame
dat4 <- dat3[,c('Shoots_mean', 'Index_viability','Area_prop','Maximum_caliper','Circ_prop','Area_diff_prop','FeretD_diff',
                'Prop_established','Area_exp','Perc_Bviable','Growth_biomassdiff')]

dat4names <- c("Shoot percentage", "Apparent viability","Area (mm?)","Length (mm)","Circularity","Wet vs dry area","Wet vs dry length",
               "Number of established shoots","Colonized surface","Viable biomass","RGR (g/month)")

rownames(dat4) <- dat3$Sample
colnames(dat4) <- dat4names


### 3.2 - Correlations ------------------------------------------------------------

# We separate here now that all the traits and indicators are paired, so we can confront them testing the correlations

dat5 <- dat4[,c("Shoot percentage", "Apparent viability","Area (mm?)","Length (mm)","Circularity","Wet vs dry area","Wet vs dry length")]
dat6 <- dat4[,c("Number of established shoots","Colonized surface","Viable biomass","RGR (g/month)")]

# We go now for the values of the correlations using the p.adjusted values for assessing the plots
# The p.adjusted values are obtained throw hypothesis driven analysis, going iteratively over the the different hypothesis
# Loops for the p.adjusted values (holm) ------------------------------------------------------

Testdf<-data.frame(x=character(),y=character(),Tholm1=character(),Tholm2=character(),stringsAsFactors=FALSE)

for (i in 1:ncol(dat5)){
  for(j in 1:ncol(dat6)){
    eqt4 <- colnames(dat5[i])
    eqt5 <- colnames(dat6[j])
    Tholm <- corr.test(dat5[,eqt4], dat6[,eqt5], method = "pearson", adjust = "holm", alpha = .05, ci= TRUE)
    testdfi<-c(eqt4,eqt5,as.numeric(Tholm$p.adj),as.numeric(Tholm$r))
    Testdf<-rbind(Testdf,testdfi,stringsAsFactors=FALSE)
    }}

colnames(Testdf)<-c("Traits","Indicators","p.adjusted","r")

Testdf$p.adjusted <- as.numeric(Testdf$p.adjusted)
Testdf$r <- as.numeric(Testdf$r)


# Creation of the matrix for the plots and the matrix for assessing the plots --------

## set up storage matrices, one for the p.adjusted values and other for the r values
# get names for row and columns

nameVals <- dat4names

# construct 0 matrices of correct dimensions with row and column names

myMattestA <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
myMattestB <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))

# fill in the matrices with matrix indexing on row and column names

myMattestA[as.matrix(Testdf[c("Traits", "Indicators")])] <- Testdf[["p.adjusted"]]
myMattestB[as.matrix(Testdf[c("Traits", "Indicators")])] <- Testdf[["r"]]

# Transform in square matrices for the plots

myMattestA.t <- t(myMattestA)

myMattestA2 <- myMattestA + myMattestA.t

myMattestB.t <- t(myMattestB)

# we will use this next matrix with the original r values for both the networkplot and the corrplot, 
# it will represent better the non-significant traits in the clustering of the network plot, given a better global picture and we can use the p.adjust matrix for 2 reasons:
# setting the plotting conditions of the network plot for a better representation of the statistically significant correlations and also assessing the significance codes for the corrplot

myMattestB2 <- myMattestB + myMattestB.t 

# Using the p.adjusted values for selecting the r values of the significant correlations. We overwrite the r values using conditions using 10 and 0
# Using 10 for erasing laterly non significant r values

NewTestdfA <- as.data.frame(myMattestA2)
NewTestdfA[NewTestdfA > 0.05] <- 10
NewTestdfA[NewTestdfA == 0] <- 10

# Using 0 for getting directly the significant r values later

NewTestdfA[NewTestdfA < 0.05] <- 0

NewTestdfA <- as.matrix(NewTestdfA)

# We add the values to select easier with the values = 10 and = 0

NewTestdfB <- myMattestB2 + NewTestdfA
NewTestdfB <- as.data.frame(NewTestdfB)

# We finally get the selection of the significant r values, erasing the ones non-significant and non altering the original values

NewTestdfB[NewTestdfB > 1] <- 0

# Using the information of this NewTestdfb dataframe information we can now create the annotations for the corrplot when plotting and set the plotting conditions for the networkplot

# We prepare the matrices for the plots

dat7 <- as.data.frame(myMattestB2)
dat7 <- dat7[-c(8:11),-c(1:7)]  # we simplify for the proper usage of the corrplot
dat7 <- as.matrix(dat7) 

# We use the original one with the new plotting conditions for using the best of each matrix, appearing only the significant correlations with also a better perspective of the non significant traits

dat8 <- as_cordf(as.data.frame(myMattestB2)) # coercing to correlation dataframe for the networkplot. 


### 3.3 - Plots -------------------------------------------------------------------
# Plot Fig 2A and export of the plot ####

col <- colorRampPalette(c('#2166ac','#67a9cf','#d1e5f0','#fddbc7','#ef8a62','#b2182b')) #colorblindsafe
png(file="myfig2a.png",width=2200,height=2700, units = "px", res = 375) # it will create the png file on the current directory. It helps to improve resolution and reproducibility when saving the file.
par(xpd = TRUE) # allow to plot in the margins
windowsFonts() # to avoid problems with the fonts of the corrplot
  corrplot(dat7,method="color",  col=col(200), cl.pos = 'r',  tl.pos = "full",
           addCoef.col = "white", # Add coefficient of correlation
           tl.col="black", tl.srt=45, is.corr = FALSE, family="sans" ,number.digits = 3, 
           na.label = "square", na.label.col = "white",tl.cex = 1.5,number.cex = 1.3, cl.cex = 1.3, mar = c(0,0,12,2),cl.align.text = "l")
# Annotations on the plot, significance codes based on NewTestdfB. They will be saved if we do dev.off later
  text(x=c(rep(1,7)),y=c(1.3,2.3,3.3,4.3,5.3,6.3,7.3), c("**","*","","","","",""),col = "white", cex = 1.3)+
  text(x=c(rep(2,7)),y=c(1.3,2.3,3.3,4.3,5.3,6.3,7.3), c("*","*","","","","",""),col = "white", cex = 1.3)+
  text(x=c(rep(3,7)),y=c(1.3,2.3,3.3,4.3,5.3,6.3,7.3), c("","","","","","","*"),col = "white", cex = 1.3)+
  text(x=c(rep(4,7)),y=c(1.3,2.3,3.3,4.3,5.3,6.3,7.3), c("*","*","","**","","*",""),col = "white", cex = 1.3)+
  
dev.off() # once device is set off, the plot with the measurements specified is exported correctly to the current directory

# Plot Fig2B and export of the plot ####

col <- colorRampPalette(c('#fddbc7','#ef8a62','#b2182b')) #colorblindsafe. We use the positive section of the previous palette, as the selection is only positive to be coherent with the previous plot 2a
png(file="myfig2b.png",width=2100,height=1800, units = "px", res = 250) # it will create the png file on the current directory. It helps to improve resolution and reproducibility when saving the file.
         network_plot.cor_df_black(dat8,
                          min_cor = 0.47,
                          colours=  col(200),# c('#0571b0','#92c5de','#f4a582','#ca0020'),
                          legend_range =c(0.4,0.7),
                          repel = TRUE)+ 
                          theme(legend.text = element_text(colour="black", size=12.5,face="bold"), plot.margin = margin(0,0.15,0,0.15, "cm")) # to prevent legend to be without margin
         
dev.off() # once device is set off, the plot with the measurements specified is exported correctly to the current directory

### 3.4 - Mosaic plot and export -------------------------------------------------------------

# In this case the mosaic plot is more difficult to assamble, as the corplot resulted is a list, so we have to call back the images as png
rl <- lapply(list("myfig2a.png", "myfig2b.png"), png::readPNG)
gl <- lapply(rl, grid::rasterGrob)

# Export

png(file="myfig2.png",width=1800,height=900, units = "px", res = 100) # it will create the png file on the current directory. It helps to improve resolution and reproducibility when saving the file.
  par(xpd = TRUE) # allow to plot in the margins
  par(mar = c(0,0,0,0), xaxs="i", yaxs="i") # another adjust of the plotting conditions for better reproducibility
  marrangeGrob(gl,nrow = 1, ncol = 2, top=NULL) # to display both plots in one line
  
# Annotations on the plot, they will be saved if we do dev.off later
  
  grid.text(label = "A", x= 0.075 , y=0.925, rot = 0, gp= gpar(fontsize=24, fontface = "bold"))
  grid.text(label = "B", x= 0.5 , y=0.925, rot = 0, gp= gpar(fontsize=24, fontface = "bold"))

dev.off() # once device is set off, the plot with the measurements specified is exported correctly to the current directory

# Cleaning for next part of the code ####

rm(list = ls()) # to clean data before start

############ 4 - Correlations of indicators with traits 2 (Figure S1) --------------------

### 4.1 - General preparation -----------------------------------------------
# Comments ----------------------------------------------------------------

# There is a big part of code which is the same of the figure 2. There are some extra packages for the plots

# Libraries ####

library("ggpubr") # for the use of ggscatter (v. 0.4.0)
library(dplyr) # for managing the data and the pipelines (v. 2.2.1)
library(psych) # for corr.test that includes adjust function (v. 2.2.5)
library(grid) # for the annotations in the figures (v. 4.1.2)
library(cowplot) # for the mosaic plot (v. 1.1.1)
library(ggtext) # for the italics in the axis labels of the plots (v. 0.1.2)
library(grDevices) # for using png for saving the images (v. 4.1.2)
library(car) # for Levene's test, Anova and boxCox (v.3.1-0)
library(stats) # for shapiro.test, ks.test, bartlett.test, fligner.test, aov for doing the models and preparing the transformation, residuals for getting the overdispersion, cor for the R2 (v. 4.1.2)

### 4.2 - Same parts of code of figure 2 ------------------------------------------


# Dataset 1 ####

# Size and shape traits measurements

path.to.data <- "https://raw.githubusercontent.com/Fernando-Hurtado/DCR/main/Establishment/"
dat <- read.csv(paste0(path.to.data,"Dataset_Traits.csv"), header=T, sep=";", dec= ",")

# Dataset 2 ####

# Indicators of establishment success measurements

path.to.data2 <- "https://raw.githubusercontent.com/Fernando-Hurtado/DCR/main/Establishment/"
dat2 <- read.csv(paste0(path.to.data2,"Dataset_Indicators.csv"), header=T, sep=";", dec= ",")

## Incorporation of the transformed indicators and RGR (non-transformed) -----------------------------

# Number of established propagules (propagules/mm2) ------------------------------------------------------

options(contrasts=c(factor="contr.sum", ordered="contr.poly")) ## We have factors in the models, so we change the contrast matrix as it do not corresponds with the theory in stadistics, by deafult R applies a different one
eqt3 <- as.formula(Shoots_established ~ Sp.*Size)
Shoots_established_aov3 <- aov(eqt3, data= dat2)
Shoots_established <- dat2$Shoots_established + min(dat2$Shoots_established[dat2$Shoots_established!=0], na.rm = TRUE) # Box-cox transformations do not deal well with 0 values, so we prepare the data for the transformation
dat2$Shoots_established <- Shoots_established
lambda.bc3 <- with(boxCox(Shoots_established_aov3, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
Shoots_established_bc3 <- ((Shoots_established^lambda.bc3)-1)/(lambda.bc3)
dat2$Shoots_established_bc3 <- Shoots_established_bc3

# Colonized surface -------------------------------------------------------

options(contrasts=c(factor="contr.sum", ordered="contr.poly")) ## We have factors in the models, so we change the contrast matrix as it do not corresponds with the theory in stadistics, by deafult R applies a different one
dat2$perc_col <- dat2$S_Shootss/dat2$S_substratum
dat2$perc_col <- dat2$perc_col*100
eqt3 <- as.formula(perc_col ~ Sp.*Size)
perc_col_aov3 <- aov(eqt3, data= dat2)
perc_col <- dat2$perc_col + min(dat2$perc_col[dat2$perc_col!=0], na.rm = TRUE) # Box-cox transformations do not deal well with 0 values, so we prepare the data for the transformation
dat2$perc_col <- perc_col
lambda.bc3 <- with(boxCox(perc_col_aov3, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
perc_col_bc3 <- ((perc_col^lambda.bc3)-1)/(lambda.bc3)
dat2$perc_col_bc3 <- perc_col_bc3

# Viable biomass ----------------------------------------------------------

options(contrasts=c(factor="contr.sum", ordered="contr.poly")) ## We have factors in the models, so we change the contrast matrix as it do not corresponds with the theory in stadistics, by deafult R applies a different one
dat2$biomVper <- dat2$per_BiomV
eqt3 <- as.formula(biomVper ~ Sp.*Size)
biomVper_aov3 <- aov(eqt3, data= dat2)
biomVper <- dat2$biomVper + min(dat2$biomVper[dat2$biomVper!=0], na.rm = TRUE) # Box-cox transformations do not deal well with 0 values, so we prepare the data for the transformation
dat2$biomVper <- biomVper
lambda.bc3 <- with(boxCox(biomVper_aov3, lambda=seq(-2,2, 1/1000)), x[which.max(y)])
biomVper_bc3 <- ((biomVper^lambda.bc3)-1)/(lambda.bc3)
dat2$biomVper_bc3 <- biomVper_bc3

# RGR ---------------------------------------------------------------------

# We create the variable RGR. It does not need transformation

dat2$RGR <- (log(dat2$Biom_total) - log(0.0040))/2

# Create the function for the merge ####

# This function is really useful, we use "Sample" as a common field for the merge for managing the datasets

MyMerge <- function(x, y){
  df <- merge(x, y, by= "Sample", all=TRUE, sort=FALSE)
  return(df)
}

# Obtaining RGR values ####

# being included and not needing to transform, we apply the formula directly

dat2$RGR <- (log(dat2$Biom_total) - log(0.0040))/2

# Creation of Area diff (wet vs. dry area) and adjust to mm2 ####

dat$Area_wet <- (dat$Area_wet)*100
dat$Area_dry <- (dat$Area_dry)*100
dat$Area_diff <- (dat$Area_wet - dat$Area_dry)/dat$Area_wet

# Creation of Feret Diameter diff (wet vs. dry length) and adjust to mm ####

dat$Feret_wet <- (dat$Feret_wet)*10
dat$Feret_dry <- (dat$Feret_dry)*10
dat$Feret_diff <- (dat$Feret_wet - dat$Feret_dry)/dat$Feret_wet

# Creation of new datasets with mean for correlations ####

# Shoots %

Shoots_mean <- as.data.frame(dat %>%
                               filter(Type == "S") %>%
                               count(Sample,Type == "S"))
Shoots_mean[,3]<-(Shoots_mean[,3]/30)*100
Shoots_mean<-Shoots_mean[,-2]
names(Shoots_mean)[names(Shoots_mean) == "n"] <- "Shoots_mean"

# Apparent viability

Index_viability <- as.data.frame(dat %>%
                                   group_by(Sample)%>%
                                   summarise(Index_viability=mean(Index.of.viability,na.rm=TRUE),
                                             Index_viability_Sd=sd(Index.of.viability,na.rm=TRUE)))

# Area (mm?)

Area_prop <- as.data.frame(dat %>%
                             group_by(Sample)%>%
                             summarise(Area_prop=mean(Area_wet,na.rm=TRUE),
                                       Area_prop_Sd=sd(Area_wet,na.rm=TRUE)))

# Length (mm)

Maximum_caliper <- as.data.frame(dat %>%
                                   group_by(Sample)%>%
                                   summarise(Maximum_caliper=mean(Feret_wet,na.rm=TRUE),
                                             Maximum_caliper_Sd=sd(Feret_wet,na.rm=TRUE)))

# Circularity

Circ_prop <- as.data.frame(dat %>%
                             group_by(Sample)%>%
                             summarise(Circ_prop=mean(Circ._wet,na.rm=TRUE),
                                       Circ_prop_Sd=sd(Circ._wet,na.rm=TRUE)))

# Wet vs. dry area

Area_diff_prop <- as.data.frame(dat %>%
                                  group_by(Sample)%>%
                                  summarise(Area_diff_prop=mean(Area_diff,na.rm=TRUE),
                                            Area_diff_prop_Sd=sd(Area_diff,na.rm=TRUE)))

# Wet vs. dry length

FeretD_diff <- as.data.frame(dat %>%
                               group_by(Sample)%>%
                               summarise(FeretD_diff=mean(Feret_diff,na.rm=TRUE),
                                         FeretD_diff_Sd=sd(Feret_diff,na.rm=TRUE)))

# Number of established propagules

Prop_established <- as.data.frame(dat2 %>%
                                    group_by(Sample)%>%
                                    summarise(Prop_established=mean(Shoots_established_bc3,na.rm=TRUE),
                                              Prop_established_Sd=sd(Shoots_established_bc3,na.rm=TRUE)))

# Colonized surface

Area_exp <- as.data.frame(dat2 %>%
                            group_by(Sample)%>%
                            summarise(Area_exp=mean(perc_col_bc3,na.rm=TRUE),
                                      Area_exp_Sd=sd(perc_col_bc3,na.rm=TRUE)))

# Viable biomass

Perc_Bviable <- as.data.frame(dat2 %>%
                                group_by(Sample)%>%
                                summarise(Perc_Bviable=mean(biomVper_bc3,na.rm=TRUE),
                                          Perc_Bviable_Sd=sd(biomVper_bc3,na.rm=TRUE)))

# RGR (g/month)

Growth_biomassdiff <- as.data.frame(dat2 %>%
                                      group_by(Sample)%>%
                                      summarise(Growth_biomassdiff=mean(RGR,na.rm=TRUE),
                                                Growth_biomassdiff_Sd=sd(RGR,na.rm=TRUE)))

# We reduce here using the function created for merging

dat3 <- Reduce(MyMerge, list(Shoots_mean, Index_viability,Area_prop,Maximum_caliper,Circ_prop,Area_diff_prop,FeretD_diff,
                             Prop_established,Area_exp,Perc_Bviable,Growth_biomassdiff))
dat3 <- as.data.frame(dat3) # we convert it in a data frame
dat4 <- dat3[,c('Shoots_mean', 'Index_viability','Area_prop','Maximum_caliper','Circ_prop','Area_diff_prop','FeretD_diff',
                'Prop_established','Area_exp','Perc_Bviable','Growth_biomassdiff')]


### 4.3 - Different part of code: final Scatterplots ####

# 1 Selection of final plots

correplot1 <- ggscatter(dat4, x = "Perc_Bviable", y = "Shoots_mean", 
                        add = "reg.line", conf.int = TRUE, 
                        xlab = "Viable biomass ", ylab = "Shoot percentage")

correplot2 <- ggscatter(dat4, x = "Growth_biomassdiff", y = "Index_viability", 
                        add = "reg.line", conf.int = TRUE, 
                        xlab = "RGR (g/month) ", ylab = "Apparent viability")

correplot3 <- ggscatter(dat4, x = "Growth_biomassdiff", y = "Maximum_caliper", 
                        add = "reg.line", conf.int = TRUE, 
                        xlab = "RGR (g/month) ", ylab = "Length (mm)")

# 2 Differences in propagules Area

correplot4 <- ggscatter(dat4, x = "Prop_established", y = "Area_diff_prop", 
                        add = "reg.line", conf.int = TRUE, 
                        xlab = "Number of established shoots ", ylab = "Wet *vs.* dry area")+
                        theme(axis.title.y = element_markdown())

correplot5 <- ggscatter(dat4, x = "Area_exp", y = "Area_diff_prop", 
                        add = "reg.line", conf.int = TRUE, 
                        xlab = "Colonized surface", ylab = "Wet *vs.* dry area")+
                        theme(axis.title.y = element_markdown())

correplot6 <- ggscatter(dat4, x = "Growth_biomassdiff", y = "Area_diff_prop", 
                        add = "reg.line", conf.int = TRUE, 
                        xlab = "RGR (g/month) ", ylab = "Wet *vs.* dry area")+
                        theme(axis.title.y = ggtext::element_markdown())

# 3 Differences in propagule length 

correplot7 <- ggscatter(dat4, x = "Prop_established", y = "FeretD_diff", 
                        add = "reg.line", conf.int = TRUE, 
                        xlab = "Number of established shoots ", ylab = "Wet *vs.* dry length")+
                        theme(axis.title.y = ggtext::element_markdown())
        

correplot8 <- ggscatter(dat4, x = "Area_exp", y = "FeretD_diff", 
                        add = "reg.line", conf.int = TRUE, 
                        xlab = "Colonized surface ", ylab = "Wet *vs.* dry length")+
                        theme(axis.title.y = ggtext::element_markdown())

correplot9 <- ggscatter(dat4, x = "Growth_biomassdiff", y = "FeretD_diff", 
                        add = "reg.line", conf.int = TRUE, 
                        xlab = "RGR (g/month) ", ylab = "Wet *vs.* dry length")+
                        theme(axis.title.y = ggtext::element_markdown())
  
# We get the R of pearson correlations and the p adjusted value to plot them properly in the mosaic plot
mycor1 <- corr.test(dat4$Perc_Bviable, dat4$Shoots_mean, method = "pearson", adjust = "holm", alpha = .05, ci= TRUE)
mycor2 <- corr.test(dat4$Growth_biomassdiff, dat4$Index_viability, method = "pearson", adjust = "holm", alpha = .05, ci= TRUE)
mycor3 <- corr.test(dat4$Growth_biomassdiff, dat4$Maximum_caliper, method = "pearson", adjust = "holm", alpha = .05, ci= TRUE)
mycor4 <- corr.test(dat4$Prop_established, dat4$Area_diff_prop, method = "pearson", adjust = "holm", alpha = .05, ci= TRUE)
mycor5 <- corr.test(dat4$Area_exp, dat4$Area_diff_prop, method = "pearson", adjust = "holm", alpha = .05, ci= TRUE)
mycor6 <- corr.test(dat4$Growth_biomassdiff, dat4$Area_diff_prop, method = "pearson", adjust = "holm", alpha = .05, ci= TRUE)
mycor7 <- corr.test(dat4$Prop_established, dat4$FeretD_diff, method = "pearson", adjust = "holm", alpha = .05, ci= TRUE)
mycor8 <- corr.test(dat4$Area_exp, dat4$FeretD_diff, method = "pearson", adjust = "holm", alpha = .05, ci= TRUE)
mycor9 <- corr.test(dat4$Growth_biomassdiff, dat4$FeretD_diff, method = "pearson", adjust = "holm", alpha = .05, ci= TRUE)

# Correction of margins for a proper visualization of the complex mosaic plot
correplot1 <- correplot1 + theme(plot.margin = unit(c(1, 0.25, 0.25, 0.25), "cm"))
correplot2 <- correplot2 + theme(plot.margin = unit(c(1, 0.25, 0.25, 0.25), "cm"))
correplot3 <- correplot3 + theme(plot.margin = unit(c(1, 0.25, 0.25, 0.25), "cm"))
correplot4 <- correplot4 + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))
correplot5 <- correplot5 + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))
correplot6 <- correplot6 + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))
correplot7 <- correplot7 + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))
correplot8 <- correplot8 + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))
correplot9 <- correplot9 + theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "cm"))

### 4.4 - Mosaic plot: Export -----------------------------------------------------

png(file="myfigS1.png",width=893,height=648, units = "px", res = 100) # it will create the png file on the current directory. It helps to improve resolution and reproducibility when saving the file.

plot_grid(correplot1, correplot2,correplot3,correplot4, correplot5,correplot6,correplot7, correplot8,correplot9, ncol=3, nrow = 3, align = "hv")

grid.text(label = paste(paste("R=",round(mycor1$r,3), sep = " "),paste("p=",round(mycor1$p.adj,3),sep = " "),sep = ", "), x= 0.26 , y=0.765, rot = 0, gp= gpar(fontsize=10, fontface = "bold"))
grid.text(label = paste(paste("R=",round(mycor2$r,3), sep = " "),paste("p=",round(mycor2$p.adj,3),sep = " "),sep = ", "), x= 0.592 , y=0.765, rot = 0, gp= gpar(fontsize=10, fontface = "bold"))
grid.text(label = paste(paste("R=",round(mycor3$r,3), sep = " "),paste("p=",round(mycor3$p.adj,3),sep = " "),sep = ", "), x= 0.925 , y=0.765, rot = 0, gp= gpar(fontsize=10, fontface = "bold"))
grid.text(label = "A", x= 0.025 , y=0.975, rot = 0, gp= gpar(fontsize=14, fontface = "bold"))


grid.text(label = paste(paste("R=",round(mycor4$r,3), sep = " "),paste("p=",round(mycor4$p.adj,3),sep = " "),sep = ", "), x= 0.26 , y=0.43, rot = 0, gp= gpar(fontsize=10, fontface = "bold"))
grid.text(label = paste(paste("R=",round(mycor5$r,3), sep = " "),paste("p=",round(mycor5$p.adj,3),sep = " "),sep = ", "), x= 0.592 , y=0.43, rot = 0, gp= gpar(fontsize=10, fontface = "bold"))
grid.text(label = paste(paste("R=",round(mycor6$r,3), sep = " "),paste("p=",round(mycor6$p.adj,3),sep = " "),sep = ", "), x= 0.925 , y=0.43, rot = 0, gp= gpar(fontsize=10, fontface = "bold"))
grid.text(label = "B", x= 0.025 , y=0.665, rot = 0, gp= gpar(fontsize=14, fontface = "bold"))


grid.text(label = paste(paste("R=",round(mycor7$r,3), sep = " "),paste("p=",round(mycor7$p.adj,3),sep = " "),sep = ", "), x= 0.26 , y=0.095, rot = 0, gp= gpar(fontsize=10, fontface = "bold"))
grid.text(label = paste(paste("R=",round(mycor8$r,3), sep = " "),paste("p=",round(mycor8$p.adj,3),sep = " "),sep = ", "), x= 0.592 , y=0.095, rot = 0, gp= gpar(fontsize=10, fontface = "bold"))
grid.text(label = paste(paste("R= 0.520"),paste("p=",round(mycor9$p.adj,3),sep = " "),sep = ", "), x= 0.925 , y=0.095, rot = 0, gp= gpar(fontsize=10, fontface = "bold")) # we got the r value before, so we substitute by 0.520 instead as 0.52 (resulted from rounded). If not done, it will be shorter than the rest (2 decimals instead of 3)
grid.text(label = "C", x= 0.025 , y=0.315, rot = 0, gp= gpar(fontsize=14, fontface = "bold"))

dev.off() # once device is set off, the plot with the measurements specified is exported correctly to the current directory

# Cleaning the code ####

rm(list = ls()) # to clean data

