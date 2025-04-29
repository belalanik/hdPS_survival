library(survival)
library(Publish)
library(jtools)
library(survey)

setwd("~/GitHub/hdPS_survival/")

################################# Data ################################
rm(list = ls())
#load("Data/simdata.RData")
load("Data/simdata_with_empirical.RData")

################################# Traditional regression adjustment ################################
vars.investigator <- c("age", "sex", "comorbidity")

# Formula
Formula.out <- as.formula(paste0("Surv(follow_up, cvd) ~ arthritis + ", paste(vars.investigator, collapse = "+")))

# Adjusted for investigator-specified confounders
fit.traditional <- coxph(Formula.out, data = dat.all)
publish(fit.traditional, print = F)$regressionTable[1:2,]


################################# hdPS ################################
# PS model specification with investigator-specified and top 200 empirical covariates using Cox-LASSO
vars <- c(vars.investigator, top200.lasso)
Formula.exp <- as.formula(paste0("I(arthritis == 'Yes') ~ ", paste(vars, collapse = "+")))

# Fitting the PS model - logistic
fit.ps <- glm(Formula.exp, data = dat.all, family = binomial)
fit.ps$coefficients[is.na(fit.ps$coefficients)] <- 0
dat.all$ps <- predict(fit.ps, type = "response")
summary(dat.all$ps)

# Stabilized weight calculation
dat.all$sipw <- ifelse(dat.all$arthritis == "Yes", mean(dat.all$arthritis == "Yes")/dat.all$ps, 
                      ( 1 - mean(dat.all$arthritis == "Yes")) / (1 - dat.all$ps))
summary(dat.all$sipw)

# Create a design
svy.design <- svydesign(id = ~1, weights = ~sipw, data = dat.all)

# Balance checking
tab1b <- svyCreateTableOne(vars = vars.investigator, strata = "arthritis", data = svy.design, test = FALSE)
print(tab1b, smd = TRUE)

# Outcome analysis - Cox PH
fit.hdps <- coxph(Surv(follow_up, cvd) ~ arthritis, data = dat.all, weights = sipw)
publish(fit.hdps, pvalue.method = "robust", confint.method = "robust", print = F)$regressionTable[1:2,]

# fit.hdps1 <- svycoxph(Surv(follow_up, cvd) ~ arthritis, design = svy.design)
# publish(fit.hdps1, print = F)$regressionTable[1:2,]

