library(simsurv)
library(tableone)
library(autoCovariateSelection)
library(survival)
library(glmnet)
library(parallel)
library(doParallel)
library(ranger)

setwd("~/GitHub/hdPS_survival/")

# Analytic data
rm(list = ls())
load("Data/simdata.RData")

#################### Random proxies from different website ######################
# Proxies from different website such as from cms.gov, gov.bc.ca
proxies <- read.csv("Data/Proxy.csv", header = T)
head(proxies)

## 3-digit diagnostic codes from hospital database
dat.diag <- subset(proxies, dim == "diag")
dat.diag$code <- substr(dat.diag$code, start = 1, stop = 3)
dat.diag$code[dat.diag$code==""] <- NA
dat.diag <- na.omit(dat.diag)

## 3-digit procedure codes from hospital database
dat.proc <- subset(proxies, dim == "proc")
dat.proc$code <- substr(dat.proc$code, start = 1, stop = 3)
dat.proc$code[dat.proc$code==""] <- NA
dat.proc <- na.omit(dat.proc)

## 3-digit icd codes from physician claim database
dat.msp <- subset(proxies, dim == "msp")
dat.msp$code <- substr(dat.msp$code, start = 1, stop = 3)
dat.msp$code[dat.msp$code==""] <- NA
dat.msp <- na.omit(dat.msp)

## DINPIN from drug dispensation database
dat.din <- subset(proxies, dim == "din")
dat.din$code[dat.din$code==""] <- NA
dat.din <- na.omit(dat.din)

# Merge all dimensions
dat.proxy <- rbind(dat.diag, dat.proc, dat.msp, dat.din)

# Drop missing proxies 
dat.proxy <- na.omit(dat.proxy)
table(dat.proxy$dim, useNA = "always") # 4 data dimensions

################### Converting proxies into empirically identified variables ################
id <- simdat$studyid

#### Generate candidate empirical covariates ####
step1 <- get_candidate_covariates(df = dat.proxy, domainVarname = "dim", 
                                  eventCodeVarname = "code", 
                                  patientIdVarname = "studyid", 
                                  patientIdVector = id, 
                                  n = 200, 
                                  min_num_patients = 20)
out1 <- step1$covars_data

#### Assessing recurrence of codes ####
all.equal(id, step1$patientIds)

step2 <- get_recurrence_covariates(df = out1, 
                                   eventCodeVarname = "code", 
                                   patientIdVarname = "studyid",
                                   patientIdVector = id)
out2 <- step2$recurrence_data

# Empirical covariates
vars.empirical <- names(out2)[-1]

#### Merge empirical covariates with the analytic dataset ####
dat.all <- merge(simdat, out2, by = "studyid", all.x = T)

################### Prioritizing empirical covariates using Bross formula ################
# Apply Bross formula
step3 <- get_prioritised_covariates(df = out2, patientIdVarname = "studyid", 
                                    exposureVector = ifelse(simdat$arthritis == "Yes", 1, 0),
                                    outcomeVector = simdat$cvd,
                                    patientIdVector = simdat$studyid,
                                    k = ncol(out2) - 1)
empvars.bross <- data.frame(vars = names(step3$multiplicative_bias), m.bias = step3$multiplicative_bias)
empvars.bross <- empvars.bross[order(empvars.bross$m.bias, decreasing = T),]
rownames(empvars.bross) <- NULL
head(empvars.bross)

## Top 200 empirical covariates section based on Bross formula
top200.bross <- empvars.bross$vars[1:200]

################### Prioritizing empirical covariates using Cox-LASSO ################
# Formula with only empirical covariates
formula.out <- as.formula(paste("Surv(follow_up, cvd) ~ ", paste(vars.empirical, collapse = " + ")))

# Model matrix for fitting Cox with LASSO regularization
X <- model.matrix(formula.out, data = dat.all)[,-1]
Y <- as.matrix(data.frame(time = dat.all$follow_up, status = dat.all$cvd))

# Detect the number of cores
n_cores <- parallel::detectCores()

# Create a cluster of cores
cl <- makeCluster(n_cores - 1)

# Register the cluster for parallel processing
registerDoParallel(cl)

# Hyperparameter tuning using 5-fold cross-validation 
set.seed(123)
fit.lasso <- cv.glmnet(x = X, y = Y, nfolds = 5, parallel = T, alpha = 1, family = "cox")
stopCluster(cl)

plot(fit.lasso)

## Best lambda
fit.lasso$lambda.min

# Variable ranking based on Cox-LASSO
#empvars.lasso <- coef(fit.lasso, s = fit.lasso$lambda.min) 
#Since proxies were random and unrelated to the simulated data, LASSO produced all zero coefficients.
#Let choose an arbitrary value as to demonstrate the process of variable selection 

empvars.lasso <- coef(fit.lasso, s = exp(-6)) 
empvars.lasso <- data.frame(as.matrix(empvars.lasso))
empvars.lasso <- data.frame(vars = rownames(empvars.lasso), coef = empvars.lasso)
colnames(empvars.lasso) <- c("vars", "coef")
rownames(empvars.lasso) <- NULL

## Rank empirical covariates based on absolute value of log hazard ratio
empvars.lasso$coef.abs <- abs(empvars.lasso$coef)
empvars.lasso <- empvars.lasso[order(empvars.lasso$coef.abs, decreasing = T),]
head(empvars.lasso)

## Top 200 empirical covariates section based on Cox-LASSO
top200.lasso <- empvars.lasso$vars[1:200]

################### Prioritizing empirical covariates using random survival forest ################
# Formula with only empirical covariates
formula.out <- as.formula(paste("Surv(follow_up, cvd) ~ ", paste(vars.empirical, collapse = " + ")))

# Hyperparameter tuning using grid search approach
# Grid with 1000 models - huge time consuming
# grid.search <- expand.grid(mtry = 1:10, node.size = 1:10,
#                           num.trees = seq(50,500,50),
#                           oob_error = 0)

# Grid with 36 models as an exercise
grid.search <- expand.grid(mtry = 5:7, node.size = 1:3, 
                           num.trees = seq(200,500,100), 
                           oob_error = 0)
head(grid.search)

## Calculate prediction error for each grid 
for(ii in 1:nrow(grid.search)) {
  # Model on training set with grid
  fit.rsf.tune <- ranger(formula = formula.out, 
                        data = dat.all, 
                        num.trees = grid.search$num.trees[ii],
                        mtry = grid.search$mtry[ii], 
                        min.node.size = grid.search$node.size[ii],
                        importance = 'impurity')
  
  # Add Out-of-bag (OOB) error to each grid
  grid.search$oob_error[ii] <- sqrt(fit.rsf.tune$prediction.error)
}
head(grid.search)

# Best hyperparameters
position <- which.min(grid.search$oob_error)
grid.search[position,]

# Model after tuning hyperparameters
set.seed(123)
fit.rsf <- ranger(formula.out, data = dat.all, 
                   num.trees = grid.search$num.trees[position], 
                   min.node.size = grid.search$node.size[position], 
                   mtry = grid.search$mtry[position], 
                   importance = "impurity")

## Rank empirical covariates based on variable importance
empvars.rsf <- data.frame(vars = names(fit.rsf$variable.importance), vim = fit.rsf$variable.importance)
rownames(empvars.rsf) <- NULL
empvars.rsf <- empvars.rsf[order(empvars.rsf$vim, decreasing = T),]
head(empvars.rsf)

## Top 200 empirical covariates section based on random survival forest
top200.rsf <- empvars.rsf$vars[1:200]

################### Save ################
save(
  # Simulated data with only investigator-specified covariates
  simdat, 
  
  # Simulated data with investigator-specified covariates and empirical covariates
  dat.all, 
  
  # Empirical covariates
  emp.vars, 
  
  # Empirical covariates ranked by Bross formula, Cox-LASSO, and RSF
  empvars.bross, empvars.lasso, empvars.rsf, 
  
  # Top 200 empirical covariates ranked by Bross formula, Cox-LASSO, and RSF
  top200.bross, top200.lasso, top200.rsf, 
  
  # Save as RData
  file = "Data/simdata_with_empirical.RData")
