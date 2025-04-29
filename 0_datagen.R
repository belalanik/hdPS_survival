library(simsurv)
library(tableone)

setwd("~/GitHub/hdPS_survival/")

#################### Generating a survival outcome with a binary exposure ######################
## studyid: Unique identifier

## Survival outcome
# cvd: CVD status
# follow_up: Time-to-cvd in years

## Binary exposure
# arthritis: Arthritis status 

## Investigator-specified confounders
# age: Age in years
# sex: Whether male or female
# comorbidity: Whether any Comorbidities

rm(list = ls())

# Sample size
n <- 3000

# Seed
set.seed(229)

# Confounders
age <- round(rnorm(n, mean = 50, sd = 10))
sex <- rbinom(n, 1, prob = 0.4)
comorbidity <- rbinom(n, 1, prob = 0.25)

# Exposure
arthritis <- rbinom(n, 1, prob = plogis(-3.35 + 0.05*age - 0.5*sex + 0.4*comorbidity))

# Exposure and confounders in a data frame
vars <- data.frame(arthritis, age, sex, comorbidity)

# Generating data
simdat <- simsurv(dist = "weibull", lambda = 0.1, gammas = 0.1, maxt = 20, x = vars,
                betas = c(arthritis = 0.7, age = 0.03, sex = -1.9, comorbidity = 0.15))
simdat <- cbind(simdat, vars)

# Rename column names
names(simdat)[names(simdat) == "id"] <- "studyid"
names(simdat)[names(simdat) == "eventtime"] <- "follow_up"
names(simdat)[names(simdat) == "status"] <- "cvd"

# Follow-up time
simdat$follow_up <- round(simdat$follow_up, 2)
simdat$follow_up[simdat$follow_up == 0] <- 0.01
head(simdat)

# Recoding
simdat$arthritis <- factor(simdat$arthritis, levels = c(0, 1), labels = c("No", "Yes"))
simdat$sex <- factor(simdat$sex, levels = c(0, 1), labels = c("Female", "Male"))
simdat$comorbidity <- factor(simdat$comorbidity, levels = c(0, 1), labels = c("No", "Yes"))

# Table 1
tab1 <- CreateTableOne(c("cvd", "follow_up", "age", "sex", "comorbidity"), strata = "arthritis", 
                       data = simdat, test = FALSE, factorVars = "cvd")
print(tab1, smd = T)

#################### Save dataset ######################
save(simdat, file = "Data/simdata.RData")
