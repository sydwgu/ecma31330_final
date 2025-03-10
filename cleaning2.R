## -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(haven)
library(glmnet)
library(knitr)
library(kableExtra)
library(caret)
library(stats)
library(fastDummies)


## -----------------------------------------------------------------------------
cex <- read_dta("~/Desktop/metrics_and_ml/ecma31330_final/113270/cexall.dta")
psid <- read_dta("~/Desktop/metrics_and_ml/ecma31330_final/113270/data.dta")


## -----------------------------------------------------------------------------
natpr <- read_dta("~/Desktop/metrics_and_ml/ecma31330_final/113270/natpr.dta")
tax9192 <- read_dta("~/Desktop/metrics_and_ml/ecma31330_final/113270/tax9192.dta")
finprc <- read_dta("~/Desktop/metrics_and_ml/ecma31330_final/113270/finprice.dta")


## -----------------------------------------------------------------------------
psid <- merge(psid, tax9192, by = c("id", "year"))
psid$ftax[psid$year %in% c(1991, 1992)] <- psid$fiitax[psid$year %in% c(1991, 1992)]
psid <- psid[, !(names(psid) %in% c("fiitax"))]


## -----------------------------------------------------------------------------
psid <- merge(psid, natpr, by = "year")


## -----------------------------------------------------------------------------
cex <- cex %>%
  mutate(hours = weekh * hourh)
cex <- cex %>%
  mutate(hourw = weekw * hourw)


## -----------------------------------------------------------------------------
cex<- cex %>%
  mutate(
    ftax = income - income_at,  # Federal tax estimate
    asset = pf + pfaw           # Total asset estimate
  )


## -----------------------------------------------------------------------------
state_to_region <- function(state_code) {
  if (state_code %in% c(9, 23, 25, 33, 34, 36, 42, 44, 50)) {
    return(1)  # Northeast
  } else if (state_code %in% c(17, 18, 19, 20, 26, 27, 29, 31, 38, 39, 46, 55)) {
    return(2)  # Midwest
  } else if (state_code %in% c(1, 5, 10, 11, 12, 13, 21, 22, 24, 28, 37, 40, 45, 47, 48, 51, 54)) {
    return(3)  # South
  } else if (state_code %in% c(2, 4, 6, 8, 15, 16, 30, 32, 35, 41, 49, 53, 56)) {
    return(4)  # West
  } else {
    return(NA)  # Unknown state code
  }
}

psid <- psid %>%
  mutate(region = sapply(state, state_to_region))

head(psid)


## -----------------------------------------------------------------------------
psid <- psid %>%
  mutate(race = case_when(
    race == 1 ~ 1,  # White
    race == 2 ~ 2,  # Black or African American
    race == 3 ~ 3,  # Native American
    race == 4 ~ 4,  # Asian
    race == 5 ~ NA, # Latino origin
    race == 6 ~ NA, # Other color
    race == 7 ~ NA, # Other
    race == 9 ~ NA, # Don't know
    TRUE ~ NA      
  ))


## -----------------------------------------------------------------------------
psid <- psid %>%
  mutate(income = tyhw + tyoth)


## -----------------------------------------------------------------------------
cex <- cex %>% select(-c(complete, cons1, cons2, cpi, educw, family, famtype, income, income_net, pf, pfaw, pfin, pfood, pfout, wageh, wagew, weekh, weekw, yb))


## -----------------------------------------------------------------------------
cex <- cex %>%
  rename(fsize = ncomp,
         kids = kid,
         educ = educh,
         income = income_at,
         )


## -----------------------------------------------------------------------------
cex <- cex %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))


## -----------------------------------------------------------------------------
cex <- cex %>%
  mutate(across(c(hours, hourw), ~ ifelse(is.na(.), 0, .)))


## -----------------------------------------------------------------------------
#psid <- dummy_cols(psid, select_columns = c("race", "region"))
#psid <- psid %>%
 # select(-c(race, region, state))


## -----------------------------------------------------------------------------
#cex <- dummy_cols(cex, select_columns = c("race", "region"))
#cex <- cex %>%
 # select(-c(race, region))

