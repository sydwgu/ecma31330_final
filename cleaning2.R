## -------------------------------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(haven)
library(glmnet)
library(knitr)
library(kableExtra)
library(caret)
library(stats)


## -------------------------------------------------------------------------------------------------------
cex <- read_dta("~/Desktop/metrics_and_ml/ecma31330_final/113270/cexall.dta")
psid <- read_dta("~/Desktop/metrics_and_ml/ecma31330_final/113270/data.dta")


## -------------------------------------------------------------------------------------------------------
psid <- dummy_cols(psid, select_columns = c("race", "state", "smsa"))

# Remove columns that represent >2 categorical variables
psid <- psid %>%
  select(-c(vet, disab, outkid, unioni, unionj, race, state, smsa))


## -------------------------------------------------------------------------------------------------------
cex <- dummy_cols(cex, select_columns = c("race", "region", "famtype"))

# Remove columns that represent >2 categorical variables
cex <- cex %>%
  select(-c(race, region, famtype))
head(cex)


## -------------------------------------------------------------------------------------------------------
natpr <- read_dta("~/Desktop/metrics_and_ml/ecma31330_final/113270/natpr.dta")
tax9192 <- read_dta("~/Desktop/metrics_and_ml/ecma31330_final/113270/tax9192.dta")
finprc <- read_dta("~/Desktop/metrics_and_ml/ecma31330_final/113270/finprice.dta")


## -------------------------------------------------------------------------------------------------------
psid <- merge(psid, tax9192, by = c("id", "year"))
psid$ftax[psid$year %in% c(1991, 1992)] <- psid$fiitax[psid$year %in% c(1991, 1992)]
psid <- psid[, !(names(psid) %in% c("fiitax"))]


## -------------------------------------------------------------------------------------------------------
psid <- merge(psid, natpr, by = "year")


## -------------------------------------------------------------------------------------------------------
cex <- cex %>%
  mutate(hours = weekh * hourh)
cex <- cex %>%
  mutate(hourw = weekw * hourw)


## -------------------------------------------------------------------------------------------------------
# for PSID 
psid <- psid %>%
  mutate(total_income = tyhw + tyoth)


## -------------------------------------------------------------------------------------------------------
psid <- psid %>%
  mutate(hw_laborinc = ly + wly)


## -------------------------------------------------------------------------------------------------------
cex <- cex %>%
  mutate(hw_wages = wageh + wagew)


## -------------------------------------------------------------------------------------------------------
psid <- psid %>%
  mutate(hwage_inc = wages * hours)


## -------------------------------------------------------------------------------------------------------
cex <- cex %>%
  rename(fsize = ncomp)
cex <- cex %>%
  rename(kids = kid)


## -------------------------------------------------------------------------------------------------------
cex <- cex %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))


## -------------------------------------------------------------------------------------------------------
cex <- cex %>%
  mutate(across(c(hourh, hourw, wageh, wagew, hours, hw_wages), ~ ifelse(is.na(.), 0, .)))

