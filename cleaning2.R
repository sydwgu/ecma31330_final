## -----------------------------------------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(haven)
library(glmnet)
library(knitr)
library(kableExtra)
library(caret)
library(stats)
library(fastDummies)


## -----------------------------------------------------------------------------------------------------------------
cex <- read_dta("~/Desktop/metrics_and_ml/ecma31330_final/113270/cexall.dta")
psid <- read_dta("~/Desktop/metrics_and_ml/ecma31330_final/113270/data.dta")
natpr <- read_dta("~/Desktop/metrics_and_ml/ecma31330_final/113270/natpr.dta")
tax9192 <- read_dta("~/Desktop/metrics_and_ml/ecma31330_final/113270/tax9192.dta")
finprc <- read_dta("~/Desktop/metrics_and_ml/ecma31330_final/113270/finprice.dta")


## -----------------------------------------------------------------------------------------------------------------
psid <- merge(psid, tax9192, by = c("id", "year"), all.x=TRUE)
psid$ftax[psid$year %in% c(1991, 1992)] <- psid$fiitax[psid$year %in% c(1991, 1992)]
psid <- psid[, !(names(psid) %in% c("fiitax"))]


## -----------------------------------------------------------------------------------------------------------------
psid <- merge(psid, natpr, by = "year")


## -----------------------------------------------------------------------------------------------------------------
cex <- cex %>%
  mutate(hours = weekh * hourh)
cex <- cex %>%
  mutate(hourw = weekw * hourw)


## -----------------------------------------------------------------------------------------------------------------
cex <- cex %>%
  mutate(
    ftax = income - income_at,  
    asset = pf + pfaw          
  )


## -----------------------------------------------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------------------------------------------
psid$yb <- psid$year - psid$age


## -----------------------------------------------------------------------------------------------------------------
psid <- psid %>%
  mutate(income = ly + wly + y) %>%
  mutate(food = food + fout)


## -----------------------------------------------------------------------------------------------------------------
cex <- cex %>% select(-c(income))


## -----------------------------------------------------------------------------------------------------------------
cex <- cex %>%
  rename(fsize = ncomp,
         kids = kid,
         educ = educh,
         income = income_at,
         )


## -----------------------------------------------------------------------------------------------------------------
cex <- cex %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))


## -----------------------------------------------------------------------------------------------------------------
cex <- cex %>%
  mutate(across(c(hours, hourw), ~ ifelse(is.na(.), 0, .)))


## -----------------------------------------------------------------------------------------------------------------
psid$kids <- replace(psid$kids, is.na(psid$kids), 0)


## -----------------------------------------------------------------------------------------------------------------
cex <- cex[cex$income >= 0, ]
psid <- psid[psid$income >= 0, ]


## -----------------------------------------------------------------------------------------------------------------
cex$hs_drop <- ifelse(cex$educ == 1, 1, 0)
cex$hs_grad <- ifelse(cex$educ > 1 & cex$educ < 7, 1, 0)

cex$"1981" <- ifelse(cex$year == 1981, 1, 0)
cex$"1982" <- ifelse(cex$year == 1982, 1, 0)
cex$"1983" <- ifelse(cex$year == 1983, 1, 0)
cex$"1984" <- ifelse(cex$year == 1984, 1, 0)
cex$"1985" <- ifelse(cex$year == 1985, 1, 0)
cex$"1986" <- ifelse(cex$year == 1986, 1, 0)
cex$"1987" <- ifelse(cex$year == 1987, 1, 0)
cex$"1988" <- ifelse(cex$year == 1988, 1, 0)
cex$"1989" <- ifelse(cex$year == 1989, 1, 0)
cex$"1990" <- ifelse(cex$year == 1990, 1, 0)
cex$"1991" <- ifelse(cex$year == 1991, 1, 0)
cex$"1992" <- ifelse(cex$year == 1992, 1, 0)

cex$onec <- ifelse(cex$kids == 1, 1, 0)
cex$twoc <- ifelse(cex$kids == 2, 1, 0)
cex$three_plusc <- ifelse(cex$kids >= 3, 1, 0)

cex$northeast <- ifelse(cex$region == 1, 1, 0)
cex$midwest <- ifelse(cex$region == 2, 1, 0)
cex$south <- ifelse(cex$region == 3, 1, 0)

cex$"55-59" <- ifelse(cex$yb >= 1955 & cex$yb <= 1959, 1, 0)
cex$"50-54" <- ifelse(cex$yb >= 1950 & cex$yb <= 1954, 1, 0)
cex$"45-49" <- ifelse(cex$yb >= 1945 & cex$yb <= 1949, 1, 0)
cex$"40-44" <- ifelse(cex$yb >= 1940 & cex$yb <= 1944, 1, 0)
cex$"35-39" <- ifelse(cex$yb >= 1935 & cex$yb <= 1939, 1, 0)
cex$"30-34" <- ifelse(cex$yb >= 1930 & cex$yb <= 1934, 1, 0)
cex$"25-29" <- ifelse(cex$yb >= 1925 & cex$yb <= 1929, 1, 0)

cex$white <- ifelse(cex$race == 1, 1, 0)


## -----------------------------------------------------------------------------------------------------------------
psid$hs_drop <- ifelse(psid$educ == 1, 1, 0)
psid$hs_grad <- ifelse(psid$educ > 1 & psid$educ < 7, 1, 0)

psid$"1981" <- ifelse(psid$year == 1981, 1, 0)
psid$"1982" <- ifelse(psid$year == 1982, 1, 0)
psid$"1983" <- ifelse(psid$year == 1983, 1, 0)
psid$"1984" <- ifelse(psid$year == 1984, 1, 0)
psid$"1985" <- ifelse(psid$year == 1985, 1, 0)
psid$"1986" <- ifelse(psid$year == 1986, 1, 0)
psid$"1987" <- ifelse(psid$year == 1987, 1, 0)
psid$"1988" <- ifelse(psid$year == 1988, 1, 0)
psid$"1989" <- ifelse(psid$year == 1989, 1, 0)
psid$"1990" <- ifelse(psid$year == 1990, 1, 0)
psid$"1991" <- ifelse(psid$year == 1991, 1, 0)
psid$"1992" <- ifelse(psid$year == 1992, 1, 0)

psid$onec <- ifelse(psid$kids == 1, 1, 0)
psid$twoc <- ifelse(psid$kids == 2, 1, 0)
psid$three_plusc <- ifelse(psid$kids >= 3, 1, 0)

psid$northeast <- ifelse(psid$region == 1, 1, 0)
psid$midwest <- ifelse(psid$region == 2, 1, 0)
psid$south <- ifelse(psid$region == 3, 1, 0)

psid$"55-59" <- ifelse(psid$yb >= 1955 & psid$yb <= 1959, 1, 0)
psid$"50-54" <- ifelse(psid$yb >= 1950 & psid$yb <= 1954, 1, 0)
psid$"45-49" <- ifelse(psid$yb >= 1945 & psid$yb <= 1949, 1, 0)
psid$"40-44" <- ifelse(psid$yb >= 1940 & psid$yb <= 1944, 1, 0)
psid$"35-39" <- ifelse(psid$yb >= 1935 & psid$yb <= 1939, 1, 0)
psid$"30-34" <- ifelse(psid$yb >= 1930 & psid$yb <= 1934, 1, 0)
psid$"25-29" <- ifelse(psid$yb >= 1925 & psid$yb <= 1929, 1, 0)

psid$white <- ifelse(psid$race == 1, 1, 0)



## -----------------------------------------------------------------------------------------------------------------
cex <- cex %>% select(-c(complete, cons1, cons2, educw, family, famtype, 
                         income_net, wageh, wagew, weekh, weekw, ndurplus, 
                         hourh, pfoodin, pfin, pfaw, pf))

