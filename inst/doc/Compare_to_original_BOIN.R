## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bfboin)

## -----------------------------------------------------------------------------
get.oc.bf(ntrial = 10000,
          seed = 9,
          target = 0.25,
          p.true = c(0.1, 0.5),
          ncohort = 10,
          cohortsize = 3,
          n.earlystop = 9,
          startdose = 1,
          titration = FALSE,
          cutoff.eli = 0.95,
          extrasafe = TRUE,
          offset = 0.1,
          boundMTD=FALSE,
          n.cap = 12,
          end.backfill = TRUE,
          n.per.month = 1,
          dlt.window = 1,
          p.response.true = c(0.00001, 0.00001))


BOIN::get.oc(target = 0.25,
             p.true = c(0.1, 0.5),
             ncohort = 10,
             cohortsize = 3,
             n.earlystop = 9,
             startdose = 1,
             titration = FALSE,
             cutoff.eli = 0.95,
             extrasafe = TRUE,
             offset = 0.1,
             boundMTD=FALSE,
             ntrial = 10000, seed = 6)

