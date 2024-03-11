## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "Input-"
)

## ----warning=FALSE, message=FALSE---------------------------------------------
library(cophescan)

## ----warning=FALSE, message=FALSE---------------------------------------------
data("cophe_multi_trait_data")
trait_dat = cophe_multi_trait_data$summ_stat$Trait_1
str(trait_dat)

## -----------------------------------------------------------------------------
trait_dat$LD = cophe_multi_trait_data$LD
str(trait_dat$LD[1:10, 1:10])

