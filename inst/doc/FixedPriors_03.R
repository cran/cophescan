## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "FP-"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(cophescan)

## ----message=FALSE, warning=FALSE---------------------------------------------
data("cophe_multi_trait_data")
names(cophe_multi_trait_data)

## ----fig.width = 4, fig.height=4, fig.align = "center"------------------------
querytrait <- cophe_multi_trait_data$summ_stat[['Trait_1']] 
querysnpid <- cophe_multi_trait_data$querysnpid
LD <- cophe_multi_trait_data$LD

## ----regManhat, fig.width = 4, fig.height=4, fig.align = "center"-------------
# Additional  field named 'position' is required for the Manahattan plot. It is a numeric vector of chromosal positions
querytrait$position <- sapply(querytrait$snp, function(x) as.numeric(unlist(strsplit(x, "-"))[2]))
plot_trait_manhat(querytrait, querysnpid)

## ----fig.width = 4, fig.height=4, fig.align = "center"------------------------
# Run cophescan under a single causal variant assumption by providing the snpid of the known causal variant for trait 1 = querysnpid
res.single <- cophe.single(querytrait, querysnpid = querysnpid, querytrait='Trait_1')
summary(res.single)

## -----------------------------------------------------------------------------
res.single.predict <- cophe.hyp.predict(res.single)
(paste0('The predicted hypothesis is: ', res.single.predict$cophe.hyp.call, ' [PP.Hc =', round(res.single.predict$PP.Hc,3), ']' ))

## ----fig.width = 4, fig.height=4, fig.align = "center", message=FALSE---------
# Run cophescan with susie (multiple variants) by providing the snpid of the known causal variant for trait 1 = querysnpid
querytrait$LD <- LD
res.susie <- cophe.susie(querytrait, querysnpid = querysnpid, querytrait='Trait_1')
summary(res.susie)

res.susie.predict <- cophe.hyp.predict(res.susie)
(paste0('The predicted hypothesis is: ', res.susie.predict$cophe.hyp.call, ' [PP.Hc =', round(res.susie.predict$PP.Hc,3), ']' ))

