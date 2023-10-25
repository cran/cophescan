## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "HP-"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(cophescan)
library(dplyr)

## -----------------------------------------------------------------------------
data("cophe_multi_trait_data")
trait_dat = cophe_multi_trait_data$summ_stat$Trait_1
str(trait_dat)
querysnpid <- cophe_multi_trait_data$querysnpid
LD <- cophe_multi_trait_data$LD

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
## Hide print messages from coloc
res.multi.lbf <- list()
for (trait_idx in seq_along(cophe_multi_trait_data$summ_stat)){
  querytrait_ss <- cophe_multi_trait_data$summ_stat[[trait_idx]]
  # Here LD is the same
  querytrait_ss$LD <- LD
  trait_variant_pair <- paste0('Trait', trait_idx, '_', querysnpid)
  res.multi.lbf[[trait_variant_pair]] <- cophe.susie.lbf(querytrait_ss, querysnpid = querysnpid, querytrait = paste0('Trait_', trait_idx))
}

res.multi.lbf.df = bind_rows(res.multi.lbf)


## -----------------------------------------------------------------------------
head(res.multi.lbf.df)

## ----warning=F, message=F, fig.width = 8, fig.height=8, out.width = "75%",  fig.ncol = 3, fig.align = "center"----
# covar=FALSE
## Set covar to TRUE to include covariates
covar=TRUE
covar_vec = cophe_multi_trait_data$covar_vec
cophe.hier.res <- run_metrop_priors(res.multi.lbf.df, avg_posterior=TRUE, avg_pik = TRUE, covar_vec = covar_vec, covar = covar, nits = 30000)
names(cophe.hier.res)

## ----mcmcDiagnostics, warning=F, message=F, fig.width=5, fig.height=5, out.width="50%", fig.align="center"----

loglik <- cophe.hier.res$ll
parameters <- cophe.hier.res$parameters
col <- rgb(red = 0.4, green = 0.7, blue = 0.5, alpha = 0.8)
### store user parameters
old_par = par(no.readonly = TRUE)

## plot diagnostics
par(mfrow=c(2,2))
plot(seq_along(loglik), loglik, main="loglik",type="l", col=col, ylab = "ll", xlab="")
plot(seq_len(ncol(parameters)), parameters[1,], main="alpha",type="l", col=col, ylab = "alpha", xlab="")
plot(seq_len(ncol(parameters)), parameters[2,], main="beta",type="l", col=col, ylab = "beta", xlab="")
if (covar)
  plot(seq_len(ncol(parameters)), parameters[3,], main="gamma",type="l", col=col, ylab = "gamma", xlab="")

### reset user parameters
par(old_par)

## -----------------------------------------------------------------------------
res.post.prob = cbind(cophe.hier.res$avg.posterior, cophe.hier.res$data)

## -----------------------------------------------------------------------------
res.hier.predict <- cophe.hyp.predict(as.data.frame(res.post.prob ))
col_disp <- c( "PP.Hn", "PP.Ha", "PP.Hc", "nsnps", "querysnp", "querytrait", "typeBF",  "grp", "cophe.hyp.call")
knitr::kable(res.hier.predict[, col_disp], row.names = FALSE, digits=3)

## ----cophePlots, message=FALSE------------------------------------------------
res.plots = cophe_plot(res.hier.predict, traits.dat = cophe_multi_trait_data$summ_stat, querysnpid = querysnpid, query_trait_names = paste0('Trait_', 1:24))

# if (!require(ggpubr)) {
  # install.packages("ggpubr") 
# }
# ggpubr::ggarrange(res.plots$pval, res.plots$ppHa, res.plots$ppHc, ncol = 2, nrow = 2)

## ---- warning=FALSE, message=FALSE, echo=FALSE, out.width = "50%"-------------
res.plots$pval
res.plots$ppHa
res.plots$ppHc

