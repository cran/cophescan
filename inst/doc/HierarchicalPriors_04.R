## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "HP-"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(cophescan)
library(dplyr)
library(ggplot2)

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
covar=FALSE
covar_vec=rep(1, nrow(res.multi.lbf.df))
## Set covar to TRUE to include covariates, uncomment the following 2 lines 
# covar=TRUE
# covar_vec = cophe_multi_trait_data$covar_vec
cophe.hier.res <- run_metrop_priors(res.multi.lbf.df, avg_posterior=TRUE, avg_pik = TRUE, covar_vec = covar_vec, covar = covar, nits = 50000, thin = 5)
names(cophe.hier.res)

## ----message=FALSE------------------------------------------------------------
cophe.hier.res.chain.list <- lapply(1:4, function(x) 
  run_metrop_priors(res.multi.lbf.df, avg_posterior=TRUE, avg_pik = TRUE, 
                    covar_vec = covar_vec, covar = covar, nits = 50000, thin = 5))

## ----mcmcDiagnostics, warning=FALSE, message=FALSE, fig.width=5, fig.height=5, out.width="80%", fig.align="center"----
# Store user parameters
old_par <- par(no.readonly = TRUE)

# chain_colors <- c("#e63946c4", "#f1faee", "#a8dadc", "#457b9d" )
chain_colors <- c("#f4f1de", "#e07a5fb2", "#3d405bb2", "#81b29aa6")

layout(matrix(c(1, 2, 3, 4, 5, 5), ncol=2, byrow = TRUE), respect = TRUE, 
       heights = c(0.9, 0.9, 0.1))


matplot(sapply(cophe.hier.res.chain.list, function(x) x$ll), type = "l", 
        col = chain_colors, 
     main ="loglik", ylab = "ll", xlab = "Iteration", lty = 1)

y_ax <- c("alpha", "beta", "gamma")
num_pars <- ifelse(covar, 3, 2) 
for (idx in 1:num_pars) {
    matplot(sapply(cophe.hier.res.chain.list, function(x) x$parameters[idx, ]),
        type = "l", col = chain_colors,
        main = paste(y_ax[idx]), ylab = y_ax[idx], xlab = "Iteration", lty = 1
    )
}

if (!covar) {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
}

par(mar=c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("top", legend = paste("Chain", 1:4), col = chain_colors, lty = 1, lwd = 2, 
       horiz = TRUE, bty = "n")

# Reset user parameters
par(old_par)


## -----------------------------------------------------------------------------
res.post.prob = cbind(cophe.hier.res.chain.list[[1]]$avg.posterior, cophe.hier.res$data)

## -----------------------------------------------------------------------------
res.hier.predict <- cophe.hyp.predict(as.data.frame(res.post.prob ))
col_disp <- c( "PP.Hn", "PP.Ha", "PP.Hc", "nsnps", "querysnp", "querytrait",
               "typeBF",  "grp", "cophe.hyp.call")
knitr::kable(res.hier.predict[, col_disp], row.names = FALSE, digits=3)

## ----cophePlots, message=FALSE------------------------------------------------
res.plots = cophe_plot(res.hier.predict, traits.dat = cophe_multi_trait_data$summ_stat, querysnpid = querysnpid, query_trait_names = paste0('Trait_', 1:24))

# if (!require(ggpubr)) {
  # install.packages("ggpubr") 
# }
# ggpubr::ggarrange(res.plots$pval, res.plots$ppHa, res.plots$ppHc, ncol = 2, 
#                 nrow = 2)

## ----warning=FALSE, message=FALSE, echo=FALSE, out.width = "40%"--------------
res.plots$pval + theme(legend.position="bottom")
res.plots$ppHa + theme(legend.position="bottom")
res.plots$ppHc + theme(legend.position="bottom")

