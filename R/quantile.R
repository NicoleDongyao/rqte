path="C:/Users/nicole/Desktop/code/rqte/R"
setwd(path)
source('causalForest.R')
source('causalTree.anova.R')
source('causalTree.anova2.R')
source('causalTree.branch.R')
source('causalTree.class.R')
source('causalTree.control.R')
source('causalTree.exp.R')
source('causalTree.matrix.R')
source('causalTree.poisson.R')
source('causalTree.R')
source('causalTreecallback.R')
source('causalTreeco.R')
source('create.data.frame.R')
source('create.sample.R')
source('est.causalTree.R')
source('estimate.causalTree.R')
source('estimate.R')
source('formatg.R')
source('honestModel.R')
source('importance.R')
source('init.causalForest.R')
source('labels.causalTree.R')
source('meanvar.causalTree.R')
source('model.frame.causalTree.R')
source('modify.estimations.R')
source('na.causalTree.R')
source('path.causalTree.R')
source('plot.causalTree.R')
source('plotcp.R')
source('pred.causalTree.R')
source('predict.causalForest.R')
source('predict.causalTree.R')
source('print.causalTree.R')
source('printcp.R')
source('propensityForest.R')
source('prune.causalTree.R')
source('prune.R')
source('reestimate.tau.R')
source('refit.causalTree.R')
source('residuals.causalTree.R')
source('snip.causalTree.mouse.R')
source('snip.causalTree.R')
source('summary.causalTree.R')
source('text.causalTree.R')
source('xpred.causalTree.R')
source('zzz.R')
################################################################################
library(usethis)
library(devtools)  #Load the installed "devtools" package
library(rpart)   #Load the installed "rparts" package
library(rpart.plot)
library(rattle)
library(ggplot2)

################################################################################

# qqplot
rm(list = ls())

sigma = 1
d = 20
n = 800

tree_mult = 4
n.test = 1000
simu.reps = 2

# m(x)
baseline = function(x) {
  2 * (x[1] - 0.5)
}

propensity = function(x) {
  0.25 + dbeta(x[1], 2, 4)/4
}

X.test = matrix(runif(n.test * d, 0, 1), n.test, d)

single.run = function(n) {
  X = matrix(runif(n * d, 0, 1), n, d) # features
  e = apply(X, 1, propensity)
  W = rbinom(n, 1, e) #treatment condition
  
  # no treatment effect
  Y = apply(X, 1, baseline) + sigma * rnorm(n)
  
  #
  # random forest
  #
  
  forest = propensityForest(X, Y, W, num.trees = round(tree_mult * n), sample.size = n^(0.8), nodesize = 1)
  predictions = predict(forest, X.test)
  predictions
}


results = lapply(1:simu.reps, function(rr)single.run(n))

save.image("figure2b_gauss.RData")

preds = Reduce(cbind, results)
preds.std = t(apply(preds, 1, function(xx) (xx - mean(xx)) / sd(xx)))
ymax = max(abs(preds.std))

gauss.quantiles = qnorm((1:simu.reps) / (simu.reps + 1))

#plot(NA, NA, xlim = range(gauss.quantiles), ylim = c(-ymax, ymax), xlab = "theoretical quantiles", ylab = "sample quantiles")
#for(iter in 1:nrow(preds.std)) {
#	points(gauss.quantiles, sort(preds.std[iter,]))
#}
#abline(0, 1, lwd = 3, col = 2)

all.data = data.frame(Reduce(rbind, lapply(1:nrow(preds.std),
                                           function(iter) cbind(Q= gauss.quantiles, T=sort(preds.std[iter,]))
)))


pdf("output/gauss_qqplot.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(NA, NA, xlim = range(c(gauss.quantiles, -2, 2)), ylim = c(-ymax, ymax), xlab = "theoretical quantiles", ylab = "sample quantiles")
abline(0, 1, lwd = 2, lty = 3)
boxplot(T~Q, all.data, xaxt="n", at = gauss.quantiles, boxwex=0.2, add=TRUE)
xaxes = seq(-2, 2, by = 1)
axis(1, at=xaxes, label=xaxes)
par=pardef
dev.off()






















