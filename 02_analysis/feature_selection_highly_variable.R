#Feature selection

#Quantify per-gene variation
library(scran)
dec.lmmp_mf_sce <- modelGeneVar(lmmp_mf_sce,
                                block=lmmp_mf_sce$treatment)
#noweight comparison
dec.noweight <- modelGeneVar(lmmp_mf_sce, 
                             block=lmmp_mf_sce$treatment,
                             density.weights=FALSE)

plot(dec.lmmp_mf_sce$mean, dec.lmmp_mf_sce$var, xlab="Log Expression Mean", 
     ylab= "Log Expression Variance")
curve(dec.lmmp_mf_sce$trend(x), col="dodgerblue", add=TRUE, lwd=2)
plot(dec.noweight$mean,dec.noweight$var, xlab="Log Expression Mean", 
     ylab= "Log Expression Variance")
curve(dec.noweight$trend(x), col="red", add=TRUE, lwd=2)
legend("topleft", col=c("dodgerblue", "red") , 
       legend=c("Default","Noweight"), lwd=2)
dev.off()

fit.mf[order(fit.mf$bio, decreasing = TRUE),]

#If needed, use w/spikein

#Account for blocking factors
par(mfrow = c(1,2))
blocked.stats <- dec.lmmp_mf_sce$per.block
for (i in colnames(blocked.stats)) {
  current <- blocked.stats[[i]]
  plot(current$mean, current$total, main=i, pch=16, cex=0.5,
       xlab="Mean of log-expression", ylab="Variance of log-expression")
  curfit <- metadata(current)
  points(curfit$mean, curfit$var, col="red", pch=16)
  curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2) 
}
#Highly Variable Approach
hvg_mf_var <- getTopHVGs(dec.lmmp_mf_sce, n=1000)
str(hvg_mf_var)
