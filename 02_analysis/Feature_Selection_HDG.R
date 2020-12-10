#Highly Deviant Approach， input must be un-normalized sce
#if(!require(scry)){
# BiocManager::install("scry")


library(scry)
sce.nonorm <- devianceFeatureSelection(sce.nonorm,
                         assay= "counts", batch=treatment)
plot(rowData(sce.nonorm)$binomial_deviance, type="l", xlab="ranked genes",
     ylab="binomial deviance", main="Feature Selection with Deviance")
abline(v=2000, lty=2, col="red")

#Dimension Reduction w/nullresiduals

sce.nonorm <- nullResiduals(sce.nonorm, assay="counts", type="deviance")
sce.nonorm <- nullResiduals(sce.nonorm, assay="counts", type="pearson")
sce.nonorm2 <- sce.nonorm[1:1000, ] #select top 1000 HDG
pca.nr <- function(Y, L=2, center=TRUE, scale=TRUE){
  res<-prcomp(as.matrix(t(Y)), center=center, scale.=scale, rank.=L)
  factors<-as.data.frame(res$x)
  colnames(factors)<-paste0("dim", 1:L)
  factors
}

pca_d<-pca(assay(sce.nonorm2, "binomial_deviance_residuals"))
pca_d$resid_type<-"deviance_residuals"
pca_p<-pca(assay(sce.nonorm2, "binomial_pearson_residuals"))
pca_p$resid_type<-"pearson_residuals"
cm<-as.data.frame(colData(sce2))
pd<-rbind(cbind(cm, pca_d), cbind(cm, pca_p))
ggplot(pd, aes(x=dim1, y=dim2, colour=phenoid)) + geom_point() +
  facet_wrap(~resid_type, scales="free", nrow=2) +
  ggtitle("PCA applied to null residuals of high deviance genes")

#