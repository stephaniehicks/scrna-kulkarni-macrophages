suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scRNAseq)
  library(scater)
  library(flexmix)
  library(splines)
  library(BiocParallel)
  library(miQC)
})
metrics <- as.data.frame(colData(sce1))
p <- ggplot(metrics,aes(x=detected,y=subsets_Mito_percent)) + geom_point()
p


mixtureModel <- function(sce1, model_type = "linear") {
  metrics <- as.data.frame(colData(sce1))
  
  if (model_type == "linear") {
    model <- flexmix(subsets_Mito_percent~detected,
                     data = metrics, k = 2)
  } else if (model_type == "spline") {
    model <- flexmix(subsets_Mito_percent~bs(detected),
                     data = metrics, k = 2)
  } else if (model_type == "polynomial") {
    model <- flexmix(subsets_Mito_percent~poly(detected, degree = 2),
                     data = metrics, k = 2)
  }
  
  model
}

model <- mixtureModel(sce1)

parameters(model)

head(posterior(model))

plotModel(sce1, model)

plotModel <- function(sce1, model = NULL, detected = "detected",
                      subsets_Mito_percent = "subsets_Mito_percent") {
  metrics <- as.data.frame(colData(sce1))
  
  if(is.null(model)) {
    warning("call 'mixtureModel' explicitly to get stable model features")
    model <- mixtureModel(sce1)
  }
  
  predictions <- fitted(model)
  Comp.1 <- predictions[,1]
  Comp.2 <- predictions[,2]
  fitted_models <- as.data.frame(cbind(detected = metrics$detected,
                                       Comp.1 = Comp.1,
                                       Comp.2 = Comp.2))
  
  
  intercept1 <- parameters(model, component = 1)[1]
  intercept2 <- parameters(model, component = 2)[1]
  if (intercept1 > intercept2) {
    compromised_dist <- 1
  } else {
    compromised_dist <- 2
  }
  
  post <- posterior(model)
  prob_compromised <- post[, compromised_dist]
  metrics <- cbind(metrics, prob_compromised = prob_compromised)
  
  p <- ggplot(metrics, aes(x = detected, y = subsets_Mito_percent,
                           colour = prob_compromised)) +
    labs(x = "Unique genes found", y = "Percent reads mitochondrial",
         color = "Probability\ncompromised") +
    geom_point() +
    geom_line(data = fitted_models, inherit.aes = FALSE,
              aes(x = detected, y = Comp.1), lwd = 2) +
    geom_line(data = fitted_models, inherit.aes = FALSE,
              aes(x = detected, y = Comp.2), lwd = 2) +
    ylim(0, max(metrics$subsets_Mito_percent))
  
  p
}


plotFiltering <- function(sce1, model = NULL, posterior_cutoff = 0.75,
                          palette = c("#999999", "#E69F00"),
                          detected = "detected",
                          subsets_Mito_percent = "subsets_Mito_percent") {
  metrics <- as.data.frame(colData(sce1))
  
  if(is.null(model)) {
    warning("call 'mixtureModel' explicitly to get stable model features")
    model <- mixtureModel(sce1)
  }
  
  intercept1 <- parameters(model, component = 1)[1]
  intercept2 <- parameters(model, component = 2)[1]
  if (intercept1 > intercept2) {
    compromised_dist <- 1
  } else {
    compromised_dist <- 2
  }
  
  post <- posterior(model)
  prob_compromised <- post[, compromised_dist]
  keep <- prob_compromised <= posterior_cutoff
  
  metrics <- cbind(metrics, prob_compromised = prob_compromised, keep = keep)
  
  p <- ggplot(metrics, aes(x = detected, y = subsets_Mito_percent,
                           colour = keep)) +
    labs(x = "Unique genes found", y = "Percent reads mitochondrial",
         color = "Keep") +
    scale_color_manual(values = palette) + geom_point()
  
  p
}
plotFiltering(sce1, model)
