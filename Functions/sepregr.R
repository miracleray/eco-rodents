# Originally written by Dr Gabriele Sansalone
# Comments and formatting by Ariel Marcy

# What format is expected of x, y, and especially of group??

# What is the output? 

sepregr <- function(x, y, group, scalevariable = 1, conf = 0, nperm = 999, labels = NULL, reorder = TRUE){
   
# Load necessary packages     
library(MASS)
library(car)
library(vegan)
library(CCA)
library(rgl)
library(calibrate)
library(compositions)

# Handle incorrrectly formatted data with warnings
if(!is.factor(group)) stop("'group' must be a factor")
if(!is.vector(x)) stop("'x' must be a vector")
if(is.null(names(x))) stop("'x' must have names")
  
# Redefine the variable "group" to one with only unique classifiers
if(reorder == TRUE){
        group <- factor(group, levels = unique(group))
}
  
# Set up for analyses
y <- as.matrix(y)
fat_species <- group
species <- as.numeric(group)
  
# Create "dati" 
dati <- cbind(y, x, species)
rownames(dati) <- rownames(y)
taglie <- tapply(x, species, max)
tagliemin <- tapply(x, species, min)

# Permutational Multivariate Analysis of Variance w/ distance matrices 
total <- adonis(y ~ x, method = "euclidean", permutations = nperm)

# Calculate CCA scores
ccascores <- cbind(rcc(as.matrix(x), y, 0.1, 0.1)$scores$xscores, rcc(as.matrix(x), y, 0.1, 0.1)$scores$yscores)
model.fullccascores <- lm(ccascores[, 2] ~ ccascores[, 1] * group)
mat.predccascores <- predict(model.fullccascores)
  
# Set up new variables for For loop
residplusmax <- NULL
obsatmin <- NULL
obsatmax <- NULL
mat.resid <- NULL
mat.pred <- NULL
mat.coef <- NULL
mat.test <- matrix(ncol = 2, nrow = max(species))

# For loop subsetting 
for(i in 1:(max(species))) {
    dat <- as.matrix(subset(as.data.frame(dati), species == i))
    
    Y <- as.matrix(dat[, 1:ncol(y)])
    X <- dat[, (ncol(y) + 1):(ncol(dati) - 1)]
    names(X) <- names(x)[species == i]
    spe <- dat[, ncol(dati)]
    vetti <- c(1, taglie[i])
    vettimin <- c(1, tagliemin[i])
    model.fulli <- lm(Y[spe == i, ] ~ X[spe == i])
    full_maxi <- as.numeric(crossprod(coef(model.fulli), vetti))
    full_mini <- as.numeric(crossprod(coef(model.fulli), vettimin))
    ado <- adonis(Y[spe == i, ] ~ X[spe == i], method = "euclidean", permutations = nperm)
    
    residi <- as.matrix(resid(model.fulli))
    rownames(residi) <- names(X)[spe == i]
    
    predicti <- as.matrix(predict(model.fulli))
    rownames(predicti) <- names(X)[spe == i]
    
    # New section...
    mat.test[i, ] <- c(ado$aov.tab$R2[1], ado$aov.tab$Pr[1])
    mat.resid <- rbind(mat.resid, residi)
    mat.pred <- rbind(mat.pred, predicti)
    mat.coef <- rbind(mat.coef, coef(model.fulli))
    obsatmax <- rbind(obsatmax, full_maxi)
    obsatmin <- rbind(obsatmin, full_mini) 
  }
 
# New section...
  mat.residfin <- as.matrix(mat.resid[match(names(x), rownames(mat.resid)), ])
  rownames(mat.residfin) <- names(x)
  mat.predfin <- as.matrix(mat.pred[match(names(x), rownames(mat.pred)), ])
  rownames(mat.predfin) <- names(x)
  
  PCsonpred <- prcomp(mat.predfin)$x
  PCsony <- prcomp(y)$x
  varpcsonpred <- summary(prcomp(mat.predfin))$importance
  colnames(mat.test) <- c("R2","p.value")
  rownames(mat.test) <- levels(fat_species)
  rownames(obsatmax) <- levels(fat_species)
  rownames(obsatmin) <- levels(fat_species)
  ttnumber <- as.numeric(group)
  tt <- tabulate(ttnumber)
  ff <- obsatmax[rep(seq(nrow(obsatmax)), tt), ]
  residplusmax <- mat.resid + ff
  residplusmaxfin <- as.matrix(residplusmax[match(names(x), rownames(residplusmax)), ])
  rownames(residplusmaxfin) <- names(x)
  
  ffmin <- obsatmin[rep(seq(nrow(obsatmin)), tt), ]
  residplusmin <- mat.resid + ffmin
  residplusminfin <- as.matrix(residplusmin[match(names(x), rownames(residplusmin)), ])
  rownames(residplusminfin) <- names(x)
  
# New section...
  pooled_correlation <- as.matrix(c(total$aov.tab$R2[1], total$aov.tab$Pr[1]))
  rownames(pooled_correlation) <- c("R-sq","p_value")
  
 # New section...
  x11()
  plot(x,y[,1],col="white")
  legend(min(x), max(y[, 1]), levels(group), cex = 1, col = c(1:nlevels(group)), pch = 20, box.col = "white")
  
# The output of the function
out <- (list(x = x, 
             y = y, 
             group = group, 
             pooled_correlation = pooled_correlation, 
             regression_results = mat.test, 
             PCsony = PCsony, 
             predicted = mat.predfin, 
             residuals = mat.residfin,
             obsatmin = obsatmin, 
             obsatmax = obsatmax,
             residplusmax = residplusmaxfin,
             residplusmin = residplusminfin,
             coefficients = mat.coef,
             PCsonpred = PCsonpred,
             varpcsonpred = varpcsonpred,
             ccascores = ccascores,
             predccascores = mat.predccascores,
        
             plot = if(ncol(PCsonpred) >= 3) { 
                        plot3D(PCsonpred[, 1:3],
                        col = as.numeric(group),
                        bbox = FALSE,
                        type = "s",
                        cex = scalevariable) 
                     } else {NULL}, 
             
             axis1 = if(ncol(PCsonpred) >= 3) { 
                     arrows3D(c(max(PCsonpred[, 1:3]), 0, 0), c(-max(PCsonpred[, 1:3]), 0, 0), length = 0.03, lwd = 100, add = TRUE, labs= "1PConpred") 
                     } else {NULL} ,
             
             axis2 = if(ncol(PCsonpred) >= 3) {
                     arrows3D(c(0, max(PCsonpred[, 1:3]), 0), c(0, -max(PCsonpred[,1:3]), 0), length = 0.03, add = TRUE, labs = "2PConpred")
                     } else {NULL},
             
             axis3 = if(ncol(PCsonpred) >= 3) {
                     arrows3D(c(0, 0, max(PCsonpred[, 1:3])), c(0, 0, -max(PCsonpred[, 1:3])), length = 0.03, add = TRUE, labs = "3PConpred")
                     } else {NULL},
             
             x11 = x11(),
             plot2 = plot(x,ccascores[, 2], cex = scalevariable, col = as.numeric(group), ylim = range(ccascores[, 2])), par = par(new=T),
             
             plot3 = plot(x, mat.predccascores, cex = scalevariable, col = as.numeric(group), pch = 20, ylim = range(ccascores[, 2])),  
             #newscene=open3d(), # WHY IS THIS COMMENTED OUT?
             
             plot6 = if(ncol(PCsonpred) >= 2) {
                     plot3dcol(cbind(x, PCsonpred[, 1:2]), group, scalevariable, reorder = reorder, legend = F, axlab = c("x", "1PConpred", "2PConpred"))
                     } else {NULL},  
             
             ellipses = for(i in 1:max(species)) {
               par <- par(new = TRUE) 
               dataEllipse(x[species == i], ccascores[, 2][species == i], col = i, add = TRUE, levels = c(conf), center.pch = 18, center.cex = scalevariable * 2, plot.points = FALSE, fill = TRUE)
             },
             x11 = x11(),
             plot5 <- plot(x, lm(mat.predccascores ~ x * group)$fitted.values, col = as.numeric(group), pch = 20, cex = scalevariable), 
             tt = if(!is.null(labels)) {
                     textxy(x, lm(mat.predccascores ~ x * group)$fitted.values, labels) 
                     } else {NULL}
  ))
  return(out)
}