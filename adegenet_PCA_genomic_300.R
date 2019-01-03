## Load vcf, run PCA, calculate K-means clustering, calculate matrix of geneic distances, run AMOVA
## Filip Kolar 2017

library(adegenet)
library(adegraphics)
library(vcfR)
library(pegas)
library(StAMPP)
library(ade4)
library(MASS)


# IMPORT SNP data from VCF
vcf <- read.vcfR("file.vcf.gz")   #read in all data

# convert to genlight 
aa.genlight <- vcfR2genlight.tetra(vcf)                           ## use the modified function vcfR2genlight.tetra at the end of the file
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")   # add real SNP.names
pop(aa.genlight)<-substr(indNames(aa.genlight),1,5)               # add pop names: here pop names are first 3 chars of ind name

#check
aa.genlight
indNames(aa.genlight)
ploidy(aa.genlight)


############
#   PCA 
############
# run PCA
pca.1 <- glPcaFast(aa.genlight, nf=300) # use the modified function glPcaFast at the end of the file

# proportion of explained variance by first three axes
pca.1$eig[1]/sum(pca.1$eig) # proportion of variation explained by 1st axis
pca.1$eig[2]/sum(pca.1$eig) # proportion of variation explained by 2nd axis
pca.1$eig[3]/sum(pca.1$eig) # proportion of variation explained by 3rd axis


# just to see pops coloured in a palette
col <- funky(10)
s.class(pca.1$scores, pop(aa.genlight),  xax=1, yax=2, col=transp(col,.6), 
        ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F)

# save nice figs
pdf ("PCA_all_SNPs_ax12.pdf", width=14, height=7)
g1 <- s.class(pca.1$scores, pop(aa.genlight),  xax=1, yax=2, col=transp(col,.6), 
              ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F, plot = FALSE)
g2 <- s.label (pca.1$scores, xax=1, yax=2, ppoints.col = "red", plabels = list(box = list(draw = FALSE), 
                                                                               optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
ADEgS(c(g1, g2), layout = c(1, 2))
dev.off()

#ploidy - differentiated plots
pdf ("PCA_all_ploidycol_SNPs_ax12.pdf", width=14, height=7)
g1 <- s.class(pca.1$scores, as.factor(as.vector(ploidy(aa.genlight))), xax=1, yax=2, col=transp(c("#FF0000", "#0000FF")), 
              ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F, plab.cex = 0 , plot = FALSE)
g2 <- s.label (pca.1$scores, xax=1, yax=2, ppoints.col = "red", plabels = list(box = list(draw = FALSE), 
                                                                               optim = TRUE), paxes.draw=T, pgrid.draw =F, plabels.cex=1, plot = FALSE)
ADEgS(c(g1, g2), layout = c(1, 2))
dev.off()

############################
#     K-means clustering
############################
## find clusters = K-means clustering
grp <- find.clusters(aa.genlight, max.n.clust=20, glPca = pca.1, perc.pca = 100, n.iter=1e6, n.start=1000)

# save the grouping
write.table(grp$grp, file="grouping_all.txt", sep="\t", quote=F, col.names=F)


###################################################
####   distance-based analyses  #######
###################################################

### Calculate Nei's distances between individuals/pops
# ---------------------------------------------------
aa.D.ind <- stamppNeisD(aa.genlight, pop = FALSE)  # Nei's 1972 distance between indivs
stamppPhylip(aa.D.ind, file="aa.indiv_Neis_distance.phy.dst") # export matrix - for SplitsTree
aa.D.pop <- stamppNeisD(aa.genlight, pop = TRUE)   # Nei's 1972 distance between pops
stamppPhylip(aa.D.pop, file="aa.pops_Neis_distance.phy.dst") # export matrix - for SplitsTree


### create the dist objects
colnames(aa.D.ind) <- rownames(aa.D.ind) 
aa.D.ind.dist<-as.dist(aa.D.ind, diag=T)
attr(aa.D.ind.dist, "Labels")<-rownames(aa.D.ind)          # name the rows of a matrix  

colnames(aa.D.pop) <- rownames(aa.D.pop) 
aa.D.pop.dist<-as.dist(aa.D.pop, diag=T)
attr(aa.D.pop.dist, "Labels")<-rownames(aa.D.pop)          # name the rows of a matrix  


### Isolation by distance 
#-----------------------
#coords file = two-column file with column names of coords belonging to each pop (lat	lon)
coords <- read.csv ("pop_coords.txt", sep ="\t")  
Dgeo <- dist(coords)

#test IBD
IBD <- mantel.randtest(Dgeo,aa.D.pop.dist)
IBD
plot(Dgeo,aa.D.pop.dist, pch=20,cex=.5)
abline(lm(aa.D.pop.dist~Dgeo))

#plot and check for denser areas in the plot indicating sub-groups
dens <- kde2d(Dgeo,aa.D.pop.dist, n=300, lims=c(-1, 22, 0, 0.04))
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, aa.D.pop.dist, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(aa.D.pop.dist~Dgeo))
title("Correlation of Genetic and Geograpgical distances")


### AMOVA
#---------
pops<-as.factor(substr(rownames(aa.D.ind),1,3))
(res <- pegas::amova(aa.D.ind.dist ~ pops))  # default nperm=1000






##################### 
# MODIFIED FUNCTIONS

# a function for conversion from vcfR object to genlight in tetraploids
vcfR2genlight.tetra <- function (x, n.cores = 1) 
{
  bi <- is.biallelic(x)
  if (sum(!bi) > 0) {
    msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
    msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
    msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object."))
    warning(msg)
    x <- x[bi, ]
  }
  x <- addID(x)
  CHROM <- x@fix[, "CHROM"]
  POS <- x@fix[, "POS"]
  ID <- x@fix[, "ID"]
  x <- extract.gt(x)
  x[x == "0|0"] <- 0
  x[x == "0|1"] <- 1
  x[x == "1|0"] <- 1
  x[x == "1|1"] <- 2
  x[x == "0/0"] <- 0
  x[x == "0/1"] <- 1
  x[x == "1/0"] <- 1
  x[x == "1/1"] <- 2
  x[x == "1/1/1/1"] <- 4
  x[x == "0/1/1/1"] <- 3
  x[x == "0/0/1/1"] <- 2
  x[x == "0/0/0/1"] <- 1
  x[x == "0/0/0/0"] <- 0
  if (requireNamespace("adegenet")) {
    x <- new("genlight", t(x), n.cores = n.cores)
  }
  else {
    warning("adegenet not installed")
  }
  adegenet::chromosome(x) <- CHROM
  adegenet::position(x) <- POS
  adegenet::locNames(x) <- ID
  return(x)
}

## -------------------------------------------------------
### a patch for MUCH MUCH faster PCA calculation on genlight objects
# see https://github.com/thibautjombart/adegenet/pull/150
glPcaFast <- function(x,
                      center=TRUE,
                      scale=FALSE,
                      nf=NULL,
                      loadings=TRUE,
                      alleleAsUnit=FALSE,
                      returnDotProd=FALSE){
  
  if(!inherits(x, "genlight")) stop("x is not a genlight object")
  # keep the original mean / var code, as it's used further down
  # and has some NA checks..
  if(center) {
    vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }
  if(scale){
    vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }
  # convert to full data, try to keep the NA handling as similar
  # to the original as possible
  # - dividing by ploidy keeps the NAs
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x)
  # handle NAs
  NAidx <- which(is.na(mx), arr.ind = T)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[,2]]
  } else {
    mx[NAidx] <- 0
  }
  # center and scale
  mx <- scale(mx,
              center = if (center) vecMeans else F,
              scale = if (scale) vecVar else F)
  # all dot products at once using underlying BLAS
  # to support thousands of samples, this could be
  # replaced by 'Truncated SVD', but it would require more changes
  # in the code around
  allProd <- tcrossprod(mx) / nInd(x) # assume uniform weights
  ## PERFORM THE ANALYSIS ##
  ## eigenanalysis
  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]
  ## scan nb of axes retained
  if(is.null(nf)){
    barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }
  ## rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  ##res$matprod <- allProd # for debugging
  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")
  ## GET LOADINGS ##
  ## need to decompose X^TDV into a sum of n matrices of dim p*r
  ## but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    myPloidy <- ploidy(x)
    for(k in 1:nInd(x)){
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if(center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }
      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
    }
    res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }
  ## FORMAT OUTPUT ##
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }
  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }
  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }
  res$call <- match.call()
  class(res) <- "glPca"
  return(res)
}

# ---------------------------------------------------------



