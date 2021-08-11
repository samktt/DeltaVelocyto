## usethis namespace: start
#' @useDynLib DeltaVelocyto, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' This reads the STARsolo Velocyto folder output
#' @export
read <- function(velocytoFolder) {
  cat(crayon::red("Reading STARsolo Velocyto folder contents..."))

  X <- Matrix::readMM(paste(velocytoFolder,"/matrix.mtx", sep=""))
  s <- Matrix::readMM(paste(velocytoFolder,"/spliced.mtx", sep=""))
  u <- Matrix::readMM(paste(velocytoFolder,"/unspliced.mtx", sep=""))
  barcodes <- read.table(paste(velocytoFolder,"/barcodes.tsv", sep=""), sep = '\t', header = FALSE)
  features <- read.table(paste(velocytoFolder,"/features.tsv", sep=""), sep = '\t', header = FALSE)
  colnames(s) <- barcodes[,1]
  rownames(s) <- features[,2]
  colnames(u) <- barcodes[,1]
  rownames(u) <- features[,2]
  colnames(X) <- barcodes[,1]
  rownames(X) <- features[,2]

  ladata = list(s = s, u = u, X = X)
  return(ladata)
}

#' Select the top (most variable) n genes with a minimum of x shared counts across cells
#' @export
filter_and_norm <- function(ladata, n = 1000, min_counts = 20) {

  #cat(crayon::red("Removing genes with less than", x, "counts"))
  ladata$s = ladata$s[rowSums(as.matrix(ladata$X))>min_counts,]
  ladata$u = ladata$u[rowSums(as.matrix(ladata$X))>min_counts,]
  ladata$X = ladata$X[rowSums(as.matrix(ladata$X))>min_counts,]

  #Normalization
  ladata$s <- Seurat::NormalizeData(ladata$s,verbose = FALSE)
  ladata$u <- Seurat::NormalizeData(ladata$u,verbose = FALSE)
  ladata$X <- Seurat::NormalizeData(ladata$X, verbose = FALSE)

  #cat(crayon::red("Ranking genes based on variance and mean of the log-expression values"))
  stats <- scran::modelGeneVar(ladata$X)
  topGenes = scran::getTopHVGs(stats,  n = n)
  ladata$s = ladata$s[intersect(intersect(rownames(ladata$s), rownames(ladata$u)), topGenes),]
  ladata$u = ladata$u[intersect(intersect(rownames(ladata$s), rownames(ladata$u)), topGenes),]

  return(ladata)
}

#' Get first order moment (mean across nearest neighbors)
#' @export
moment <- function(ladata, n = 30) {
  cat(crayon::red("Computing PCA (this may take a while...)"))
  #redundunt 1:1000
  pca = DeltaVelocyto::getPCA(t(as.matrix(ladata$X))[,1:1000])
  cat(crayon::red("Computing ", n, " nearest neighbors."))
  #neighbors = Seurat::FindNeighbors(pca$x[,1:30], k.param = 30, nn.method = "annoy")
  dim(pca)
  rownames(pca) =  colnames(ladata$X)
  neighbors = Seurat::FindNeighbors(pca[,1:30], k.param = 30, nn.method = "annoy")

  NN = neighbors$nn
  diag(NN) = 0
  #n*n boolean matrix telling us the neighbors of each cell
  conn = as.matrix(NN)>0

  ladata$Ms = t(conn %*% t(as.matrix(ladata$s)))
  ladata$Mu = t(conn %*% t(as.matrix(ladata$u)))

  return(ladata)
}

#' Find the linear model of spliced and unspliced moments in a gene
#' @export
lm <- function(ladata, gene) {

  ppMs = ladata$Ms[gene, which(ladata$Ms[gene,] != 0)]
  ppMu = ladata$Mu[gene, which(ladata$Ms[gene,] != 0)]

  #Fitting a linear model
  lmGene = stats::lm(ppMu~ppMs, data =as.data.frame(cbind(ppMs, ppMu)))
  summary(lmGene)
  slope = lmGene$coefficients[2]
  plot(ppMs, ppMu)
  abline(lmGene)

  return(slope)
}


#' compare two slope of unspliced vs spliced moments in a gene of two samples
#' @export
delta <- function(ladata1, ladata2, gene) {

  #filter out counts where the spliced count is 0
  ppMs1 = ladata1$Ms[gene, which(ladata1$Ms[gene,] != 0)]
  ppMu1 = ladata1$Mu[gene, which(ladata1$Ms[gene,] != 0)]
  data1 = as.data.frame(cbind(ppMu1, ppMs1))

  ppMs2 = ladata2$Ms[gene, which(ladata2$Ms[gene,] != 0)]
  ppMu2 = ladata2$Mu[gene, which(ladata2$Ms[gene,] != 0)]
  data2 = as.data.frame(cbind(ppMu2, ppMs2))

  colnames(data1) = c("y", "x")
  colnames(data2) = c("y", "x")

  combined = rbind(data1, data2)

  A = rep("A", length(ppMs1))
  B = rep("B", length(ppMs2))
  g = c(A,B)
  d2 = cbind(combined,g)

  m2 <- lm(y ~ x * g, data = d2)
  m2.null <- lm(y ~ x + g, data = d2)
  anova(m2.null, m2)
  coef(m2)

  return(anova(m2.null, m2))
}






