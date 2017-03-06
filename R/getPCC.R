getPCC <- function(A){
  # Normalization
  A.list <-
    as.list(as.data.frame(t(A)))
  A.pair.list <-
    combn(A.list, 2, simplify = F)
  # A.pair.list.Name <-
  #   lapply(A.pair.list, function(x) names(x))
  M <-
    ncol(A)
  N <-
    nrow(A)

  # 1. Noise model correlation scores ----
  rept <- 10 # Taks about 4.4 mins
  for(i in 1:rept){
    A.rpoisson <-
      apply(A, c(1, 2), function(x){rpois(1, lambda = x)})
    C.rpoisson <-
      A.rpoisson + 1/M
    B.rpoisson <-
      C.rpoisson/rowSums(C.rpoisson)
    B.cor <-
      cor(t(B.rpoisson))
    B.cor[is.na(B.cor)] <-
      0
    if( i == 1){
      PCC.mat <- B.cor
    }else{
      PCC.mat <- PCC.mat + B.cor
    }
  }
  PCC.mat.avg <-
    PCC.mat/rept
  B.PCC <-
    PCC.mat.avg[lower.tri(PCC.mat.avg, diag = FALSE)]
  return(B.PCC)
}
