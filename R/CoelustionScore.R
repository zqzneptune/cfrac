CoelustionScore <- function(A){
  library(dplyr)
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
  PCC <-
    PCC.mat.avg[lower.tri(PCC.mat.avg, diag = FALSE)]
  # cbind(t(combn(rownames(PCC.mat.avg), 2)),
  #       PCC.mat.avg[lower.tri(PCC.mat.avg, diag = FALSE)])
  # 2. Weighted Cross-Correlation ----
  library(wccsom)
  A.pair.list.WCC <-
    lapply(A.pair.list, function(x) wcc( x[[1]] , x[[2]] , trwdth = 1))
  WCC <- unlist(A.pair.list.WCC)
  WCC[is.na(WCC)] <- 0
  
  # 3. Co-APEX scores ----
  A.mxpep <-
    apply(A, 1, which.max)
  Apex_A <-
    data.frame(InteractorA = names(A.mxpep),
               ApexA = A.mxpep,
               stringsAsFactors = F)
  Apex_B <-
    data.frame(InteractorB = names(A.mxpep),
               ApexB = A.mxpep,
               stringsAsFactors = F)
  C.df <- 
    data.frame(t(combn(names(A.mxpep), 2)), stringsAsFactors = F)
  colnames(C.df) <-
    c("InteractorA", "InteractorB")
  C.df <-
    left_join(C.df, Apex_A, by = "InteractorA")
  C.df <-
    left_join(C.df, Apex_B, by = "InteractorB")
  C.df[, "Coapex"] <- 
    ifelse(C.df$ApexA == C.df$ApexB, 1, 0)
  A.apex <-
    C.df[, c("InteractorA", "InteractorB", "Coapex")]
  # A.apex <-
  #   unlist(lapply(A.pair.list, function(x){
  #     ifelse(A.mxpep[[x[[1]]]] == A.mxpep[[x[[2]]]], 1, 0)}))
  # Cbind all scores ----
  Out_score <-
    cbind(A.apex, PCC, WCC)
  return(Out_score)
}