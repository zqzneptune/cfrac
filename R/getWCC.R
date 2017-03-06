getWCC <- function(A){
  # Normalization
  A.list <-
    as.list(as.data.frame(t(A)))
  # 2. Weighted Cross-Correlation ----
  library(wccsom)
  A.pair.list.WCC <-
    lapply(A.pair.list, function(x)wcc(x[[1]], x[[2]], trwdth = 1))
  A.WCC <- unlist(A.pair.list.WCC)
  A.WCC[is.na(A.WCC)] <- 0
  return(A.WCC)
}
