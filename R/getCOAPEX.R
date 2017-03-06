getCOAPEX <- function(A){

  # 3. Co-APEX scores ----
  A.mxpep <-
    as.list(apply(A, 1, which.max))
  A.apex <-
    unlist(lapply(A.pair.list.Name, function(x){
      ifelse(A.mxpep[[x[[1]]]] == A.mxpep[[x[[2]]]], 1, 0)}))
  return(A.apex)
}
