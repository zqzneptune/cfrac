runScoring <- function(Raw){
  #### Remove Sum cnt equals == 0 ###
  Dat <- 
    as.data.frame(Raw[, -c(1)])
  rownames(Dat) <- 
    as.character(unlist(Raw[, 1]))
  
  ID_rv <- # Remove proteins appearing less than two fractions with non-zero peptide counts
    rownames(Dat)[as.vector(apply(Dat, 1, function(x){sum(x==0)}) >= ncol(Dat)-1)]
  if(length(ID_rv) != 0){
    print(paste("Removing ", 
                length(ID_rv), 
                " proteins with appear in less than two fractions ...",
                sep = ""))
    Dat <-
      Dat[!(rownames(Dat) %in% ID_rv), ]
  }
  
  Sc_dat <-
    CoelustionScore(Dat)
  s <- t(apply(Sc_dat[, c("InteractorA", "InteractorB")], 1, sort))
  Sc_dat[, "PPI"] <- 
    paste(s[, 1], s[, 2], sep = "~"); rm(s)
  Out_dat <-
    Sc_dat[, c("PPI", "Coapex", "PCC", "WCC")]
  Out_dat <-
    separate(Out_dat, `PPI`, sep = "~",
             into = c("InteractorA", "InteractorB"),
             remove = F)
  return(list(`Proteins` = rownames(Dat),
              `RvPrt` = ID_rv,
              `ScoredPPI` = Out_dat))
}