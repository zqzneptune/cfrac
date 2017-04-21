runScoring <- function(Raw){
  Dat <- 
    as.data.frame(Raw[, -c(1)])
  rownames(Dat) <- 
    as.character(unlist(Raw[, 1]))
  protein_rs <- 
    rowSums(Dat)
  DatID <- # Remove proteins that have Zero peptide counts
    names(protein_rs[protein_rs != 0])
  Dat <-
    Dat[rownames(Dat) %in% DatID, ]
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
  return(list(`Proteins` = DatID,
              `ScoredPPI` = Out_dat))
}