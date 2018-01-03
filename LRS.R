LRS <- function(incl, total, group){

  coln <- length(group)
  
  inc0 <- sum(incl)
  tot0 <- sum(total)
  
  inc1 <- 0
  tot1 <- 0
  
  inc2 <- 0
  tot2 <- 0
  
  for(i in 1:coln){
    inci <- incl[i]
    toti <- total[i]
    
    if(group[i] == 1){
      inc1 <- inc1 + inci
      tot1 <- tot1 + toti
    }
    else{
      inc2 <- inc2 + inci
      tot2 <- tot2 + toti    
    }

  }
  
  p0 <- inc0/tot0
  p1 <- inc1/tot1
  p2 <- inc2/tot2
  
  llreduced <- 0
  llfull <- 0
  
  for(j in 1:coln){
    incj <- incl[j]
    totj <- total[j]
    llreduced <- llreduced + log(dbinom(incj,totj,p0))
    
    if(group[j] == 1){
      llfull <- llfull + log(dbinom(incj,totj,p1))
    }
    else{
      llfull <- llfull + log(dbinom(incj,totj,p2))   
    }
  }
  
  lrs <- -2 * (llreduced - llfull)
  
  pvalue <- pchisq(lrs, df=1, lower.tail=FALSE)
  
  result <- list(p0 = p0, p1 = p1, p2=p2, llreduced = llreduced, llfull = llfull, lrs = lrs, pvalue = pvalue)
  return(result)
}

diff.splicing <- function(psi.files, group){
  
  filecount <- length(psi.files)

  file <- read.table(psi.files[1],header=TRUE, sep="\t")
  incl <- matrix(file$num_incl_reads)
  total <- matrix(file$num_total_reads)
  output <- data.frame(gene = file$gene,exon = file$exon)
  
  exonCount <- length(incl)
  
  for(i in 2:filecount){
    file <- read.table(psi.files[i],header=TRUE, sep="\t")
    incl <- cbind(incl, file$num_incl_reads)
    total <- cbind(total, file$num_total_reads)
  }
  
  lrs <- data.frame(LRS(incl[1,],total[1,],group))
  
  for(i in 2:exonCount){
    lrs <- rbind(lrs,LRS(incl[i,],total[i,],group))
  }
  
  output <- cbind(output,lrs)
  
  pValue <- lrs$pvalue
  padj <- p.adjust(pValue)
  output <- cbind(output,padj)

  return(output)
}