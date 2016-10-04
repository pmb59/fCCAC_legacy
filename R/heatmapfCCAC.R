heatmapfCCAC <- function(fc){
  
  M <- c() #store sample names
  for (i in 1:length(fc$samples) ){
    M<- c(M, strsplit( as.character(fc$samples[i]), split="_vs_")[[1]] )
  }
  M <- unique(M)
  
#plot heatmap
  Fv <- matrix(NA, nrow=length(M), ncol=length(M))
  for (i in 1:length(M)){
    for (j in 1:length(M)){
      o1 <- which(as.character(fc$samples) == paste(M[i], M[j], sep="_vs_" ) )
      o2 <- which(as.character(fc$samples) == paste(M[j], M[i], sep="_vs_" ) )
      if (length(o1)==1){ Fv[i,j] <- fc[o1,2] }  #fc$F[o1]
      if (length(o2)==1){ Fv[i,j] <- fc[o2,2] }  #fc$F[o2]
      if (length(o1)==0 & length(o2)==0) {Fv[j,i] <- 100 }
    }
  }
  colnames(Fv) <- M
  rownames(Fv) <- M


  Heatmap(Fv, name='F (%)')



}
