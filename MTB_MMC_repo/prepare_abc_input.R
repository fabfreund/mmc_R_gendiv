#' Prepare ALL TB data sets for ABC analysis


library(pegas)

#' Read in data sets

#data_folder <- "data_call90"
data_folder <- "data_call75"

loc1 <- list.files(path=data_folder,pattern = "fasta")
data_names <- unlist(strsplit(loc1,".fasta"))
loc1 <- paste0(data_folder,"/",loc1)


data_main <- vector("list",length(data_names))
names(data_main) <- data_names

for (i in 1:length(data1)){
  seq1 <- read.FASTA(loc1[i])
  seq1 <- t(sapply(seq1,function(x){as.numeric(x!=seq1[[1]])}))  
  seq1 <- seq1[-1,]
  #To kick out SNPs with counts 0 or sample size
  SNP_nonseg <- (colSums(seq1)%in% c(0,nrow(seq1)))
  seq1 <- seq1[,!(SNP_nonseg)]
  data_main[[i]] <- vector("list",3)
  names(data_main[[i]]) <- c("n_ind","n_mut","seq_0_1")
  data_main[[c(i,1)]] <- nrow(seq1)
  data_main[[c(i,2)]] <- ncol(seq1)
  data_main[[c(i,3)]] <- seq1}



save(data_main,file=paste0("mtb_",data_folder,".RData"))
