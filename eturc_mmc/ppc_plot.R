#' Plot posterior predictive checks for the best fitting genealogy model
#' Using scaled mutation rate theta equaling the generalized Watterson estimator
#' This is our script for performing ALL PPCs at once,
#' for own runs modify names1, names2 to match the PPC simulations performed.

#' To plot histograms, also if there is only a single value   
hist2 <- function(sims2,target2a,main_a,main_b){
  par(mfrow=c(4,5))
  for (j in 1:length(target2a)){
    main1 <- ""
    if (j==1){main1 <- main_a}
    if (j==2){main1 <- main_b}
    xrange <- range(c(sims2[j,],target2a[j]),na.rm=TRUE,finite=TRUE)
    noneg_stat <- TRUE
    if (names(target2a)[j]=="Tajima's D"){noneg_stat <- FALSE}
    plot_range <- function(p,noneg=FALSE){
              out1 <- xrange + c(-1,1)*p*mean(abs(xrange))
              if(noneg){out1 <- sapply(out1,function(x){max(0,x)})}
              return(out1)}
    if (sd(sims2[j,],na.rm = TRUE)>0){
      hist(sims2[j,],xlab=rownames(sims2)[j],main=main1,breaks = 20,
           xlim = plot_range(0.2,noneg_stat),freq=FALSE)} else {
             bar1 <- mean(sims2[j,],na.rm=TRUE)
             range1 <- plot_range(0.2,noneg_stat)
             bar2 <- abs(range1[2]-range1[1])
             plot(c(bar1,bar1),
                  c(0,1),xlim=range1,ylim=c(0,1.1),type = "l",
                  xlab = rownames(sims2)[j],
                  ylab = "Density",main=main1,bty="n")
             #box(col = "white")
             #hist(sims2[j,],xlim=plot_range(5,noneg_stat),
             #    ylim=c(0,1),breaks=1,
             #    xlab = rownames(sims2)[j],
             #     ylab = "Density",main=main1)
             
           }
    abline(v=target2a[j],col="red",lty=4)
  }
}

#' Output of ppc_abcrf_mod1.R
load("ppc_eturc_TajD_rev2.RData")



give_quant <- function(pop1,err=.05){flags <- NULL
                                 for (scaff1 in 1:5){
                                  dat1 <- ppc_sims_fit1[[pop1]][[scaff1]]
                                  obs1 <- target1[[pop1]][[scaff1]]
                                  for (i in seq(along=obs1)){ 
                  if (ecdf(dat1[i,])(obs1[i])>1-err/2 | ecdf(dat1[i,])(obs1[i])<err/2){
                              flags <- c(flags,paste("Scaff",scaff1,names(obs1)[i],
                                                     ecdf(dat1[i,])(obs1[i])))
                  }
                                  }}
                       cat(flags,"\n")
                                  }
                                                    
  
pdf("ppc_eturc_rev2.pdf")
#for (lineage in c("green","red","lightblue","pink","kenya")){ 
for (lineage in c("green","red","lightblue","kenya")){ 
  pop1 <- switch(lineage,"green"="Small Clonal","kenya"="Kenya","lightblue"="French Clonal",
                 "pink"="Diverse","red"="Big Clonal")
  for (scaff in 1:5){
  hist2(ppc_sims_fit1[[lineage]][[scaff]],target1[[lineage]][[scaff]],
        main_a = pop1,main_b = paste("Scaffold",scaff))
}}
dev.off()

load("ppc_eturc_TajD_rev_mod2.RData")

pdf("ppc_eturc_revmod2.pdf")
for (lineage in c("green","red","lightblue")){ 
  pop1 <- switch(lineage,"green"="Small Clonal","kenya"="Kenya","lightblue"="French Clonal",
                 "pink"="Diverse","red"="Big Clonal")
  for (scaff in 1:5){
    hist2(ppc_sims_fit1[[lineage]][[scaff]],target1[[lineage]][[scaff]],
          main_a = pop1,main_b = paste("Scaffold",scaff))
  }}
dev.off()
