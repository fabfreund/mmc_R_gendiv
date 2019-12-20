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
    if (names(target2a)[j]=="Taj.'s D"){noneg_stat <- FALSE}
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
load("ppc_mod1.RData")

names1 <- c("Bainomug. 2018","Bainomug. 2018","Bjorn-Mort. 2016",
            "Comas 2015","Eldholm 2016",rep("Eldholm 2015",4),
            rep("Folkvardsen 2017",4),rep("Lee 2015",5),
            "Roetzer 2013","Stucki 2015",rep("Stucki 2016",3),"Shitikov 2017")
names2 <- c("","sampled in 2014","","","","",
            paste("sampled in",c("2001","2003","1998")),"",
            paste("sampled in",c("2010","2012","2009")),
            "Clade A","sampled in 2012","Clade B","","Clade C","",
            "","","Clade A","Clade B","")
           
                                                    
  
pdf("ppc_fit1.pdf")
for (i in seq(along=ppc_sims_fit1)){ 
  hist2(ppc_sims_fit1[[i]],target1[[i]],main_a = names1[i],main_b = names2[i])
}
dev.off()

 






