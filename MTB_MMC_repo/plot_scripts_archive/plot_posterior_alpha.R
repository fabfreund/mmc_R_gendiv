files= list.files(path = "resdata/.", pattern = "unig_beta1redstat_het90")
pdf("alpha_B_1-2_spline.pdf")


names=c("Bainomugisa 2018 (2014)" , "Bainomugisa 2018" , "Bjorn-Mortensen 2016" , "Comas 2015" , "Eldholm 2015 (1998)" , "Eldholm 2015 (2001)" , "Eldholm 2015 (2003)" , "Eldholm 2015" , "Eldholm 2016" , "Folkvardsen 2017 (2009)" , "Folkvardsen 2017 (2010)" , "Folkvardsen 2017 (2012)" , "Folkvardsen 2017" , "Lee 2015 (2012 clade A)" , "Lee 2015 (2012)" , "Lee 2015 (clade A)" , "Lee 2015 (clade B)" , "Lee 2015 (clade C)" , "Lee 2015" , "Roetzer 2013" , "Shitikov 2017" , "Stucki 2015" , "Stucki 2016")

par(mfrow=c(4,2))

i=0

for (file in files){
i= i+1
par(mar=c(3,3,1.5,1),mgp=c(1.5,0.5,0))


load (paste0("resdata/",file))    
    name = gsub("res2","",file)
    name = gsub("unig_beta1redstat_het90.RData","",name)




x=c()
y=c()

for (t in seq(from = 1, to = 37, by = 4 )){
    
	
	if (t ==1){
		x=c(betafit_rf$quantiles[2])

		y=c(0)

	}else{
		if (t==37){	
			x=c(x,betafit_rf$quantiles[40])
			y=c(y,0)
		}else{
			x=c(x,betafit_rf$quantiles[t+2])
			y=c(y,(0.05)/(betafit_rf$quantiles[t+4]-betafit_rf$quantiles[t]))

		}
	}

}

#spline_int <- as.data.frame(spline(x, y,n = 1*length(x)))

if (name=="Shitikov2017"){ymax=10}else{

	if (name=="Bjorn-Mortensen2016"){ymax=15}else{

	        if (name=="Lee2015_2012_cladeA"){ymax=7.5}else{

			ymax=5
		}
	}
}


xdec=c(betafit_rf$quantiles[1],betafit_rf$quantiles[5],betafit_rf$quantiles[9],betafit_rf$quantiles[13],betafit_rf$quantiles[17],betafit_rf$quantiles[21],betafit_rf$quantiles[25],
	betafit_rf$quantiles[29],betafit_rf$quantiles[33],betafit_rf$quantiles[37],betafit_rf$quantiles[41])

plot(x,y,xlim=c(0,2),col="blue",main=names[i],type="l",ylim=c(-ymax/30,ymax),xlab="alpha",ylab="Density",yaxt='n')

rug(xdec,col="blue",ticksize=0.06)

abline(v=betafit_rf$quantiles[21],lty=3,col="blue")
abline(v=1,lty=2)
abline(h=0,lty=1)


load(paste0("resdata/res3",name,"logg_beta0_redstat_het90.RData"))

x=c()
y=c()

for (t in seq(from = 1, to = 37, by = 4 )){


        if (t ==1){
                x=c(betafit_rf$quantiles[2])

                y=c(0)

        }else{
                if (t==37){
                        x=c(x,betafit_rf$quantiles[40])
                        y=c(y,0)
                }else{
                        x=c(x,betafit_rf$quantiles[t+2])
                        y=c(y,((0.05)/(betafit_rf$quantiles[t+4]-betafit_rf$quantiles[t])))

                }
        }

}

#spline_int <- as.data.frame(spline(x, y,n = 5*length(x)))

xdec=c(betafit_rf$quantiles[1],betafit_rf$quantiles[5],betafit_rf$quantiles[9],betafit_rf$quantiles[13],betafit_rf$quantiles[17],betafit_rf$quantiles[21],betafit_rf$quantiles[25],
        betafit_rf$quantiles[29],betafit_rf$quantiles[33],betafit_rf$quantiles[37],betafit_rf$quantiles[41])   

lines(x,y,col="red")                          
rug(xdec,col="red",ticksize=0.06)
abline(v=betafit_rf$quantiles[21],lty=3,col="red")


legend("topleft",legend=c("alpha[1,2]", "alpha[0,2]"), col=c("blue","red"),lty=1, bty="n",cex=0.75)

}

dev.off()

