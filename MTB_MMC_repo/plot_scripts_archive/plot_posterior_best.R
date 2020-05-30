pdf("best_post.pdf")
par(mfrow=c(4,2))

res=NULL

table = read.table("res_unig_90")   ### # cat res/*unig*beta1*het90.txt | grep v dataset > res_unig_90

names=c("Bainomugisa 2018 (2014)" , "Bainomugisa 2018" , "Bjorn-Mortensen 2016" , "Comas 2015" , "Eldholm 2015 (1998)" , "Eldholm 2015 (2001)" , "Eldholm 2015 (2003)" , "Eldholm 2015" , "Eldholm 2016" , "Folkvardsen 2017 (2009)" , "Folkvardsen 2017 (2010)" , "Folkvardsen 2017 (2012)" , "Folkvardsen 2017" , "Lee 2015 (2012 clade A)" , "Lee 2015 (2012)" , "Lee 2015 (clade A)" , "Lee 2015 (clade B)" , "Lee 2015 (clade C)" , "Lee 2015" , "Roetzer 2013" , "Shitikov 2017" , "Stucki 2015" , "Stucki 2016")


for (i in (1:nrow(table))){

	par(mar=c(3,3,1.5,1),mgp=c(1.5,0.5,0))

	if (table[i,4] == 1){scale=5000}else{scale=1}


	for (t in seq(from = 1, to = 37, by = 4 )){
		z=t+6		# match table'scolumn

        	if (t ==1){
                	x=c(table[i,8])

               	 y=c(0)

        	}else{
                	if (t==37){
                        	x=c(x,table[i,46])
                        	y=c(y,0)
                	}else{
                        	x=c(x,table[i,z+2])
                        	y=c(y,(scale/(table[i,z+4]-table[i,z])))

	                }
       		 }

	}

	ymax=max(y)	

	if (table[i,4] == 1){

		plot(x,y,xlim=c(0,5000),main=names[i],type="l",xlab="g",ylab="Density",yaxt='n',ylim=c(-(ymax*1.1)/30,ymax*1.1))
	}

	if (table[i,4] == 2){


		plot(x,y,xlim=c(1,2),main=names[i],type="l",xlab="alpha",ylab="Density",yaxt='n',ylim=c(-(ymax*1.1)/30,ymax*1.1))

	}

	if (table[i,4] == 3){


                plot(x,y,xlim=c(0,1),main=names[i],type="l",xlab="psi",ylab="Density",yaxt='n',ylim=c(-(ymax*1.1)/30,ymax*1.1))


	}


	xdec=c(table[i,7],table[i,11],table[i,15],table[i,19],table[i,23],table[i,27],table[i,31],
                table[i,35],table[i,39],table[i,43],table[i,47])

	abline(h=0)
	abline(v=table[i,27],lty=3)
	rug(xdec,ticksize=0.065)
}


 
dev.off()

