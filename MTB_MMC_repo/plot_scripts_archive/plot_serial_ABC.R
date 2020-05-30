library(lattice)
library(gridExtra)

# Set whether you plot misclassification or parameter estimation errors
data1 <- c("modsec","pe")[2]
logg <- c("logg","contg")[1] 

datasets <- c("Eldholm2015","Lee2015","Roetzer2013")
datasets_plot <- c("Eldholm 2015","Lee 2015","Roetzer 2013")
plot_l <- vector("list",3*length(datasets))

for (i in seq(along=datasets)){
dataset1 <- datasets[i]
dataname <- datasets_plot[i]
load(paste0("resdata/serial_",data1,"_",dataset1,logg,".RData"))
if (data1=="modsec"){temp1 <- errmmc_rf_serial} else {temp1 <- err_rf_serial}
#load(paste0("resdata/serial_res/serial_",data1,"_",dataset1,"_largec.RData"))
plot_pos <- 3*(i-1)
for (j in 1:3){
main1 <- paste(dataname,": ",switch(j,"Exp. growth","Beta coalescent","Dirac coalescent"))
xlab1 <- switch(j,"Growth parameter g","Coalescent parameter alpha","Coalescent parameter p")
ylab1 <- "Rescaling factor c'"
if (data1=="modsec"){
#errmmc_rf_serial[[j]] <- rbind(temp1[[j]],errmmc_rf_serial[[j]])
if (j==1){colnames(errmmc_rf_serial[[j]])[7:8] <- c("1K","2K")}  
plot_l[[plot_pos+j]] <- levelplot(t(errmmc_rf_serial[[j]]),
                         xlab=xlab1,ylab=ylab1,main=main1,at=seq(0,1,0.1))} else {
#err_rf_serial[[j]] <- rbind(temp1[[j]],err_rf_serial[[j]])
if (j==1){colnames(err_rf_serial[[j]])[7:8] <- c("1K","2K")}
plot_l[[plot_pos+j]] <- levelplot(t(err_rf_serial[[j]]),
xlab=xlab1,ylab=ylab1,main=main1)}                           
}
}

png(paste0("pics/miscl_heatmap_",data1,logg,"_full2.png"),
    height=1700,width=switch(data1,"pe"=2000,"modsec"=2000))
orig_graph <- trellis.par.get()
trellis.par.set("layout.widths",list(left.padding=0.4,right.padding=0.4))
trellis.par.set("fontsize",list(text=24))
trellis.par.set("par.xlab.text",list(cex=1.2))
trellis.par.set("par.ylab.text",list(cex=1.2))
trellis.par.set("par.xaxis.text",list(cex=1.2))
trellis.par.set("par.yaxis.text",list(cex=1.2))
trellis.par.set("par.main.text",list(cex=1.2))                                                                            
                                                                                                                                           
grid.arrange(plot_l[[1]],plot_l[[2]],plot_l[[3]],
             plot_l[[4]],plot_l[[5]],plot_l[[6]],
             plot_l[[7]],plot_l[[8]],plot_l[[9]],
             nrow = length(datasets),ncol=3)
dev.off()


