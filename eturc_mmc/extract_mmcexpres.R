kenya_res <- NULL
BC_res <- NULL
for (i in 1:5){load(paste0("resdata/abcout_kenya",i,"_mmcg.RData"))
   temp <- c("selmodel"=as.numeric(levels(
     data2model$allocation)[data2model$allocation]),
             "postprob"=round(data2model$post.prob,3),
     "meanoob"=round(mean(conf.matrix[,4]),3))
  kenya_res <- rbind(kenya_res,temp) 
  load(paste0("resdata/abcout_red",i,"_mmcg.RData"))
  temp <- c("selmodel"=as.numeric(levels(
    data2model$allocation)[data2model$allocation]),
    "postprob"=round(data2model$post.prob,3),
    "meanoob"=round(mean(conf.matrix[,4]),3))
  BC_res <- rbind(BC_res,temp) }