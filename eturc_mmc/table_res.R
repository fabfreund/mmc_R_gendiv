full_files <- list.files(path = "res",pattern=".txt",
                         full.names = TRUE)

res_report <- NULL

for (name1 in full_files){
res_report <- rbind(res_report,read.table(name1,header = TRUE))  
}



write.table(res_report,file = "res_m123_table.txt",
            row.names = FALSE,quote = FALSE)


