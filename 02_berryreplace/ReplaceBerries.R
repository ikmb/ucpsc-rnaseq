library("data.table")
berrytable<-fread("../00_RawData/results/berries.txt")
berrytable[,berry:=gsub("berry","",berry)]
fwrite(berrytable,file="local/replaced_berries_result.txt")
