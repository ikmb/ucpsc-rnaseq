library("data.table")
berrytable<-fread("../00_RawData/berries.txt")
berrysum<-sum(berrytable$number)
write(berrysum,file="results/berrysum_result.txt")
