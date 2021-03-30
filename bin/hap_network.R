#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
nloci = as.numeric(args[2])
target = paste("X",args[3],sep='')
suppressPackageStartupMessages(library(pegas))
regions = read.table(args[1],row.names=1,header = TRUE)

# pie <- countHap(h,data)

regDNA <- as.DNAbin(ifelse(regions==1,"a", "t"))

print(regDNA,details=TRUE)
h <- pegas::haplotype(regDNA)
print(h)

h <- sort(h, what = "label")
net <- pegas::haploNet(h)

ind.hap<-with(stack(setNames(attr(h, "index"), rownames(h))),table(hap=ind, pop=rownames(regDNA)[values])) #(colours by ids in first col
# print(ind.hap[1:4,1:4])
mydata <- as.data.frame(ind.hap)
good <- mydata[mydata$Freq == 1,]
# print(head(good))
locations <- strsplit(as.character(good$pop), '_')
# print(locations)
locations <- sapply(locations, `[[`, 2)
# print(locations)
new.hap <- table(good$hap, locations)
# print(new.hap[1:4,1:4])

#adapted from iain's version:
snp.i <- which(colnames(regDNA)==target)     #which snp is the focal snp
# print(snp.i)
# print(str(regions))
# print(head(regions))
# print(regions[snp.i])
inc.i <- which(regions[,snp.i]==1)      #Which individuals have the focal snp
# print(inc.i)
hap.has.snp <- rep(FALSE, dim(h)[1])
for(i in 1:dim(h)[1]){
    if(attr(h, "index")[[i]][1] %in% inc.i){
        hap.has.snp[i] <- TRUE
    }
}

# print(new.hap[1:4,1:4])
pdf("hap_network.pdf")
# plot(net, size=log2(attr(net, "freq")+1), scale.ratio=.5, cex = 0.8,labels=FALSE,
#       show.mutation=FALSE,bg=ifelse(hap.has.snp, "red", "blue"),pie=new.hap)
plot(net, size=log2(attr(net, "freq")+1), scale.ratio=.5, cex = 0.8,labels=FALSE,
      show.mutation=FALSE,bg=terrain.colors(ncol(new.hap)),pie=new.hap)
legend("bottomright", colnames(new.hap), col=terrain.colors(ncol(new.hap)), pch=20)
dev.off()
