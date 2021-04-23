#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
nloci = as.numeric(args[2])
target = paste("X",args[3],sep='')
suppressPackageStartupMessages(library(pegas))
regions = read.table(args[1],row.names=1,header = TRUE)

regDNA <- as.DNAbin(ifelse(regions==1,"a", "t"))

# print(regDNA,details=TRUE)
h <- pegas::haplotype(regDNA)
print(h)

h <- sort(h, what = "labels")
net <- pegas::haploNet(h)

ind.hap<-with(stack(setNames(attr(h, "index"), rownames(h))),table(hap=ind, pop=rownames(regDNA)[values])) #(colours by ids in first col
# print(ind.hap[1:4,1:4])
mydata <- as.data.frame(ind.hap)

#tie to population to colour by that
good <- mydata[mydata$Freq == 1,]
locations <- strsplit(as.character(good$pop), '_')
locations <- sapply(locations, `[[`, 2)
pop.hap <- table(good$hap, locations)

pdf("hap_network_pops.pdf")
plot(net, size=log2(attr(net, "freq")+1), scale.ratio=.5, cex = 0.8,labels=FALSE,
      show.mutation=FALSE,bg=terrain.colors(ncol(pop.hap)),pie=pop.hap)
legend("bottomright", colnames(pop.hap), col=terrain.colors(ncol(pop.hap)), pch=20)
dev.off()

#adapted from iain's version:
countHap <- function(hap = h, dna = data){
    with(
        stack(setNames(attr(hap, "index"), rownames(hap))),
        table(hap = ind, pop = attr(dna, "dimnames")[[1]][values])
    )
}

snp.pie <- countHap(h,regDNA)
snp.i <- which(colnames(regDNA)==target)     #which snp is the focal snp
inc.i <- which(regions[,snp.i]==1)      #Which individuals have the focal snp
hap.has.snp <- rep(FALSE, dim(h)[1])
for(i in 1:dim(h)[1]){
    if(attr(h, "index")[[i]][1] %in% inc.i){
        hap.has.snp[i] <- TRUE
    }
}

# print(pop.hap[1:4,1:4])
pdf("hap_network_snp.pdf")
plot(net, size=log2(attr(net, "freq")+1), scale.ratio=.5, cex = 0.8,labels=FALSE,
      show.mutation=FALSE,bg=ifelse(hap.has.snp, "red", "blue"))#,pie=snp.pie)
dev.off()
