#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
nloci = as.numeric(args[2])
suppressPackageStartupMessages(library(pegas))
regions = read.loci(args[1], col.loci = 2:(1+nloci), row.names = 1)
# print(typeof(regions))

# pie <- countHap(h,data)

# print(regions,details=TRUE)
regions <- as.DNAbin(ifelse(t(regions==1),"a", "t"))

# print(regions,details=TRUE)
h <- pegas::haplotype(regions)
# print(h)

h <- sort(h, what = "label")
net <- pegas::haploNet(h)

ind.hap<-with(stack(setNames(attr(h, "index"), rownames(h))),table(hap=ind, pop=rownames(regions)[values])) #(colours by frequencies, but want to colour by presence of coding snp vs high-weighted nc snp)
pdf("hap_network.pdf")
plot(net, size=attr(net, "freq"), scale.ratio=0.2, pie=ind.hap)
legend(50,50, colnames(ind.hap), col=rainbow(ncol(ind.hap)), pch=20)
dev.off()
