library("RColorBrewer")
library("ape")

setwd("~/Documents/ExploringTEsinCannabisGenome/genome_all/")

repFamilies <- read.table("repeat_content_4genomes_filtered.csv", sep = ",", header = T)
rownames(repFamilies) <- repFamilies$X

CommonTree <- read.tree("cladogram_4genomes_NCBI.phy")
EudicotRepFamilies <- repFamilies[,CommonTree$tip.label]

families <- rownames(EudicotRepFamilies)
families[1] <- "Simple repeat"
colors <- brewer.pal(length(families),"Set3")
cex.plot <- 1.2


pdf("~/Google Drive/Pisupati et al. 2017/version2/Figure1.pdf")
layout(matrix(c(1,1,1,1,1,1,1,2,2,3,3), 11, 1, byrow = TRUE))
par(mar=c(2, 5, 1.2, 0.8))
barplot(as.matrix(EudicotRepFamilies), ylim = c(0,76),width = 10,col = colors, border = T, las = 1, space = 0.54, xlim = c(0,80), ylab = "Percentage of genome (%)", cex.axis = cex.plot, cex.names = cex.plot, cex.lab = cex.plot )
par(mar=c(2.4, 11, 1.3, 5.5))
plot(CommonTree, cex = 0.9,direction = "upwards", use.edge.length = F, node.depth = 2, show.tip.label = F)
plot.new()
legend("center", families, fill = colors, ncol=3, bty = "n", cex = cex.plot, pt.cex = 1, text.width = rep(0.2,length(families)), x.intersp=0.5,xjust=0,yjust=0)
dev.off()




