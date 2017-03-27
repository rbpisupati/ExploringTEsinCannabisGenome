## ANOVA between the different genomes.
library(RColorBrewer)
setwd("~/Documents/ExploringTEsinCannabisGenome/")
#setwd("~/mygit/ExploringTEsinCannabisGenome/")
filter.rep <- "LTR"
## Environment files are from a new project
mor_tes <- read.csv("./genome_mor/TEfamilies_CLwise_final_table.csv", as.is = T)
#mor_tes_sub <- subset(mor_tes, mor_tes$Estimated.Half.life > 0 & mor_tes$Estimated.Half.life < 1)
mor_tes_sub <- subset(mor_tes, mor_tes$Estimated.Half.life > 0 & grepl(filter.rep, mor_tes$Top.Hit))
hum_tes <- read.csv("./genome_hum/TEfamilies_CLwise_final_table.csv", as.is = T)
#hum_tes_sub <- subset(hum_tes, hum_tes$Estimated.Half.life > 0 & hum_tes$Estimated.Half.life < 1)
hum_tes_sub <- subset(hum_tes, hum_tes$Estimated.Half.life > 0 & grepl(filter.rep, hum_tes$Top.Hit))
uso_tes <- read.csv("./genome_uso/TEfamilies_CLwise_final_table.csv", as.is = T)
#uso_tes_sub <- subset(uso_tes, uso_tes$Estimated.Half.life > 0 & uso_tes$Estimated.Half.life < 1)
uso_tes_sub <- subset(uso_tes, uso_tes$Estimated.Half.life > 0 & grepl(filter.rep, uso_tes$Top.Hit))
pur_tes <- read.csv("./genome_pk/TEfamilies_CLwise_final_table.csv", as.is = T)
#pur_tes_sub <- subset(pur_tes, pur_tes$Estimated.Half.life > 0 & pur_tes$Estimated.Half.life < 1)
pur_tes_sub <- subset(pur_tes, pur_tes$Estimated.Half.life > 0 & grepl(filter.rep, pur_tes$Top.Hit))
fin_tes <- read.csv("./genome_fin/TEfamilies_CLwise_final_table.csv", as.is = T)
#fin_tes_sub <- subset(fin_tes, fin_tes$Estimated.Half.life > 0 & fin_tes$Estimated.Half.life < 1)
fin_tes_sub <- subset(fin_tes, fin_tes$Estimated.Half.life > 0 & grepl(filter.rep, fin_tes$Top.Hit))


#_______________________________________________
## ANOVA

req.repeat.class <- c("LTR")


## now only taking the LTR elements into consideration
all_tes_row <- rbind(cbind(HF = pur_tes_sub$Estimated.Half.life, GEN = "PUR"), cbind(HF = fin_tes_sub$Estimated.Half.life, GEN = "FIN"), cbind(HF = uso_tes_sub$Estimated.Half.life, GEN = "USO"), cbind(HF = hum_tes_sub$Estimated.Half.life, GEN = "HUM"), cbind(HF = mor_tes_sub$Estimated.Half.life, GEN = "MOR"))

tes.aov <- aov(all_tes_row[,1] ~ all_tes_row[,2])
summary(tes.aov)
#tes.aov$assign

par(mar=c(4.5, 6, 3, 3))
plot(TukeyHSD(tes.aov), las = 1)

### 
meansd <- function(x){return(paste(round(median(x),digits = 4), "; s.d. = ", round(sd(x),digits = 3), sep =""))}

meansd(pur_tes_sub$Estimated.Half.life)
meansd(uso_tes_sub$Estimated.Half.life)
meansd(fin_tes_sub$Estimated.Half.life)
meansd(hum_tes_sub$Estimated.Half.life)
meansd(mor_tes_sub$Estimated.Half.life)


## Difference in the means between the genomes
## We do the Wilcox-Mann Whitney test

getwilcoxmann <- function(x,y){t = wilcox.test(x,y,conf.int=T); return(paste(round(as.numeric(t$estimate),digits = 3), ", p = ", round(t$p.value,digits = 4), sep = ""))}

#t.test(pur_tes_sub$Estimated.Half.life, fin_tes_sub$Estimated.Half.life)
getwilcoxmann(pur_tes_sub$Estimated.Half.life, fin_tes_sub$Estimated.Half.life)

#t.test(pur_tes_sub$Estimated.Half.life, uso_tes_sub$Estimated.Half.life)
getwilcoxmann(pur_tes_sub$Estimated.Half.life, uso_tes_sub$Estimated.Half.life)

#t.test(pur_tes_sub$Estimated.Half.life, mor_tes_sub$Estimated.Half.life)
getwilcoxmann(pur_tes_sub$Estimated.Half.life, mor_tes_sub$Estimated.Half.life)

#t.test(pur_tes_sub$Estimated.Half.life, hum_tes_sub$Estimated.Half.life)
getwilcoxmann(pur_tes_sub$Estimated.Half.life, hum_tes_sub$Estimated.Half.life)

#t.test(mor_tes_sub$Estimated.Half.life, hum_tes_sub$Estimated.Half.life)
getwilcoxmann(mor_tes_sub$Estimated.Half.life, hum_tes_sub$Estimated.Half.life)

#t.test(mor_tes_sub$Estimated.Half.life, fin_tes_sub$Estimated.Half.life)
getwilcoxmann(mor_tes_sub$Estimated.Half.life, fin_tes_sub$Estimated.Half.life)

#t.test(mor_tes_sub$Estimated.Half.life, uso_tes_sub$Estimated.Half.life)
getwilcoxmann(mor_tes_sub$Estimated.Half.life, uso_tes_sub$Estimated.Half.life)


## Adding significance to the box plot

addSegment <- function(x1, x2, y, label, lwd, cex){
  segments(x0 = x1, x1 = x2, y0 = y, y1 = y, lwd = lwd)
  segments(x0 = x1, x1 = x1, y0 = y -1, y1 = y, lwd = lwd)
  segments(x0 = x2, x1 = x2, y0 = y -1, y1 = y, lwd = lwd)
  text(x = (x1 + x2)/2, y = y + 5, labels = label, cex = cex)
}

cex.plot = 1.2
line.width = 1.2
cex.low = 0.8
cex.high = 2.5
start.bars <- 58 

pdf("~/Downloads/Figure4.pdf")
figure3 <- boxplot(100 * pur_tes_sub$Estimated.Half.life, 100 * fin_tes_sub$Estimated.Half.life, 100 * uso_tes_sub$Estimated.Half.life, 100 * hum_tes_sub$Estimated.Half.life, 100 * mor_tes_sub$Estimated.Half.life, names = c("PK", "FIN", "USO", "HUM", "MOR"), ylim = c(0, 80), col = brewer.pal(5, "Pastel1"), notch = F, outline = T, cex.axis = cex.plot, xlab = "Estimated half life for LTR elements", cex.lab = cex.plot, ylab = "% divergence")
addSegment(1, 4, start.bars+5, "*", lwd = line.width, cex = cex.high)
addSegment(1, 5, start.bars+15, "**", lwd = line.width, cex = cex.high)
dev.off()
#### Figure 3



hist(as.numeric(subset(all_tes_row, all_tes_row[,2] == "PUR")[,1]), col = brewer.pal(5, "Pastel1")[1], xlim = c(0,1), breaks = 10, ylim = c(0,20))
hist(as.numeric(all_tes_row[,1][which(all_tes_row[,2] == "USO")]), col = brewer.pal(5, "Pastel1")[2], xlim = c(0,1), breaks = 10, add = T)

req.repeat.table <- mor_tes_sub
req.repeat.class <- unique(req.repeat.table$Top.Hit[grep("LTR", req.repeat.table$Top.Hit)])
color_rep <- brewer.pal(length(req.repeat.class), "Set2")

hist(-100, plot=T, xlim = c(0, max(req.repeat.table$Estimated.Half.life)), ylim = c(0,15), xlab = "Estimated half-life (% divergence)", main = "Histogram of TEs in PUR")
ind = 1
for (r in req.repeat.class){
  hist(add = T, req.repeat.table$Estimated.Half.life[which(req.repeat.table$Top.Hit == r)], col = color_rep[ind])
  ind = ind + 1
}

### Try out a phylogenetic tree for the LTR elements in all the genomes

# join all the LTR elements in the genomes

library("seqinr")
filter.rep <- "LTR"

genome_dir <- "~/Documents/ExploringTEsinCannabisGenome/genome_mor/"
genome_id <- "MOR"
## Inputting data
contigsfasta <- read.fasta(file = file.path(genome_dir, "Repeat_contigs_RE_filtered_min500bp.fa"), seqtype = "DNA",as.string = TRUE, set.attributes = FALSE, forceDNAtolower= F)
genome_tes <- read.csv(file.path(genome_dir, "TEfamilies_CLwise_final_table.csv"), as.is = T)
genome_tes_cls <- subset(genome_tes, genome_tes$Estimated.Half.life > 0 & grepl(filter.rep, genome_tes$Top.Hit))$Cluster.ID

filterIDs <- names(contigsfasta)[which(matrix(unlist(strsplit(names(contigsfasta), split = "Contig")), ncol = 2, byrow = T)[,1] %in% genome_tes_cls)]
reqFasta <- contigsfasta[which(names(contigsfasta) %in% filterIDs)]

write.fasta(reqFasta, names = paste(genome_id, names(reqFasta), sep = "_"), file.out = file.path(genome_dir, "Repeat_contigs_RE_LTR_elements.fa"), as.string = T)




