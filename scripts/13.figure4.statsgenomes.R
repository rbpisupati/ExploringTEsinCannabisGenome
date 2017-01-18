## ANOVA between the different genomes.
library(RColorBrewer)
setwd("~/Documents/Projects/Summer'14/")
## Environment files are from a new project
mor_tes <- read.csv("./TEs_times_genomes/MOR.final.TEfam.csv", as.is = T)
mor_tes_sub <- subset(mor_tes, mor_tes$Estimated.Half.life > 0 & mor_tes$Estimated.Half.life < 1)
#mor_tes_sub <- subset(mor_tes, mor_tes$Estimated.Half.life > 0)
hum_tes <- read.csv("./TEs_times_genomes/HUM.final.TEfam.csv", as.is = T)
hum_tes_sub <- subset(hum_tes, hum_tes$Estimated.Half.life > 0 & hum_tes$Estimated.Half.life < 1)
#hum_tes_sub <- subset(hum_tes, hum_tes$Estimated.Half.life > 0)
uso_tes <- read.csv("./TEs_times_genomes/USO.final.TEfam.csv", as.is = T)
uso_tes_sub <- subset(uso_tes, uso_tes$Estimated.Half.life > 0 & uso_tes$Estimated.Half.life < 1)
#uso_tes_sub <- subset(uso_tes, uso_tes$Estimated.Half.life > 0)
pur_tes <- read.csv("./TEs_times_genomes/PUR.final.TEfam.mod.csv", as.is = T)
pur_tes_sub <- subset(pur_tes, pur_tes$Estimated.Half.life > 0 & pur_tes$Estimated.Half.life < 1)
#pur_tes_sub <- subset(pur_tes, pur_tes$Estimated.Half.life > 0)
fin_tes <- read.csv("./TEs_times_genomes/FIN.final.TEfam.csv", as.is = T)
fin_tes_sub <- subset(fin_tes, fin_tes$Estimated.Half.life > 0 & fin_tes$Estimated.Half.life < 1)
#fin_tes_sub <- subset(fin_tes, fin_tes$Estimated.Half.life > 0)


#_______________________________________________
## ANOVA

all_tes_row <- rbind(cbind(HF = pur_tes_sub$Estimated.Half.life, GEN = "PUR"), cbind(HF = fin_tes_sub$Estimated.Half.life, GEN = "FIN"), cbind(HF = uso_tes_sub$Estimated.Half.life, GEN = "USO"), cbind(HF = hum_tes_sub$Estimated.Half.life, GEN = "HUM"), cbind(HF = mor_tes_sub$Estimated.Half.life, GEN = "MOR"))

tes.aov <- aov(all_tes_row[,1] ~ all_tes_row[,2])
summary(tes.aov)
#tes.aov$assign

plot(TukeyHSD(tes.aov), las = 1)


#Difference in the means between the genomes

pur_tes_sub$Estimated.Half.life
t.test(pur_tes_sub$Estimated.Half.life, fin_tes_sub$Estimated.Half.life)
t.test(pur_tes_sub$Estimated.Half.life, uso_tes_sub$Estimated.Half.life)
t.test(pur_tes_sub$Estimated.Half.life, mor_tes_sub$Estimated.Half.life)
t.test(pur_tes_sub$Estimated.Half.life, hum_tes_sub$Estimated.Half.life)

t.test(mor_tes_sub$Estimated.Half.life, hum_tes_sub$Estimated.Half.life)
t.test(mor_tes_sub$Estimated.Half.life, fin_tes_sub$Estimated.Half.life)
t.test(mor_tes_sub$Estimated.Half.life, uso_tes_sub$Estimated.Half.life)


## We do the Wilcox-Mann Whitney test
wilcox.test(pur_tes_sub$Estimated.Half.life, fin_tes_sub$Estimated.Half.life)
wilcox.test(pur_tes_sub$Estimated.Half.life, uso_tes_sub$Estimated.Half.life)
wilcox.test(pur_tes_sub$Estimated.Half.life, mor_tes_sub$Estimated.Half.life)
wilcox.test(pur_tes_sub$Estimated.Half.life, hum_tes_sub$Estimated.Half.life)

wilcox.test(mor_tes_sub$Estimated.Half.life, hum_tes_sub$Estimated.Half.life)
wilcox.test(mor_tes_sub$Estimated.Half.life, fin_tes_sub$Estimated.Half.life)
wilcox.test(mor_tes_sub$Estimated.Half.life, uso_tes_sub$Estimated.Half.life)


## Adding significance to the box plot

addSegment <- function(x1, x2, y, label, lwd, cex){
  segments(x0 = x1, x1 = x2, y0 = y, y1 = y, lwd = lwd)
  segments(x0 = x1, x1 = x1, y0 = y -1, y1 = y, lwd = lwd)
  segments(x0 = x2, x1 = x2, y0 = y -1, y1 = y, lwd = lwd)
  text(x = (x1 + x2)/2, y = y + 1.5, labels = label, cex = cex)
}

display.brewer.all(5)
figure3 <- boxplot(100 * pur_tes_sub$Estimated.Half.life, 100 * fin_tes_sub$Estimated.Half.life, 100 * uso_tes_sub$Estimated.Half.life, 100 * hum_tes_sub$Estimated.Half.life, 100 * mor_tes_sub$Estimated.Half.life, names = c("PK", "FIN", "USO", "HUM", "MOR"), ylim = c(0, 50), col = brewer.pal(5, "Pastel1"), ylab = "Estimated half life (% divergence)", notch = F, outline = T, cex = 0.7)
addSegment(1, 2, 30, "NS", lwd = 1, cex = 0.6)
addSegment(1, 3, 35, "NS", lwd = 1, cex = 0.6)
addSegment(1, 5, 45, "**", lwd = 1, cex = 1.7)
addSegment(1, 4, 40, "*", lwd = 1, cex = 1.5)

#####  ARCHIVE commands



all_tes_sub <- rbind(pur_tes_sub, mor_tes_sub, hum_tes_sub, fin_tes_sub, uso_tes_sub)

qqplot(100 * pur_tes_sub$Estimated.Half.life, 100 * hum_tes_sub$Estimated.Half.life, pch = 20)
abline(0,1)

all_tes_sub <- cbind(PUR = pur_tes_sub$Estimated.Half.life, FIN = fin_tes_sub$Estimated.Half.life, USO = uso_tes_sub$Estimated.Half.life, MOR = mor_tes_sub$Estimated.Halfw.life, HUM = hum_tes_sub$Estimated.Half.life)

all_tes_row <- rbind(cbind(HF = pur_tes_sub$Estimated.Half.life, GEN = "PUR"), cbind(HF = fin_tes_sub$Estimated.Half.life, GEN = "FIN"), cbind(HF = uso_tes_sub$Estimated.Half.life, GEN = "USO"), cbind(HF = hum_tes_sub$Estimated.Half.life, GEN = "HUM"), cbind(HF = mor_tes_sub$Estimated.Half.life, GEN = "MOR"))

all_tes_row <- rbind(cbind(HF = pur_tes$Estimated.Half.life, GEN = "PUR"), cbind(HF = fin_tes$Estimated.Half.life, GEN = "FIN"), cbind(HF = uso_tes$Estimated.Half.life, GEN = "USO"), cbind(HF = hum_tes$Estimated.Half.life, GEN = "HUM"), cbind(HF = mor_tes$Estimated.Half.life, GEN = "MOR"))

c(length(pur_tes_sub$Estimated.Half.life), length(mor_tes_sub$Estimated.Half.life), length(hum_tes_sub$Estimated.Half.life), length(fin_tes_sub$Estimated.Half.life), length(uso_tes_sub$Estimated.Half.life))


100 * c(mean(pur_tes_sub$Estimated.Half.life, na.rm = T), mean(fin_tes_sub$Estimated.Half.life, na.rm = T), mean(uso_tes_sub$Estimated.Half.life, na.rm = T), mean(mor_tes_sub$Estimated.Half.life, na.rm = T), mean(hum_tes_sub$Estimated.Half.life, na.rm = T))

100 * c(sd(pur_tes_sub$Estimated.Half.life, na.rm = T), sd(fin_tes_sub$Estimated.Half.life, na.rm = T), sd(uso_tes_sub$Estimated.Half.life, na.rm = T), sd(mor_tes_sub$Estimated.Half.life, na.rm = T), sd(hum_tes_sub$Estimated.Half.life, na.rm = T))

c("PUR","FIN","USO","MOR","HUM")



