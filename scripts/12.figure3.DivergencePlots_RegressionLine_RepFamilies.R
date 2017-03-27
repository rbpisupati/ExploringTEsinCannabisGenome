# Plot for divergence
library(RColorBrewer)
setwd("~/Documents/Projects/Summer'14/")
setwd("/lustre/scratch/users/rahul.pisupati/exploringRE/genome_uso/")
cex.plot <- 1.3
colors.div <- brewer.pal(4, "Set2")

subcov5k.pi.file <- "/lustre/scratch/users/rahul.pisupati/exploringRE/genome_uso/SRR351494.subcov100.mod.pi"
subcov5k.pi.file <- "~/mygit/ExploringTEsinCannabisGenome/genome_uso/windPivalues_meanDepth.txt"
subcov5k.pi.file <- "/lustre/scratch/users/rahul.pisupati/exploringRE/genome_pk/SRR352164.subcov100.mod.pi"

ylim <- c(1,3)

color.addalpha <- function(reqCol, alpha = 50){
  alpha = (100-alpha)*255/100
  newCol = rgb(col2rgb(reqCol)[1], col2rgb(reqCol)[2], col2rgb(reqCol)[3], alpha = alpha, max = 255) 
  return(newCol)
}

#_________________________
#Plot the graph for different families
#
library("data.table")
subcov5k.pi.file <- "~/Documents/ExploringTEsinCannabisGenome/genome_pk/SRR352164.subcov100.mod.pi"
pur_tes.file <- "~/Documents/ExploringTEsinCannabisGenome/genome_pk/TEfamilies_CLwise_final_table.csv"
pur_tes <- read.csv(pur_tes.file, as.is = T)
pi_values <- read.csv(file=subcov5k.pi.file,header=FALSE,sep="\t", as.is = T)
pi_values <- cbind(pi_values, V8 = log10(x = pi_values$V7))
colnames(pi_values) <- c("cl", "contig", "step", "V4", "V5", "divergence", "coverage", "logcov")

line.width <- 1.6
color_rep <- brewer.pal(5, "Set1")
color_rep <- as.character(sapply(color_rep, function(x){color.addalpha(x, alpha = 20)}))
cex.point <- 1.2

check.family <- c("Simple_repeat", "rRNA", "LTR/ERV1", "LTR/Gypsy")

plot(pi_values$divergence, pi_values$logcov,xlab = "Divergence",ylab = "Logarithmic Coverage", pch=19, cex.lab=cex.plot, cex.axis=cex.plot, type="n")
ind = 1
for (i in check.family){
  req.family.ind <- which(paste(pi_values$cl, pi_values$contig, sep = "") %in% pur_tes$query_sequence[which(pur_tes$repeat_class == i)])
  if (ind > 2){
    points(pi_values$divergence[req.family.ind],pi_values$logcov[req.family.ind], col=color_rep[ind], pch=19, cex=cex.point)
    abline(lm(pi_values$logcov[req.family.ind] ~ pi_values$divergence[req.family.ind]), col = color_rep[ind], lwd=line.width) 
  }
  ind = ind + 1
}






#_________________________
#Histograms for different LTR elements
divergenceCLwise <- read.csv(file="PUR.final.TEfam.csv", header = T)
repwiseCL  <- file("SRR352164_500bprepeats.sub5kcov.mod.pi", open = "r")
count=1
CLs <- c()
for (i in 1:nrow(divergenceCLwise)){
  if(divergenceCLwise$Top.Hit[i] == "Simple_repeat"){
    CLs[count] <- as.character(divergenceCLwise$Cluster.ID[i])
    count <- count + 1
  }
}
divergenceGypsy <- c()
count <- 1
repwiseCL  <- read.csv("SRR352164_500bprepeats.sub5kcov.mod.pi", sep = "\t", header = F)
for (i in 1:nrow(repwiseCL)){
  if(as.character(repwiseCL$V1[i]) %in% CLs){
    divergenceGypsy[count] <- repwiseCL$V6[i]
#    print(repwiseCL$V1[i])
    count <- count + 1
  }
}
boxplot(divergenceGypsy)


#----------------------------
# Color the divergence plot based on the TE family of top hit element
# Input: Divergence clusterwise TEfamiies.csv file
# 
setwd("~/Documents/Projects/Summer'14/genome_pk")
divergenceCLwise <- read.csv(file="PUR.final.TEfam.mod.csv", header = T)
repFamilies <- unique(divergenceCLwise$Top.Hit)
colors <- rainbow(length(repFamilies))
par(mfrow = c(4, 4))
repwiseCL  <- file("SRR352164_500bprepeats.sub5kcov.mod.pi", open = "r")
while (length(oneLine <- readLines(repwiseCL, n = 1)) > 0) {
  myLine <- unlist((strsplit(oneLine, ",")))
  clno <- strsplit(myLine, "\t")[[1]]
  colcl <- colors[match(divergenceCLwise$Top.Hit[match(clno[1],divergenceCLwise$Cluster.ID)], repFamilies)]
  logval <- log10(x = as.numeric(clno[7]))
  points(clno[6], logval, col=colcl, pch=19, cex=0.8)
} 
close(repwiseCL)

points(0.11, 5.4, col="#FF00FFFF", pch=19)

#_______________________________
  
  
cl41 <- subset(pi_values, pi_values$V1 == "CL41")
log_cl41 <- log10(x = cl41$V7)
points(cl41$V6, log_cl41, col="green", pch=15, cex=1.5)
abline(lm(log_cl41~cl41$V6), lwd=2, col="green")


### ______________
### Figure 3
### ______________

library(RColorBrewer)
#display.brewer.all(5)
color_rep <- brewer.pal(5, "Set1")
color_rep <- as.character(sapply(color_rep, function(x){color.addalpha(x, alpha = 20)}))
color_backgroud <- color.addalpha("darkgrey", 40)

getLMlines <- function(sepCluster, color.lines, line.width){
  clusters <- unique(sepCluster$cl)
  for (i in 1:length(clusters)){
    req <- subset(sepCluster,  sepCluster$cl == clusters[i] )
    if (nrow(req) > 1){
      abline(lm(req$logcov~req$divergence), col = color.lines, lwd=line.width)
    }
  }
}

pdf("~/Downloads/Figure3.pdf")
layout(matrix(c(1,1,2,2,0,3,3,0), 2, 4, byrow = T))
cex.plot = 1.3
cex.point = 0.8
cex.high.point=1.8
line.width=1.6


plot(pi_values$divergence,pi_values$logcov,xlab = "Divergence",ylab = "Logarithmic   Coverage", pch=19, cex=cex.point, main = "(a)", col = color_backgroud, cex.lab=cex.plot, cex.axis=cex.plot, cex.main = 1.5*cex.plot)
#_________________________
rep_maximaCopia <- subset(pi_values, pi_values$cl %in% c("CL41", "CL77", "CL100"))
points(rep_maximaCopia$divergence,rep_maximaCopia$logcov, col=color_rep[5], pch=19, cex=cex.high.point)
getLMlines(rep_maximaCopia, color_rep[5], line.width)


#_________________________
plot(pi_values$divergence,pi_values$logcov,xlab = "Divergence",ylab = "Logarithmic   Coverage", pch=19, cex=cex.point, main = "(b)", col = color_backgroud, cex.lab=cex.plot, cex.axis=cex.plot, cex.main = 1.5*cex.plot)

rep_chromoGypsy <- subset(pi_values, pi_values$cl %in% c("CL19", "CL20", "CL28", "CL40", "CL67", "CL83", "CL84", "CL125", "CL126", "CL129", "CL161"))
points(rep_chromoGypsy$divergence,rep_chromoGypsy$logcov, col=color_rep[3], pch=19, cex=cex.high.point)
getLMlines(rep_chromoGypsy, color_rep[3], line.width)

rep_athilaGypsy <- subset(pi_values, pi_values$cl %in% c("CL29", "CL109", "CL138"))
points(rep_athilaGypsy$divergence,rep_athilaGypsy$logcov, col=color_rep[4], pch=19, cex=cex.high.point)
getLMlines(rep_athilaGypsy, color_rep[4], line.width)

#_________________________
plot(pi_values$divergence,pi_values$logcov,xlab = "Divergence",ylab = "Logarithmic   Coverage", pch=19, cex=cex.point, main = "(c)", col = color_backgroud,cex.lab=cex.plot, cex.axis=cex.plot, cex.main = 1.5*cex.plot)

rRNA <- subset(pi_values,  pi_values$cl %in% c("CL52", "CL56", "CL62", "CL76", "CL80", "CL82"))
points(rRNA$divergence,rRNA$logcov, col=color_rep[2], pch=19, cex=cex.high.point)

rep_simple <- subset(pi_values, pi_values$cl %in% c("CL6", "CL46", "CL57", "CL133", "CL211", "CL216", "CL218", "CL232", "CL241"))
points(rep_simple$divergence,rep_simple$logcov, col=color_rep[1], pch=19, cex=cex.high.point)
#_________________________
dev.off()
