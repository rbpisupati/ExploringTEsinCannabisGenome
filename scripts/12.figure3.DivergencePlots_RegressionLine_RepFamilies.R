# Plot for divergence
library(RColorBrewer)
setwd("~/Documents/Projects/Summer'14/")
setwd("/lustre/scratch/users/rahul.pisupati/exploringRE/genome_uso/")
cex.plot <- 1.3
colors.div <- brewer.pal(4, "Set2")

subcov5k.pi.file <- "/lustre/scratch/users/rahul.pisupati/exploringRE/genome_uso/SRR351494.subcov100.mod.pi"
subcov5k.pi.file <- "~/mygit/ExploringTEsinCannabisGenome/genome_uso/windPivalues_meanDepth.txt"
subcov5k.pi.file <- "./GitHub/ExploringTEsinCannabisGenome/genome_pk/SRR352164.subcov100.mod.pi"
subcov5k.pi.file <- "/lustre/scratch/users/rahul.pisupati/exploringRE/genome_pk/SRR352164.subcov100.mod.pi"

pi_values <- read.csv(file=subcov5k.pi.file,header=FALSE,sep="\t")
log_cove <- log10(x = pi_values$V7)

ylim <- c(1,3)

color.addalpha <- function(reqCol, alpha = 50){
  alpha = (100-alpha)*255/100
  newCol = rgb(col2rgb(reqCol)[1], col2rgb(reqCol)[2], col2rgb(reqCol)[3], alpha = alpha, max = 255) 
  return(newCol)
}

plot(pi_values$V6,log_cove,xlab = "Divergence",ylab = "Logarithmic   Coverage", pch=19, cex.lab=cex.plot, cex.axis=cex.plot, col = color.addalpha(colors.div[2], 50))

#_________________________
#Plot the graphs in a layout 
#
par(fig=c(0,0.8,0,0.8), new=TRUE)
plot(pi_values$V6,log_cove,xlab = "Divergence",ylab = "Logarithmic   Coverage", pch=19, cex.lab=cex.plot, cex.axis=cex.plot, col = colors.div[1],ylim = ylim)
par(fig=c(0,0.8,0.55,1), new=TRUE)
boxplot(pi_values$V6, horizontal = T, axes = F)
par(fig=c(0.65,1,0,0.8),new=TRUE)
boxplot(log_cove, axes = F)
mtext("Enhanced Scatterplot", side=3, outer=TRUE, line=-3)

## Check for the specific sequences
xlim = 0.01
ylim = 2.2
abline(v=xlim)
abline(h=ylim)
pi_values[which((pi_values$V6 < xlim) & (log_cove > ylim)), c("V1","V2")]

check.seq <- c(608,1)
check.ind <- which((pi_values$V1 == paste("CL", check.seq[1], sep="")) & (pi_values$V2 == paste("Contig", check.seq[2], sep="")))
points(pi_values$V6[check.ind], log_cove[check.ind], col = "red", pch = 19)

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

layout(matrix(c(1,1,2,2,0,3,3,0), 2, 4, byrow = T))
cex.plot = 1.3
cex.point = 0.8
cex.high.point=1.8
line.width=1.6

plot(pi_values$V6,log_cove,xlab = "Divergence",ylab = "Logarithmic   Coverage", pch=19, cex=cex.point, main = "(a)", col = color_backgroud, cex.lab=cex.plot, cex.axis=cex.plot, cex.main = 1.5*cex.plot)
#_________________________
rep_maximaCopia <- subset(pi_values, pi_values$V1 =="CL41" | pi_values$V1 =="CL77" | pi_values$V1 =="CL100")
log_maximCopia <- log10(x = rep_maximaCopia$V7)
points(rep_maximaCopia$V6, log_maximCopia, col=color_rep[5], pch=19, cex=cex.high.point)
sepCluster <- rep_maximaCopia
clusters <- unique(sepCluster$V1)
for (i in 1:length(clusters)){
  req <- subset(sepCluster,  sepCluster$V1 == clusters[i] )
  log_req <- log10(x = req$V7)
  if (nrow(req) > 1){
    abline(lm(log_req~req$V6), col = color_rep[5], lwd=line.width)
  }
}
#_________________________
plot(pi_values$V6,log_cove,xlab = "Divergence",ylab = "Logarithmic   Coverage", pch=19, cex=cex.point, main = "(b)", col = color_backgroud, cex.lab=cex.plot, cex.axis=cex.plot, cex.main = 1.5*cex.plot)
#_________________________

rep_chromoGypsy <- subset(pi_values, pi_values$V1 =="CL19" | pi_values$V1 =="CL20" | pi_values$V1 =="CL28" | pi_values$V1 =="CL40" | pi_values$V1 =="CL67" | pi_values$V1 =="CL83" | pi_values$V1 =="CL84" | pi_values$V1 =="CL125" | pi_values$V1 =="CL126" | pi_values$V1 =="CL129" | pi_values$V1 =="CL161")
log_chromoGypsy <- log10(x = rep_chromoGypsy$V7)
points(rep_chromoGypsy$V6, log_chromoGypsy, col=color_rep[3], pch=19, cex=cex.high.point)
clusters <- unique(rep_chromoGypsy$V1)
for (i in 1:length(clusters)){
  req <- subset(rep_chromoGypsy,  rep_chromoGypsy$V1 == clusters[i] )
  log_req <- log10(x = req$V7)
  if (nrow(req) > 1){
    abline(lm(log_req~req$V6), col = color_rep[3], lwd=line.width)
  }
}
#_________________________
rep_athilaGypsy <- subset(pi_values, pi_values$V1 =="CL29" | pi_values$V1 =="CL109" | pi_values$V1 =="CL138")
log_athila <- log10(x = rep_athilaGypsy$V7)
points(rep_athilaGypsy$V6, log_athila, col=color_rep[4], pch=19, cex=cex.high.point)
sepCluster <- rep_athilaGypsy
clusters <- unique(sepCluster$V1)
for (i in 1:length(clusters)){
  req <- subset(sepCluster,  sepCluster$V1 == clusters[i] )
  log_req <- log10(x = req$V7)
  if (nrow(req) > 1){
    abline(lm(log_req~req$V6), col = color_rep[4], lwd=line.width)
  }
}
#_________________________
plot(pi_values$V6,log_cove,xlab = "Divergence",ylab = "Logarithmic   Coverage", pch=19, cex=cex.point, main = "(c)", col = color_backgroud,cex.lab=cex.plot, cex.axis=cex.plot, cex.main = 1.5*cex.plot)
#_________________________
rRNA <- subset(pi_values,  pi_values$V1 =="CL52" | pi_values$V1 =="CL56" | pi_values$V1 =="CL62" | pi_values$V1 =="CL76" | pi_values$V1 =="CL80" | pi_values$V1 =="CL82")
log_rrna <- log10(x = rRNA$V7)
points(rRNA$V6, log_rrna, col=color_rep[2], pch=19, cex=cex.high.point)
#sepCluster <- rRNA
#clusters <- unique(sepCluster$V1)
# for (i in 1:length(clusters)){
#   req <- subset(sepCluster,  sepCluster$V1 == clusters[i] )
#   log_req <- log10(x = req$V7)
#   if (nrow(req) > 1){
#     abline(lm(log_req~req$V6), col = color_rep[2], lwd=line.width)
#   }
# }
#_________________________
rep_simple <- subset(pi_values, pi_values$V1 == "CL6" | pi_values$V1 =="CL46" | pi_values$V1 =="CL57" | pi_values$V1 =="CL133" | pi_values$V1 =="CL211" | pi_values$V1 =="CL216" | pi_values$V1 =="CL218" | pi_values$V1 =="CL232" | pi_values$V1 =="CL241")
log_simple <- log10(x = rep_simple$V7)
points(rep_simple$V6, log_simple, col=color_rep[1], pch=19, cex=cex.high.point)
#sepCluster <- rep_simple
#clusters <- unique(sepCluster$V1)
#for (i in 1:length(clusters)){
#  req <- subset(sepCluster,  sepCluster$V1 == clusters[i] )
#  log_req <- log10(x = req$V7)
#  if (nrow(req) > 1){
#    abline(lm(log_req~req$V6), col = color_rep[1], lwd = line.width)
#  }
#}
#_________________________


clusters <- unique(rep_maximaCopia$V1)
colors <- rainbow(length(clusters))
for (i in 1:length(clusters)){
  req <- subset(rep_maximaCopia,  rep_maximaCopia$V1 == clusters[i] )
  if (nrow(req) > 1){
    abline(lm(req$V7~req$V6), col = colors[i])
  }
}

summary(pi_values$V6)
abline(v = mean(pi_values$V6)+1.5*sd(pi_values$V6))
length(which(pi_values$V6-mean(pi_values$V6) < 0))

clusters <- unique(pi_values$V1)
colors <- rainbow(length(clusters))
for (i in 1:length(clusters)){
  req <- subset(pi_values,  pi_values$V1 == clusters[i] )
  if (nrow(req) > 1){
    abline(lm(req$V7~req$V6), col = colors[i])
  }
}

cor.test(pi_values$V6,pi_values$V7)
hist(pi_values$V6)
length(clusters)