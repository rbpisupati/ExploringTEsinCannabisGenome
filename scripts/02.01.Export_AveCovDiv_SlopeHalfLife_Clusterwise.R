# Print average coverage, divergence, slope and halfLife for the cluster with many contigs
library(data.table)
args  <- commandArgs(TRUE)

#subcov5k.pi.file <- "/lustre/scratch/users/rahul.pisupati/exploringRE/genome_finola/SRR351929.subcov100.mod.pi"
#re.summary.file <- "~/mygit/ExploringTEsinCannabisGenome/genome_fin/summary_RE.csv" 

subcov5k.pi.file <- args[1]
re.summary.file <- args[2]

out.file <- "TEfamilies_CLwise_final_table.csv"
re.summary <- fread(re.summary.file, header = F)
pi_values <- read.csv(file=subcov5k.pi.file,header=FALSE,sep="\t")

clusters <- as.character(unique(pi_values$V1))

final.tefam <- as.data.frame(matrix(nrow=length(clusters),ncol=11))
final.tefam.cols <- c("Cluster ID",	"Average Coverage",	"SD Coverage",	"Average Divergence",	"SD Divergence",	"Slope",	"Estimated Half-life",	"Top Hit",	"No. of Hits",	"%cluster",	"%genome")

for (cl in clusters){
  cl.data <- pi_values[which(pi_values$V1 == cl),]
  cl.ind <- which(clusters == cl)
  final.tefam[cl.ind,1] <- cl
  if(nrow(cl.data) > 1){
    # Average copy number
    final.tefam[cl.ind,2] <- mean(cl.data$V7)
    final.tefam[cl.ind,3] <- sd(cl.data$V7)
    # Average divergence
    final.tefam[cl.ind,4] <- mean(cl.data$V6)
    final.tefam[cl.ind,5] <- sd(cl.data$V6)
    # Printing slope of regression line
    slope <- lm(log10(x = cl.data$V7)~cl.data$V6)
    final.tefam[cl.ind,6] <- as.character(slope$coefficients[2])
    final.tefam[cl.ind,7] <- as.numeric(-log(2)/slope$coefficients[2])
    top.hit.row <- re.summary[which(re.summary$V2 == cl),]
    top.hit.family <- gsub("\\)", "",gsub("\\(","",sub(",","",unlist(strsplit(top.hit.row$V7, split = " ")))))
    try(final.tefam[cl.ind,c(8,9,10)] <- top.hit.family, silent = T)
    try(final.tefam[cl.ind,11] <- top.hit.row$V5, silent = T)
  }
}

colnames(final.tefam) <- final.tefam.cols

write.csv(x = final.tefam, file = out.file)

