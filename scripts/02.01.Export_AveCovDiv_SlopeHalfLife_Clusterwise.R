# Print average coverage, divergence, slope and halfLife for the cluster with many contigs

args  <- commandArgs(TRUE)
data <- read.csv(args[1], header=F, sep="\t")
logCove <- log (x = data$V7)

if(nrow(data) == 1){
  outLine <- paste(data$V7, "NA", data$V6, "NA", "NA", "NA", sep = "\t" )
  sink(file = "temp.fin")
  cat(outLine)
  sink()
}

else {
  # Average copy number
  aveCov <- mean(data$V7)
  sdCov <- sd(data$V7)
  
  # Average divergence
  aveDiv <- mean(data$V6)
  sdDiv <- sd(data$V6)
  
  # Printing slope of regression line
  slope <- lm(logCove~data$V6)
  halfLife <- -log(2)/slope$coefficients[2]
  
  # Variable for line
  outLine <- paste( aveCov, sdCov, aveDiv, sdCov, slope$coefficients[2], halfLife, sep = "\t")
  
  # Printing whole line into a temp file
  sink(file = "temp.fin")
  cat(outLine)
  sink()
}