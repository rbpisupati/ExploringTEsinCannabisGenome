
## Generate a csv file from the blast out single file JSON from NCBI blast
library(jsonlite)

jsonFile <- "~/Downloads/blast_uso31.json"

blastjson <- fromJSON(jsonFile)
blast_title <- character()
blast_hit <- character()
blast_id <- character()
blast_num <- numeric()
for (i in 1:length(blastjson$BlastOutput2$report$results$search$query_title)){
  if (!is.null(blastjson$BlastOutput2$report$results$search$hits[[i]]$num[1])){
    title = blastjson$BlastOutput2$report$results$search$query_title[i]
    hitname = blastjson$BlastOutput2$report$results$search$hits[[i]]$description[[1]]$title[1]
    hitnum = sum(blastjson$BlastOutput2$report$results$search$hits[[i]]$hsps[[1]]$score)
    hitid = blastjson$BlastOutput2$report$results$search$hits[[i]]$description[[1]]$id[1]
    blast_title <- c(blast_title, title)
    blast_hit <- c(blast_hit, hitname)
    blast_id <- c(blast_id, hitid)
    blast_num <- c(blast_num, hitnum)
  }
}
df <- data.frame(query_title=blast_title, hit_id=blast_id, hit_num=blast_num,hit_name=blast_hit)
write.csv(df, file=paste(strsplit(jsonFile, "\\.")[[1]][1], ".df.csv", sep = ""))

#### Filter the RE contigs from the blast hits csv, added a column filter with 1's which need to get filtered
#### https://docs.google.com/spreadsheets/d/1uQOHwS5Wnwgpu4_I0Z3zWkmAYxrf947ULPpWLAq9HUU/edit#gid=1793392466
##

library("seqinr")
genome_dir <- "/vol/HOME/mygit/ExploringTEsinCannabisGenome/genome_uso/"

## Inputting data
contigsfasta <- read.fasta(file = file.path(genome_dir, "Repeat_contigs_RE_min500bp.fa"), seqtype = "DNA",as.string = TRUE, set.attributes = FALSE, forceDNAtolower= F)
blast_hits <- read.csv(file.path(genome_dir, "blast_hits.csv"))

## Filtering
filterIDs <- as.character(blast_hits$query_title[which(blast_hits$filter == 1)])
reqFasta <- contigsfasta[which(!names(contigsfasta) %in% filterIDs)]

write.fasta(reqFasta, names = names(reqFasta), file.out = file.path(genome_dir, "Repeat_contigs_RE_filtered_min500bp.fa"), as.string = T)
