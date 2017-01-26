library(jsonlite)

blastjson <- fromJSON("~/Downloads/blast_uso31.json")

num_hits = length(blastjson$BlastOutput2$report$results$search$query_title)
blast_title <- character()
blast_hit <- character()
blast_id <- character()
blast_num <- numeric()
blast_per <- numeric()
for (i in 1:num_hits){
  title = blastjson$BlastOutput2$report$results$search$query_title[i]
  hitname = blastjson$BlastOutput2$report$results$search$hits[[i]]$description[[1]]$title
  ## check  hitper = blastjson$BlastOutput2$report$results$search$hits[[i]]$hsps
  hitid = blastjson$BlastOutput2$report$results$search$hits[[i]]$description[[1]]$id
  hitnum = blastjson$BlastOutput2$report$results$search$hits[[i]]$num[1]
  blast_title <- c(blast_title, title)
  blast_hit <- c(blast_hit, hitname)
  blast_id <- c(blast_id, hitid)
  blast_num <- c(blast_num, hitnum)
}

df <- data.frame(query_title=blast_title, hit_id=blast_id,hit_name=blast_hit, hit_num=blast_num)

write.csv(df, file="blast_out.csv")

