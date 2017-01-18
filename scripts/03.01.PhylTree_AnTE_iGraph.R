# Script to generate trees from AnTE

args  <- commandArgs(TRUE)

workDir <- args[1]
#   paste(args[1], ".filter.txt", sep="")
setwd(workDir)

library("igraph")
it<-0 #iteration
tree<-read.csv(paste("tree_obs_0_0_",it,sep=""),sep=" ",row.names=1) #change filenames as neccessary
candidate_tree<-read.csv(paste("candidate_tree_0_0_",it,sep=""),sep=" ",row.names=1)
cand_netw<-candidate_tree[1000,] #use the last step of MCMC (or others as desired)
obs_netw<-tree[1000,]
G<-graph.empty(length(cand_netw))
V(G)$size=14
for(i in 1:length(cand_netw))
{
  if(cand_netw[i]>-1)
  {
    G<-add.edges(G, c(from=i,to=cand_netw[i]+1))
  }
}
E(G)$curved<-0.25
E(G)$width<-4
E(G)$arrow.mode<-1
GY<-G
for(i in 1:length(obs_netw))
{
  new_v_index<-length(V(GY))+1
  GY<-add.vertices(GY,1)
  V(GY)[new_v_index]$name=""
  V(GY)[new_v_index]$size=3
  GY<-add.edges(GY,c(from=new_v_index,to=obs_netw[i]+1))
}

GX<-delete.vertices(GY,which(degree(GY)==0))
E(GX)$curved[is.na(E(GX)$curved)]<-0
E(GX)$width[is.na(E(GX)$width)]<-.5
E(GX)$arrow.mode[is.na(E(GX)$arrow.mode)]<-0

epsOut  <- paste(args[2],".eps",sep="")
# jpeg(filename = "phytree_igraph")
setEPS()
postscript(epsOut)
plot.igraph(GX,vertex.size=V(GX)$size,vertex.label=V(GX)$name)
dev.off()

