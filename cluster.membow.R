#clustering analysis of membow trees

clust.data = read.table(paste(file.path,"clusters.membow.data.txt",sep=""),header=T,sep="\t")

less.than100 = clust.data$colony<100
more.than100 = clust.data$colony>=100
clust.data$colony[less.than100] = paste("results_0",clust.data$colony[less.than100] ,sep=""  )
clust.data$colony[more.than100] = paste("results_",clust.data$colony[more.than100] ,sep=""  )

#parse clusters
who.is.here = list()
clusters = clust.data$clusters[3]
clusters = strsplit(toString(clusters),"_")[[1]]
for (c in 1:length(clusters)){

  #for this cluster we know which cells were supposed to be classified here
  who.is.here[[c]] = as.numeric(strsplit(clusters[c],",")[[1]])

}


#plot histogram of fraction correct
ggplot(clust.data,aes(x =1-missclassified.cells/cells)) +
geom_histogram(bins=5,color="darkblue",fill="lightblue") +
labs(x = "Fraction of correct cells",y="Number of colonies") +
theme(text = element_text(size=21))
