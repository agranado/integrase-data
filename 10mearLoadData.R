setwd("/Users/alejandrog/MEGA/Caltech/trees/GIT")

source("simulation2.R")
source("simMemoirStrDist3.R")
file.dir = "/Users/alejandrog/MEGA/Caltech/trees/GraceData/10mer/"

fileName = paste(file.dir,"manual10mer",sep="")



nCells = 19
nUnits = 10 
#For this particular tree :
estimG = 4; estimAlpha = 2/3; estimMu = 0.4; 

posInfo=read.csv(paste(fileName,".csv",sep=""))
posInfo=posInfo[1:nCells,1:nUnits+1]
barcodes=array()
for(i in 1:dim(posInfo)[1]){
  barcodes[i]=paste(posInfo[i,],collapse="")
}
names(barcodes)<-as.character(c(1:9,19,10:18))


alphabet = c("x","u","r")
char.num = c("0","1","2")

for(j in 1:length(alphabet)){
  barcodes=gsub(char.num[j],alphabet[j],barcodes)
}





fastaBarcodes = convertSimToFasta(barcodes)


fasIN <-paste(fileName,".fas",sep="")


rand.barcode.order<-sample(1:length(fastaBarcodes))
#We re-order the barcodes using a fixed (but random) order
fastaBarcodes<-fastaBarcodes[rand.barcode.order]
barcodes<-barcodes[rand.barcode.order]
write(fastaBarcodes,file=fasIN)




convertMemoirToDNA(fasIN)
#sequences have a "c" at the end
#now we can use the phyDat

memoirfas<-read.phyDat(fasIN,format="fasta",type="DNA")

#calculate distance methods

dm<-dist.ml(memoirfas)
dm.ham=dist.hamming(memoirfas)
print("distance calculated")
treeUPGMA<-upgma(dm)
treeUPGMA.ham <-upgma(dm.ham)

hc=as.hclust(reverseLabels(treeUPGMA))
x11() 
fviz_dend(hc, k = k, # Cut in four groups
          cex = 0.9, # label size
          k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE, # Add rectangle around groups
          rect_border = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"), 
          rect_fill = TRUE,
          main="reconstructed lineage (membow)",
          xlab="cells",ylab="time")

#horizontal tree
x11()
fviz_dend(hc, k = k, cex = 0.8, horiz = TRUE,  k_colors = "jco", 
          rect = TRUE, rect_border = "jco", rect_fill = TRUE,xlab="time",ylab="cells")



#manualTree = upgma(as.dist(t(manualDist(as.character(barcodes),estimMu,estimAlpha,estimG ))));
#TEST manualTree using the ML function (should work slightly better)
matdist_=manualDistML(as.character(barcodes),estimMu,estimAlpha,estimG )
manualTree = upgma(as.dist(t(matdist_)));

 #manual tree using the hclust methods from R heatmap
 h=heatmap.2(matdist_+t(matdist_),trace="none",dendrogram = 'column')
 #alternative w/o plotting the actual heatmap, only hclust method
 hclust.tree=as.phylo(hclust(as.dist(t(matdist_))))
 hclust.tree$tip.label = treeUPGMA$tip.label
 manualTree=hclust.tree
 




manualTree$tip.label<- paste(names(barcodes),barcodes,sep="_")
#sometimes there are edges with negative values (-6e-18) which should not happen
manualTree$edge.length[manualTree$edge.length<0]=0
#hc.manual=as.hclust(reverseLabels(manualTree))
hc.manual=as.hclust(manualTree)
# for saving ggsave("membow_31_2tree.pdf", device=CairoPDF)
x11()
fviz_dend(hc.manual, k = k, cex = 0.8, horiz = TRUE,  k_colors = "jco", 
          rect = TRUE, rect_border = "jco", rect_fill = TRUE,xlab="time",ylab="cells")

#plotting the real tree
pos21tree=read.tree(file=paste(fileName,".nwk",sep=""))





