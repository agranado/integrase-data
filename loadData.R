

#useful commands for terminal: 

#generate lineage based on the manual distance measure 
# manualTree = upgma(as.dist(t(manualDist(as.character(barcodes),0.3,2/3,2 ))));manualTree$tip.label<- paste(names(barcodes),barcodes,sep="_")


rm(list=ls())
library("ape")
library("phangorn")
library(factoextra)

source("simulation2.R")
source("simMemoirStrDist3.R")

#load functions
#source("../simulation/simulation2.R")
#source("../simulation/simMemoirStrDist.R")

#choose the tree from the list:
file.idx=8


files = c("pos31_2","pos37","pos1_cont","pos4_cont","pos5_cont","pos6_cont","pos8_cont","pos33")
ks= c(8,5,4,3,3,3,5,6)
estimMu.s=c(0.3,0.5,0.4,0.4,0.4,0.4,0.4,0.5)
alphas=c(2/3,2/3,2/3,2/3,2/3,2/3,2/3,2/3)
generations=c(2,2,4,4,3,3,4,2)


file.path="/Users/alejandrog/MEGA/Caltech/trees/GraceData/"
fileName=files[file.idx]
fileName=paste(file.path,fileName,sep="")
k=ks[file.idx]
estimMu = estimMu.s[file.idx]
estimAlpha=alphas[file.idx]
estimG =generations[file.idx]

#save plot for continous recording trees 
plot.path=paste(file.path,"plots/contRecording",sep="")


#read csv as dataframe
 posInfo=read.csv(paste(fileName,".csv",sep=""))

 #features are in
 #posInfo$Sumary   includes the barcodes for each cell
 #posInfo$Movie.ID includes the actual ID of the cell, to compare with the nwk tree
 #posInfo$FISH.ID includes the ID from the FISH, not necessary
 
 #Sometimes there is no ID bc the cells did not appear in the field of view or bc they have to be excluded for diverse reasons. 
 #When there is no ID, the value becomes NA
 
 
 #barcodes =posInfo$Summary
 barcodes = posInfo$Summary[!is.na(posInfo$Movie.ID) & !posInfo$Summary=="xxxxxx"]
 names(barcodes) =posInfo$Movie.ID[!is.na(posInfo$Movie.ID)  & !posInfo$Summary=="xxxxxx"]
 
 #nuclear.cfp = posInfo$Nuclear.CFP[!is.na(posInfo$Movie.ID) & !posInfo$Summary=="xxxxxx"] 
 #names(nuclear.cfp) = posInfo$Movie.ID[!is.na(posInfo$Movie.ID)  & !posInfo$Summary=="xxxxxx"]
 
 #names(barcodes)=posInfo$Movie.ID
#convert data to fasta format
 #this function inserts a big R at the end of the sequence so it has at least one "c"
 #which will be constant. We do this because some calculations of base frequency need all actg bases to be present.
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

# #manual tree using the hclust methods from R heatmap
# h=heatmap.2(matdist_+t(matdist_),trace="none",dendrogram = 'column')
# 
# #alternative w/o plotting the actual heatmap, only hclust method
# hclust.tree=as.phylo(hclust(as.dist(t(matdist_))))
# hclust.tree$tip.label = treeUPGMA$tip.label
# manualTree=hclust.tree
# 




manualTree$tip.label<- paste(names(barcodes),barcodes,sep="_")
#sometimes there are edges with negative values (-6e-18) which should not happen
manualTree$edge.length[manualTree$edge.length<0]=0
#hc.manual=as.hclust(reverseLabels(manualTree))
hc.manual=as.hclust(manualTree)
# for saving ggsave("membow_31_2tree.pdf", device=CairoPDF)
fviz_dend(hc.manual, k = k, cex = 0.8, horiz = TRUE,  k_colors = "jco", 
          rect = TRUE, rect_border = "jco", rect_fill = TRUE,xlab="time",ylab="cells")

#plotting the real tree
pos21tree=read.tree(file=paste(fileName,".nwk",sep=""))


#split the name of the tips //Grace named them as 180_xx so we split using _
#we need to take the barcodes fromt the posInfo
true.tips=strsplit(pos21tree$tip.label,"_")
true.tips.mat=do.call(rbind,true.tips)

upgma.tips=strsplit(treeUPGMA$tip.label,"_")
upgma.tips.mat=do.call(rbind,upgma.tips)


nCells= dim(upgma.tips.mat)[1]; 
#cells that appear in the reconstructed tree
#we need to look for them in the true tips, edit the name and at the end, collapse and rename the tips in the tree
#only those cells that appear in the upgma tree will have a barcode
for (bc in 1:nCells){
  thisCell =upgma.tips.mat[bc,1]
  #this line basically matches the barcodes in the upgma tree to the movie IDs in the real tree. 
  #needs to be simplified
  true.tips.mat[which(true.tips.mat[,2]==thisCell),2]=paste(true.tips.mat[which(true.tips.mat[,2]==thisCell),2],upgma.tips.mat[which(upgma.tips.mat[,1]==thisCell),2])
}
#update the tree with the barcodees for those cells that we are considering
pos21tree$tip.label =  apply(true.tips.mat,1,paste,collapse="_")

#plot barcode statistics 
#including all barcodes in the original data (includes xxxxx, and cells with no Movie.ID)
x11()
barplot(prop.table(table(posInfo$Summary)),las=2,ylab="Freq",main="Barcode dist (all cells)")






