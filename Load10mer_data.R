# 29 Aug
# Load data for reconstruction of 10mer data for integrase paper


#useful commands for terminal:

#generate lineage based on the manual distance measure
# manualTree = upgma(as.dist(t(manualDist(as.character(barcodes),0.3,2/3,2 ))));manualTree$tip.label<- paste(names(barcodes),barcodes,sep="_")


rm(list=ls())
library("ape")
library("phangorn")
library(factoextra)

#erase chace files
system("rm /Users/alejandrog/MEGA/Caltech/trees/GraceData/10mer/editRate/*")

source("simulation2.R")
source("simMemoirStrDist3.R")

#choose the tree from the list:

fasta=F
clust.method=1
plot.all = 0

file.path="/Users/alejandrog/MEGA/Caltech/trees/GraceData/10mer/"

#save plot for continous recording trees
plot.path=paste(file.path,"plotsReconstruction/",sep="")

plots.edit.rate = paste(file.path,"plotsEditRate/",sep="")


existing.files = list.files("../GraceData/10mer/")
#read the xls with all the files info (some are not in the folder)
all.files=read.csv(paste(file.path,"results.guide.csv",sep=""))
#read.newick(text=toString(all.files$newick[1]))

    #this file information
##for(file.idx in c(10)){
for(file.idx in 1:length(all.files$file.name)){
#for(file.idx in 21){
    files.extension = all.files$file.extension[file.idx]
    file.name=all.files$file.name[file.idx]
    stacked.rate.file=paste(file.path,"editRate/",toString(all.files$file.name[file.idx]),"_all_rates.txt",sep="",collapse="")
    fileName=paste(file.path,file.name,sep="")

    if(length(grep(file.name,existing.files))!=0){ #only if the corresponding txt file is in the directory
      ########
        #how many groups for the colors in the dendrogram
        ks= rep(4,length(all.files$Position))


        k=ks[file.idx]
        #estimMu = estimMu.s[file.idx]
        #estimAlpha=alphas[file.idx]
        estimG =all.files$generations[file.idx]


    #read csv as dataframe
        posInfo=read.csv(paste(fileName,files.extension,sep=""),sep="\t")

        #estimate the edit rate from the fraction of unedited sites
        lengthBarcode = 10
        avgEditRate = array()
        avgEditRate_r = array()
        avgEditRate_x = array()
        nCells_ =  length(posInfo$state)

       #features are in
       #Sometimes there is no ID bc the cells did not appear in the field of view or bc they have to be excluded for diverse reasons.
       #When there is no ID, the value becomes NA

       barcodeLength = nchar(as.character(posInfo$state[1]))
       #barcodes =posInfo$Summary
       #FILTERING:::::
       # £ £ £  # # # # # # 
       barcodes = posInfo$state[!is.na(posInfo$cell) & !posInfo$state==paste(rep("x",barcodeLength),collapse="")]
       names(barcodes) =posInfo$cell[!is.na(posInfo$cell)  & !posInfo$state==paste(rep("x",barcodeLength),collapse="")]

      #After filtering the NA and XXXX.. we can calculte the edit rate per cell
      nCells_=length(as.character(barcodes))
      avgEditRate = array(); avgEditRate_r=array(); avgEditRate_x=array()
      for(i in 1:nCells_){
         avgEditRate[i] = sum(strsplit(toString(barcodes[i]),"")[[1]]=="u")/lengthBarcode
         avgEditRate_r[i] = sum(strsplit(toString(barcodes[i]),"")[[1]]=="r")/lengthBarcode
         avgEditRate_x[i] = sum(strsplit(toString(barcodes[i]),"")[[1]]=="x")/lengthBarcode
       }
       pos.vals=!(avgEditRate_x==0 & avgEditRate_r==0)
       estimMu = 1-mean(avgEditRate[pos.vals])^(1/estimG)
       estimAlpha = 1-mean(avgEditRate_x[pos.vals] / (avgEditRate_r[pos.vals] + avgEditRate_x[pos.vals]))

       # # # paper plot
       #Lets calculate the edit rate per site (in theory should be even)
       tt=apply(as.matrix(barcodes),2,strsplit,"")
       mat.sites = t(do.call(cbind,do.call(cbind,tt)))
       rate.per.site = apply(mat.sites!="u",2,sum)/dim(mat.sites)[1]


       rate.per.site = t(rate.per.site)
       #save to file for further processing and storage
       #FILE 1
       rate.file=paste(file.path,"editRate/All_rates.txt",sep="",collapse="")
       write.table(rate.per.site,rate.file, append=T, row.names =F,col.names=F)
       #order: x, r, u
       #save rates as appended file independently: WORKS
       alphabet =c("u","r","x")
       #FILE 2
       for (letter in alphabet){
         this.file=paste(file.path,"editRate/",letter,"_rates.txt",sep="",collapse="")
         letter.per.site = t(apply(mat.sites==letter,2,sum)/dim(mat.sites)[1])
         write.table(letter.per.site,this.file, append=T, row.names =F,col.names=F)
       }
       #FILE 3
       #individual edit rate (average for each tree)
       write.table(mean(rate.per.site),paste(file.path,"editRate/distribution_rates_perTree.txt",sep="",collapse="") ,
        append=T,row.names = F,col.names = F)

      #FILE 4
      #Write all barcodes for all cells for all trees together in a single file
      all.bc.file = paste(file.path,"editRate/","allBarcodes.txt",sep="")
      write.table(mat.sites,all.bc.file,append=T,row.names=F,col.names =F)

       #names(barcodes)=posInfo$Movie.ID
       #convert data to fasta format
       #this function inserts a big R at the end of the sequence so it has at least one "c"
       #which will be constant. We do this because some calculations of base frequency need all actg bases to be present.
        if(fasta==T){
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


        }


        if(clust.method==1 && plot.all == 1){
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
          manualTree$tip.label<- paste(names(barcodes),barcodes,sep="_")
          #sometimes there are edges with negative values (-6e-18) which should not happen
          manualTree$edge.length[manualTree$edge.length<0]=0
          #hc.manual=as.hclust(reverseLabels(manualTree))
          hc.manual=as.hclust(manualTree)
          # for saving ggsave("membow_31_2tree.pdf", device=CairoPDF)
      #    x11()
          fviz_dend(hc.manual, k = k, cex = 1.2, horiz = TRUE,  k_colors = "jco",
                    rect = TRUE, rect_border = "jco", rect_fill = TRUE,xlab="time",ylab="cells")
          pdf.path = paste(plot.path,all.files$file.name[file.idx],".pdf",sep="")
          scaling=0.7
          ggsave(pdf.path, device=cairo_pdf,width = 9.32*scaling,height = 10.4*scaling,units="in")

          ground.truth = all.files$newick[file.idx]
          true.tree=read.newick(text=toString(ground.truth))
          pdf.path.truth = paste(plot.path,all.files$file.name[file.idx],"_groundTruth.pdf",sep="")
          pdf(pdf.path.truth)
            plot.phylo(true.tree)
          dev.off()

          #join the pdfs
          system(paste("\'/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py\' -o ",plot.path,file.name,"_merged.pdf ",
                pdf.path," ",pdf.path.truth, sep=""))

          #remove individual files:
          system(paste("rm ",pdf.path,sep=""))
          system(paste("rm ",pdf.path.truth,sep=""))
        }
    # # # # # # # # # # #
     # # # # # # # # # # #
    #Alternative: use the hclust method from R

      if(clust.method ==2){
        matdist_names = matdist_
        colnames(matdist_names)<- paste(names(barcodes),barcodes,sep="_")
        dend<-hclust(as.dist(t(matdist_names)))
        fviz_dend(dend, k = k, cex = 0.8, horiz = TRUE,  k_colors = "jco",
                  rect = TRUE, rect_border = "jco", rect_fill = TRUE,xlab="time",ylab="cells")
      }

    }



}

#FIGURES
#plot edit rates per site as a barplot with error bars WORKS

barplot.edit.rate<-function(rate.file="/Users/alejandrog/MEGA/Caltech/trees/GraceData/10mer/editRate/All_rates.txt"){

  all.rates = as.matrix(read.table(rate.file))
  data<-data.frame(
    name= as.character(1:10), value=apply(all.rates,2,mean), sd=apply(all.rates,2,sd),order=1:10
  )
  x11()
  ggplot(data) + geom_bar( aes(x=reorder(name,order),y=value), stat="identity",fill="skyblue",alpha=0.7 )  +
    geom_errorbar( aes(x=name, ymin = value-sd,ymax = value+sd ),width=0.4,colour="orange",alpha=0.9,size=1.3 ) +
    theme(text = element_text(size=20)) + labs(x = "Element in array",y="Freq edited sites")
            #axis.text.x = element_text(angle=90, hjust=1))

}


#paper plot
#plot the frequency of each edit per site as a stacked barplot
#this plots uses ggplot2
stackedplot.edit.rate<-function(rate.file="/Users/alejandrog/MEGA/Caltech/trees/GraceData/10mer/editRate/All_rates.txt"){

  alphabet = c("u","r","x")
  N.sites=10
  all.rates = matrix(0,length(alphabet), N.sites)
  for(i in 1:length(alphabet)){
    letter = alphabet[i]
    this.file=paste(file.path,"editRate/",letter,"_rates.txt",sep="",collapse="")
    letter.rates = as.matrix(read.table(this.file))
    all.rates[i,] = apply(letter.rates,2,mean)
  }

  #since it is the average, not all colums add up to 1:
  norm.rates = t(t(all.rates) / apply(all.rates,2,sum))
  colnames(norm.rates)<-as.character(1:10)
  row.names(norm.rates)<-alphabet

  coul = brewer.pal(length(alphabet),"Pastel2")
  all.rates=norm.rates
  #make a data frame. Note that we can specify the order of factor to be plotted (order of stacks in the plot)
  data<-data.frame(
    letter = factor(c(rep("u",N.sites),rep("r",N.sites),rep("x",N.sites)), levels=c("r","x","u")),
    siten=factor(names(c(all.rates[1,],all.rates[2,],all.rates[3,])), levels=as.character(1:10)  ),
    value=c(all.rates[1,],all.rates[2,],all.rates[3,])
  )

  #making a stacked barplot with ggplot2
  ggplot() +
    geom_bar(aes(y=value,x=siten,fill=letter),data=data,stat="identity") +
    scale_fill_brewer(palette="Pastel2")

  #NOTE
  #This would be the old-school way with native R plotting:
  #barplot(norm.rates,col=coul,ylab = "freq",xlab="site number",
          #names.arg = as.character(1:10))

  #legend("topright", inset=c(0.8,0.8),
          #      fill=coul,
          #      legend=alphabet,bg="white")


}

#plot 3
histogram.plot.editrate<-function(file.path=file.path){

  rate.dist.file = paste(file.path,"editRate/distribution_rates_perTree.txt",sep="",collapse="")
  rate.dist = read.table(rate.dist.file)
}

#plot 4
#read ALL files in the results/ folder and calculate a global rate of edits, including "w" & "R"
stacked.plot.allTrees<-function(which.sites=1:10){
    results.folder = "/Users/alejandrog/Downloads/results/"
    all.results = list.files(results.folder)

    alphabet = c("u","r","x")
    N.sites=10

    all.mat.sites = matrix(, nrow=0,ncol=10)
    rate.per.site = c()
    for(i in 1:length(all.results)){
        posInfo=read.table(paste(results.folder,all.results[i],sep=""),header=T)
        barcodes = posInfo$state[!is.na(posInfo$cell) & !posInfo$state==paste(rep("x",barcodeLength),collapse="")]
        names(barcodes) =posInfo$cell[!is.na(posInfo$cell)  & !posInfo$state==paste(rep("x",barcodeLength),collapse="")]
        tt=apply(as.matrix(barcodes),2,strsplit,"")
        mat.sites = t(do.call(cbind,do.call(cbind,tt)))
        #save the average edit rate for each tree:
        #choose which sites to include in the analysis (to look for biases in edit rates)
        rate.per.site[i] = mean(apply(mat.sites[,which.sites]!="u",2,sum)/dim(mat.sites)[1])

        all.mat.sites= rbind(all.mat.sites,mat.sites)
    }

    freq.per.site=matrix(0,length(alphabet),dim(all.mat.sites)[2])
    j=1
    for (letter in alphabet){
      freq.per.site[j,] = t(apply(all.mat.sites==letter,2,sum)/dim(all.mat.sites)[1])
      j=j+1
    }

    coul = brewer.pal(length(alphabet),"Pastel2")
    #all.rates=norm.rates
    #make a data frame. Note that we can specify the order of factor to be plotted (order of stacks in the plot)
    all.rates  = freq.per.site
    all.rates = t(t(all.rates) / apply(all.rates,2,sum))
    colnames(all.rates)<-as.character(1:10)
    row.names(all.rates)<-alphabet
    data<-data.frame(
      letter = factor(c(rep("u",N.sites),rep("r",N.sites),rep("x",N.sites)), levels=c("r","x","u")),
      siten=factor(names(c(all.rates[1,],all.rates[2,],all.rates[3,])), levels=as.character(1:10)  ),
      value=c(all.rates[1,],all.rates[2,],all.rates[3,])
    )

    #making a stacked barplot with ggplot2

    ggplot() +
      geom_bar(aes(y=value,x=siten,fill=letter),data=data,stat="identity") +
      scale_fill_brewer(palette="Pastel2")

    return(list(rate.per.site,data))
}
#use the object returned above as input for these two functions:
plot.stacked.frequencies<-function(a){

  geom_bar(aes(y=value,x=siten,fill=letter),data=a[[2]],stat="identity") +
+     scale_fill_brewer(palette="Pastel2") + labs(x="Site N",y="Frequency",title="Frequency of edits per site (48hrs)")
}
#histogram
plot.histogram.edit.rates<-function(a){

  qplot(a[[1]],geom="histogram",main=paste("Edit rate distribution, n=",toString(length(a[[1]])),sep=""),
  binwidth=0.12,xlab="Fraction of edited sites",ylab="N of colonies",
  fill=I("grey"),col=I("black"))

}

# # alternative plotting method for dendrogram
# library(dendextend)
# library(ggdendro)
# aa<-hclust(as.dist(t(matdist_)))
# dend<-as.dendrogram(aa) %>%
#   set("branches_k_color",k=6) %>% set("branches_lwd",1.2) %>%
#   set("labels_colors") %>% set("labels_cex",c(0.9)) %>% set("leaves_pch",15) %>%
#   set("leaves_col",c("black"))
# plot(dend)
#
#
# #plotting the real tree
# pos21tree=read.tree(file=paste(fileName,".nwk",sep=""))
#
#
#
#
#
# #split the name of the tips //Grace named them as 180_xx so we split using _
# #we need to take the barcodes fromt the posInfo
# true.tips=strsplit(pos21tree$tip.label,"_")
# true.tips.mat=do.call(rbind,true.tips)
#
# upgma.tips=strsplit(treeUPGMA$tip.label,"_")
# upgma.tips.mat=do.call(rbind,upgma.tips)
#
#
# nCells= dim(upgma.tips.mat)[1];
# #cells that appear in the reconstructed tree
# #we need to look for them in the true tips, edit the name and at the end, collapse and rename the tips in the tree
# #only those cells that appear in the upgma tree will have a barcode
# for (bc in 1:nCells){
#   thisCell =upgma.tips.mat[bc,1]
#   #this line basically matches the barcodes in the upgma tree to the movie IDs in the real tree.
#   #needs to be simplified
#   true.tips.mat[which(true.tips.mat[,2]==thisCell),2]=paste(true.tips.mat[which(true.tips.mat[,2]==thisCell),2],upgma.tips.mat[which(upgma.tips.mat[,1]==thisCell),2])
# }
# #update the tree with the barcodees for those cells that we are considering
# pos21tree$tip.label =  apply(true.tips.mat,1,paste,collapse="_")
#
# #plot barcode statistics
# #including all barcodes in the original data (includes xxxxx, and cells with no Movie.ID)
# x11()
# barplot(prop.table(table(posInfo$Summary)),las=2,ylab="Freq",main="Barcode dist (all cells)")
