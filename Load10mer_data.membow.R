# 29 Aug
# Load data for reconstruction of 10mer data for integrase paper


#useful commands for terminal:

#generate lineage based on the manual distance measure
# manualTree = upgma(as.dist(t(manualDist(as.character(barcodes),0.3,2/3,2 ))));manualTree$tip.label<- paste(names(barcodes),barcodes,sep="_")


rm(list=ls())
library("ape")
library("phangorn")
library(factoextra)
library("phytools")

#erase chace files
system("rm /Users/alejandrog/MEGA/Caltech/trees/GraceData/10mer_membow/editRate/*")

source("simulation2.R")
source("simMemoirStrDist3.R")

#choose the tree from the list:

fasta=F
clust.method=1
plot.all = 1
#directory for memoir trees
#file.path="/Users/alejandrog/MEGA/Caltech/trees/GraceData/10mer/"
file.path="/Users/alejandrog/MEGA/Caltech/trees/GraceData/10mer_membow/"
#read the xls with all the files info (some are not in the folder)
all.files=read.csv(paste(file.path,"results.guide.membow.csv",sep=""))

#use this line to copy files :
#for(i in 1:length(all.files$file.name)){ system(paste("cp /Users/alejandrog/Downloads/results/",all.files$file.name[i], ".txt ", file.path ,sep=""))    }

#all membow trees were copied to the folder already


#save plot for continous recording trees
plot.path=paste(file.path,"plotsReconstruction/",sep="")

plots.edit.rate = paste(file.path,"plotsEditRate/",sep="")


existing.files = list.files("../GraceData/10mer_membow/")

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
          #ADD a 0 if the number of the cell has only 1 digit
          #Substitute the u,r,x for 0,1,2
          tip.names = paste(names(barcodes),barcodes,sep="_")
          for(tt in 1:length(tip.names)){
            this.cell.n = str_extract(tip.names[tt],"\\d+")
            if(nchar(this.cell.n)==1){
              tip.names[tt]=paste("0",tip.names[tt],sep="")
            }

          }

          tip.names=gsub("r","2",tip.names)
          tip.names=gsub("x","1",tip.names)
          tip.names=gsub("u","0",tip.names)

          manualTree$tip.label<- tip.names
          #sometimes there are edges with negative values (-6e-18) which should not happen
          manualTree$edge.length[manualTree$edge.length<0]=0
          #hc.manual=as.hclust(reverseLabels(manualTree))
          hc.manual=as.hclust(manualTree)
          # for saving ggsave("membow_31_2tree.pdf", device=CairoPDF)
      #    x11()
      k=all.files$groups[file.idx]
          fviz_dend(hc.manual, k = k, cex = 1.2, horiz = TRUE,  k_colors = "jco",
                    rect = F, rect_border = "jco", rect_fill = TRUE,xlab="cells",ylab="time",
                  labels_track_height=2)
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

    } #end IF FILE EXISTS



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
stacked.plot.allTrees<-function(which.sites=1:10,which.folder="membow"){



    if(which.folder =="membow"){
        results.folder = "/Users/alejandrog/Downloads/results/"
        all.results = list.files(results.folder)
    #Correction for membow mode:
        membow.trees= as.character(85:158)
        membow.files = array()
        mm=1
        for(m in 1:length(membow.trees)){

          if(length(all.results[grep(membow.trees[m],all.results)])){
               membow.files[mm] = all.results[grep(membow.trees[m],all.results)]
               mm = mm +1
             }
        }
        all.results = membow.files ####

  }else if(which.folder=="memoir"){
    results.folder = "/Users/alejandrog/Downloads/results_48/"
    all.results = list.files(results.folder)  ####
  }


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

    return(list(rate.per.site,data,all.mat.sites))
}
#use the object returned above as input for these two functions:
plot.stacked.frequencies<-function(a){
  ggplot()+
  geom_bar(aes(y=value,x=siten,fill=letter),data=a[[2]],stat="identity") +
     scale_fill_brewer(palette="Pastel2") + labs(x="Site N",y="Frequency",title="Frequency of edits per site (48hrs)")
}
#histogram
plot.histogram.edit.rates<-function(a){

  qplot(a[[1]],geom="histogram",main=paste("Edit rate distribution, n=",toString(length(a[[1]])),sep=""),
  binwidth=0.12,xlab="Fraction of edited sites",ylab="N of colonies",
  fill=I("grey"),col=I("black"))

}


number.effective.states <-function(norm.rates){
#norm.states comes if you execute the inside of stackedplot.edit.rate

#>norm.rates[,1]
#      u         r         x
#0.4823883 0.3387156 0.1788961
  effective.states=array()
  for (i in 1:dim(norm.rates)[2]){
      pi = norm.rates[,i]
      effective.states[i] = 3 ^ - sum( pi * log(pi,base=3)  )

  }
  return(effective.states)
}

boxplot.ratio.compare <-function(){
  memoir=stacked.plot.allTrees(which.sites=1:10,which.folder = "memoir")
  membow=stacked.plot.allTrees(which.sites=1:10,which.folder = "membow")

  ratio.membow=membow[[2]]$value[membow[[2]]$letter=="r"]/membow[[2]]$value[membow[[2]]$letter=="x"]
  ratio.memoir=memoir[[2]]$value[memoir[[2]]$letter=="r"]/memoir[[2]]$value[memoir[[2]]$letter=="x"]

  d<-data.frame(
     time = factor(c(rep("24hrs" , 7) , rep("48hrs",7))  ,levels=c("24hrs","48hrs")),
     ratio = c(ratio.membow[c(1,3,4,5,7,8,9)],ratio.memoir[c(1,3,4,5,7,8,9)] )
     )



     ggplot(d,aes(time,ratio,fill=time)) + geom_boxplot(outlier.size=0) + scale_fill_manual(values=wes_palette(n=2, name="GrandBudapest2")) +
     geom_jitter(aes(time,ratio), position = position_jitter(width=0.1,height=0),
                 alpha=0.6,size=3,show_guide =F) +
     xlab("Induction time") +
     ylab("Ratio r/x") +

     theme(axis.title.x = element_text(size=20,hjust=0.5),
          axis.title.y = element_text(size=20,vjust=1),
          axis.text.x = element_text(size=14,color='black'),
          axis.text.y = element_text(size=14,color='black'),
          legend.text=element_text(size=16),
          legend.title=element_text(size=18),
          legend.key.size = unit(1.5,"cm")

      )



   }



#meeting October 3
#test the ratio of x /r for 24hrs and 48hrs
