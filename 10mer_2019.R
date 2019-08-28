# 29 Aug
# Load data for reconstruction of 10mer data for integrase paper


#useful commands for terminal:

#generate lineage based on the manual distance measure
# manualTree = upgma(as.dist(t(manualDist(as.character(barcodes),0.3,2/3,2 ))));manualTree$tip.label<- paste(names(barcodes),barcodes,sep="_")


rm(list=ls())
library("ape")
library("phangorn")
library(factoextra)
library(stringr)
library(phytools)
library(data.table)
library(cluster)
library(dplyr)
library(adephylo)
library(phylobase)



source("simulation2.R")
source("simMemoirStrDist3.R")

#ONLY when re-running the whole analysis
#erase chace files
# system("rm /Users/alejandrog/MEGA/Caltech/trees/GraceData/10mer_2019/editRate/*")

#choose the tree from the list:

#source("../GraceData/RF.experiment.R")
source("../integrase-data/RF.experiment.R")
fasta=F
clust.method=1
plot.all = 0

#vars 2019
translate = T # T if files come with 201 notation

# This indicates whether we are in mac or linux
# Folder where the integrase-specific scripts are stored
integrase_folder= "integrase-data/"

#file.path="/Users/alejandrog/MEGA/Caltech/trees/GraceData/10mer_2019/"
file.path=paste("/home/agranado/MEGA/Caltech/trees/",integrase_folder,"10mer_2019/",sep="")
#save plot for continous recording trees
plot.path=paste(file.path,"plotsReconstruction/",sep="")

plots.edit.rate = paste(file.path,"plotsEditRate/",sep="")


#existing.files = list.files("../GraceData/10mer_2019/FISH/")
existing.files = list.files(  paste("../",integrase_folder,"10mer_2019/FISH/",sep="") )
#read the xls with all the files info (some are not in the folder)
all.files=fread(paste(file.path,"results.guide2019.tsv",sep=""))
#read.newick(text=toString(all.files$newick[1]))

###########################
#Function WORKS REVISION
###########################

#works, goes through all trees in the folder and uses the selected clustering method
#returns a data frame
runAllTrees<-function(clust.method = "diana"){
    res.list=lapply(1:114,inspect.tree, plot.all = F,globalG=4,global = T,clust.method = clust.method)
    res.mat = do.call(rbind, lapply(res.list,unlist) )
    res.mat = as.data.frame(res.mat)
    names(res.mat)<-c("colony","idx","ncells","memoir","membow","genotypes","edits","cophe","entropy","dc")
    res.mat %>% arrange(desc(membow),desc(memoir)) -> res.mat

    return(res.mat)
}

#ttake the MEMOIR metric calculations and compare them across clustering methods
compareClustering<-function(methods = c("diana","complete","ward.D2","single","average"),pdf_file = ""){


  res.mat.list = list()
  colors = RColorBrewer::brewer.pal(12,"Set2") #max number of colors

  for(i in 1:length(methods)){
    #Function call
    #cat( paste (toString(i)," ", methods[i],"\n",sep="\t" ))
    res.mat.list[[i]] = runAllTrees(clust.method = methods[i])
  }

  if(pdf_file==""){
    x11();
  }else{
    pdf(pdf_file)
  }

  par(mfrow = c(1,2))
  #plot memoir distance
  for(i in 1:length(methods)){
    res.mat = res.mat.list[[i]]
    cc =  ecdf(res.mat$memoir)
    if(i==1){
      plot(seq(0,1,0.01),cc(seq(0,1,0.01)),type  ="l",col = colors[i],main = "Memoir",
      ylab = "Fraction of colonies",xlab = "Tree score",
      lwd =2,ylim = c(0,1))
    }else{
      lines(seq(0,1,0.01),cc(seq(0,1,0.01)),type  ="l",col = colors[i])
    }
  }

  #plot Membow distance
  for(i in 1:length(methods)){
    res.mat = res.mat.list[[i]]
    cc =  ecdf(res.mat$membow)
    if(i==1){
      plot(seq(0,1,0.01),cc(seq(0,1,0.01)),type  ="l",col = colors[i],main ="Membow",
      ylab = "Fraction of colonies",xlab = "Tree score",
      lwd =2,ylim = c(0,1))
    }else{
      lines(seq(0,1,0.01),cc(seq(0,1,0.01)),type  ="l",col = colors[i])
    }
  }

  if(pdf_file!=""){ dev.off()}

}

inspect.tree<-function(file.idx,plot.all = T,return.tree = F,global = T,globalG = 4,clust.method = "complete"){


  files.extension = all.files$file.extension[file.idx]
  file.name=all.files$file.name[file.idx]
  stacked.rate.file=paste(file.path,"editRate/",toString(all.files$file.name[file.idx]),"_all_rates.txt",sep="",collapse="")
  fileName=paste(file.path,file.name,sep="")

  if(length(grep(file.name,existing.files))!=0 & all.files$n.cells[file.idx]>2){ #only if the corresponding txt file is in the directory
    ########
      #how many groups for the colors in the dendrogram
      ks= rep(4,length(all.files$Position))
      k=ks[file.idx]
      #estimMu = estimMu.s[file.idx]
      #estimAlpha=alphas[file.idx]
      estimG =all.files$generations[file.idx]

  #read csv as dataframe
      posInfo = fread(paste(fileName,files.extension,sep=""),sep="\t",colClasses = 'character')
      names(posInfo)<-c("cell","state")
      if(translate){
          posInfo$state = str_replace_all(posInfo$state, c("2" = "r", "1" = "u","0"="x"))
          posInfo$state = str_replace_all(posInfo$state, c("3" = "u"))
        }

      posInfo$cell = substr(posInfo$cell,5,6)

      posInfo$cell = as.character(as.numeric(posInfo$cell))
      #new format ID for cells, take only the last two digits

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

    # 1.1 estimate parameters from the tree
    if(!global){
        param.list =estim.params(barcodes)
       #
        estimMu = param.list[[1]]
        estimAlpha = param.list[[2]]

    }else{
       #GLOBAL parameters , per site
       param.list = estim.params.global(estimG = globalG,file = paste("../",integrase_folder,"10mer_2019/editRate/allBarcodes.txt",sep="" ))
       estimMu = param.list[[1]] #these are vectors!
       estimAlpha = param.list[[2]]
    }

     #Lets calculate the edit rate per site (in theory should be even)
     tt=apply(as.matrix(barcodes),2,strsplit,"")
     mat.sites = t(do.call(cbind,do.call(cbind,tt)))
     rate.per.site = apply(mat.sites!="u",2,sum)/dim(mat.sites)[1]
     rate.per.site = t(rate.per.site)

     #order: x, r, u
     #save rates as appended file independently: WORKS
     alphabet =c("u","r","x")

        # 1.2- calculate distance matrix
        ground.truth = all.files$newick[file.idx]
        true.tree=read.newick(text=toString(ground.truth))

        if(global){
          matdist_=manualDistML_2(as.character(barcodes),estimMu,estimAlpha,estimG )
          #colnames(matdist_)<-barcodes
          #row.names(matdist_)<-barcodes
        }else{
          matdist_=manualDistML__(as.character(barcodes),estimMu,estimAlpha,estimG )
        }
        #translate tip names

                tip.names = paste(names(barcodes),barcodes,sep="_")
                for(tt in 1:length(tip.names)){
                  this.cell.n = str_extract(tip.names[tt],"\\d+")
                  if(nchar(this.cell.n)==1){
                    tip.names[tt]=paste("0",tip.names[tt],sep="")
                  }
                }

                tip.names=gsub("r","2",tip.names)
                tip.names=gsub("u","1",tip.names)
                tip.names=gsub("x","0",tip.names)


        if(clust.method =="diana"){
          #TEST FOR DIANA method # # # # # # 
          matdist_1 = matdist_
          colnames(matdist_1)<-tip.names;row.names(matdist_1)<-tip.names

          hclust.tree = as.phylo(as.hclust( diana(as.dist(t(matdist_1)))))
        }else{

        #alternative w/o plotting the actual heatmap, only hclust method
          hclust.tree=as.phylo(hclust(as.dist(t(matdist_)),method = clust.method))

        }

        manualTree = hclust.tree
        manualTree$tip.label<- tip.names
        #sometimes there are edges with negative values (-6e-18) which should not happen
        manualTree$edge.length[manualTree$edge.length<0]=0
        #hc.manual=as.hclust(reverseLabels(manualTree))
        hc.manual=as.hclust(manualTree)
        # for saving ggsave("membow_31_2tree.pdf", device=CairoPDF)

        # 2.- calculate the MEMOIR full distance RF
        # Also depurate the ground truth tree from dead cells

        ### EDIT FOR RF DISTRIBUTION 17th October 2018
        #this function calculates the RF distance between the reconstruction and the groun truth.
        #this should be comparable to previous measures of RF for the simulations.
        #RF =0 means random clustering
        #RF =1 means perfect reconstruction
        RF.list = normalized.RF.experiment(true.tree, manualTree , barcodes,posInfo)
        this.score = RF.list[[1]]
        alive.tree = RF.list[[2]]
        rand.dist = RF.list[[3]]

        #KF dist from the treedist package
        d6 =KF.dist(alive.tree,manualTree)
        # 3.- Alternative clustering method?
        #we can make a vector of distances :
        colnames(matdist_)<-tip.names;row.names(matdist_)<-tip.names
        d1=this.score
         diana_res = diana(matdist_+t(matdist_) )
         dc = diana_res$dc #divisive coefficient
        d2=RF.dist(as.phylo(as.hclust( diana_res ) ),alive.tree,normalize=T)
      #  if(rand.dist>0){d2 = 1-d2/rand.dist}
        d2 = 1-d2
        #distance 2


        # 4.- MEMbow distance:
        if(length(barcodes)>4){
          d3.list = normalized.RF.membow(true.tree,manualTree,barcodes,posInfo,make.plots=plot.all,
                  file.name = file.name, clust.method = clust.method, global = global,globalG=globalG,file.path=file.path)
          d3 =d3.list[[1]]
          n.colony=d3.list[[2]]

        }else { d3 = d1;alive.tree_ = alive.tree;n.colony = length(barcodes)}


        if(plot.all==T){
          #  x11()
            par(mfrow = c(2,2))
            par(family="mono")
            par(cex = 1.1)
            plot.phylo(alive.tree,main ="ground truth")
            plot.phylo(manualTree,main = paste("reconstruction:",toString(round(d1,digits = 2))))

            plot.phylo(d3.list[[3]],main = "reduced ground truth")
            plot.phylo(d3.list[[4]],main = paste("reduced reconstruction:", toString( round(d3,digits = 2)  )))
        }


        # 5.- Fraction of edited sites
        b=as.matrix(do.call(cbind,strsplit(substr(tip.names,4,14),"")))
        d4=sum(b=="1")/prod(dim(b)) #fraction of unedited sites (to filter bad colonies)
        #this.score is a colum vector of many clustering methods, based on the same distance matrix
        this.score = t( c (d1,d2,d3,d4))
        #just print the first score

        # 6.- Print results to terminal

        this.score = t(c(toString(all.files$file.name[file.idx]),d1,d3,d2,length(alive.tree$tip.label),d4))
      cat( paste("id:",file.idx,"name:", file.name, toString(file.idx),"Full:", toString(round(d1,digits =3)), "N=",toString(length(barcodes)) , "--- MEMbow:", toString(round(d3,digits = 3)),"N =",toString(n.colony),"Pr_edit",toString(1-d4),"\n",sep="\t" ))
        #WRITE FILE RF distance

        cophe = cor(cophenetic(as.hclust(manualTree)),as.dist(matdist_1 + t(matdist_1)))
        entro = barcodeEntropy(manualTree)
        #WRITE FILE NEWICK reconstructed

        if(!return.tree){
          return( list(file.name,file.idx,length(barcodes),d1,d3,n.colony, 1-d4,cophe,sum(entro),dc))
        }else{
          return(list(file.name,file.idx,length(barcodes),d1,d3,n.colony, 1-d4, alive.tree,manualTree,matdist_1))
        }

  }

}






#Clonal accuracy
# for a given tree it will return the genotypes that appeared in more than 2 cells along with their frequencies
# The returned object is a table (array with names)
get.all.clones<-function(file.idx){

  res<-inspect.tree(file.idx,plot.all = F,return.tree = T,global = F,globalG = 4,clust.method = "diana")

  if(length(res)>0){

      ground = res[[8]]
      reconstruction = res[[9]]

      genotypes = substr(ground$tip.label,4,14 )
      #table of frequencies for all genotypes
      a<-table(genotypes)
      #a[a>1]

    #      this_clone<-names(a[a>1])[7]
    #      this_size<-a[a>1][7]
    #      which_cells<-grep(this_clone , res[[8]]$tip.label)
    #Return the part of the table which includes the genotypes that appeared more than once
    #in this colony
      return(a[a>1])
    }else{
      return(NULL)
    }
}

# New method for clonal accuracy
# Aug 2019
# uses MRCA to calculate sub-trees and scores
clonal.score<-function(id,control=F){


      r = inspect.tree(id,plot.all = F, return.tree = T)
      # if ground truth is not good enough, then returns NULL so check for that
      if(length(r)>0){
        # tree is good quality then continue
        # inspect.tree returns a list where element 8 is the ground_truth tree
        ground_truth = r[[8]]

        # to calculate the random guess. negative CONTROL
        if(control)
          ground_truth$tip.label = sample(ground_truth$tip.label)




        dist_tips = distTips(ground_truth)
        dist_tips = as.matrix(dist_tips)
        barcodes = row.names(dist_tips)

        simple_barcodes= str_split(barcodes,"_",simplify=T)[,2]

        clones_here = get.all.clones(id)
        if(length(clones_here)>0){
          #for this barcodes, where are all the cells
          clone_score = c()
          #for this BARCode
          for(i in 1:length(clones_here)){
              #i = 1
              # how many cells
              clone_size = length( which(simple_barcodes==names(clones_here[i])) )
              masked_barcodes = simple_barcodes
              masked_barcodes[ which(simple_barcodes==names(clones_here[i])) ] <- rep("A",clone_size)
              masked_barcodes[masked_barcodes != "A"] = "X"
              # > masked_barcodes
              # [1] "X" "X" "X" "X" "X" "X" "X" "A" "A" "A" "X" "X" "X"
              #rename ground truth
              ground_truth_masked = ground_truth
              ground_truth_masked$tip.label = masked_barcodes
              #get most common ancestor of this clade (includes ALL A)
              clone_mrca = MRCA( ground_truth ,which(simple_barcodes==names(clones_here[i])))
              # this is the largest tree that includes all A's
              clone_tree = extract.clade(ground_truth_masked, clone_mrca)
              # Here we need a function that iterates over all mrca partitions
              # Get all subtrees
              tree_partitions = subtrees(clone_tree)
              # Save score for all partitions
              partition_score = c()
              # Iterate over all partitions:
              Tot_A = clone_size # How many A's in the whole MRCA
              for(t in 1:length(tree_partitions)){
                  clone_subtree = tree_partitions[[t]] #this is a list
                  # how many A's in this partition
                    S_a = length(grep("A",clone_subtree$tip.label))
                  partition_score[t] = S_a / ( Tot_A + length(clone_subtree$tip.label) - S_a)
                   # S_a  /  Tot_A   +  cells_considered - S_a
                   # S_a  / Tot_A + S_b
                   # length(grep("A",clone_tree$tip.label)) / ( clone_size + length(clone_tree$tip.label) - length(grep("A",clone_tree$tip.label)) )
              }

              clone_score[i] = max(partition_score)

              #masked_barcodes = paste(as.character(1:length(masked_barcodes)),masked_barcodes,sep="" )
              #clone_score[i] = length(grep("A",clone_tree$tip.label)) / ( clone_size + length(clone_tree$tip.label) - length(grep("A",clone_tree$tip.label)) )
            }
          return(clone_score)
        }else{
          return(NULL)
        }
      }else{
        return(NULL)
      }


}

#run for control and all data and make overlapping histogram
# Final plots (?)
make.clonal.plot<-function(){
  # Aug 11th 2019 works well
  #takes a minute
  a<-lapply(1:110,clonal.score)
  a_rand<-lapply(1:110,clonal.score,control=T)

  #make the plots
  hist(unlist(a),col=rgb(0,0,1,0.5)  ,xlim=c(0,1), ylim=c(0,250), main="Clonal accuracy for all clones", xlab="Score")
  hist(unlist(a_rand), col=rgb(1,0,0,0.5), add=T)
  legend("topright", c("All clones", "Random guess"), col=c("blue", "red"), lwd=2)

  #plot ECDF
  d<-ecdf(unlist(a))
  d_rand <-ecdf(unlist(a_rand))
  plot(seq(0,1,0.001) , d(seq(0,1,0.001)) ,type = "l",lwd = 2,
      main = "Clonal accuracy of all clones",xlab = "Score",
        ylab="Fraction of clones", col = "blue")
  lines(seq(0,1,0.001), d_rand(seq(0,1,0.001)),col ="red",lwd = 2)

}







old.forloop <-function(){
    #this file information
        #for(file.idx in c(3)){
    for(file.idx in 1:length(all.files$file.name)){
        #for(file.idx in 9){
        files.extension = all.files$file.extension[file.idx]
        file.name=all.files$file.name[file.idx]
        stacked.rate.file=paste(file.path,"editRate/",toString(all.files$file.name[file.idx]),"_all_rates.txt",sep="",collapse="")
        fileName=paste(file.path,file.name,sep="")

        if(length(grep(file.name,existing.files))!=0 & all.files$n.cells[file.idx]>2){ #only if the corresponding txt file is in the directory
          ########
            #how many groups for the colors in the dendrogram
            ks= rep(4,length(all.files$Position))


            k=ks[file.idx]
            #estimMu = estimMu.s[file.idx]
            #estimAlpha=alphas[file.idx]
            estimG =all.files$generations[file.idx]


        #read csv as dataframe
            posInfo = fread(paste(fileName,files.extension,sep=""),sep="\t",colClasses = 'character')
            names(posInfo)<-c("cell","state")
            if(translate) {
              posInfo$state = str_replace_all(posInfo$state, c("2" = "r", "1" = "u","0"="x"))
              posInfo$state = str_replace_all(posInfo$state, c("3" = "u"))
            }

            posInfo$cell = substr(posInfo$cell,5,6)

            posInfo$cell = as.character(as.numeric(posInfo$cell))
            #new format ID for cells, take only the last two digits


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
             avgEditRate[i] = sum(strsplit(toString(barcodes[i]),"")[[1]]=="u")/lengthBarcode #fraction of "u" in the array
             avgEditRate_r[i] = sum(strsplit(toString(barcodes[i]),"")[[1]]=="r")/lengthBarcode
             avgEditRate_x[i] = sum(strsplit(toString(barcodes[i]),"")[[1]]=="x")/lengthBarcode
           }
           pos.vals=!(avgEditRate_x==0 & avgEditRate_r==0) #which cells to take into accoutn

           estimMu = 1-mean(avgEditRate[pos.vals])^(1/estimG)
           estimAlpha = 1-mean(avgEditRate_x[pos.vals] / (avgEditRate_r[pos.vals] + avgEditRate_x[pos.vals]))
           #trying a small hack, 1-alpha
         #  estimAlpha = 1-estimAlpha
          #estimG = estimG * 10
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


          write.table(estimMu,paste(file.path,"editRate/distribution_Mu_perTree.txt",sep="",collapse="") ,
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

              #cmds[1]="sed -i.bak 's/R/c/g'"
              #cmds[2]="sed -i.bak 's/x/t/g'"
              #cmds[3]="sed -i.bak 's/u/g/g'"
              #cmds[4]="sed -i.bak 's/r/a/g'"

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


            if(clust.method==1){



              ground.truth = all.files$newick[file.idx]
              true.tree=read.newick(text=toString(ground.truth))


              #manualTree = upgma(as.dist(t(manualDist(as.character(barcodes),estimMu,estimAlpha,estimG ))));
              #TEST manualTree using the ML function (should work slightly better)
              matdist_=manualDistML__(as.character(barcodes),estimMu,estimAlpha,estimG )
              #colnames(matdist_)<-barcodes
              #row.names(matdist_)<-barcodes


              #alternative w/o plotting the actual heatmap, only hclust method
              hclust.tree=as.phylo(hclust(as.dist(t(matdist_))))

              manualTree = hclust.tree





              tip.names = paste(names(barcodes),barcodes,sep="_")
              for(tt in 1:length(tip.names)){
                this.cell.n = str_extract(tip.names[tt],"\\d+")
                if(nchar(this.cell.n)==1){
                  tip.names[tt]=paste("0",tip.names[tt],sep="")
                }
              }

              tip.names=gsub("r","2",tip.names)
              tip.names=gsub("u","1",tip.names)
              tip.names=gsub("x","0",tip.names)

              manualTree$tip.label<- tip.names
              #sometimes there are edges with negative values (-6e-18) which should not happen
              manualTree$edge.length[manualTree$edge.length<0]=0
              #hc.manual=as.hclust(reverseLabels(manualTree))
              hc.manual=as.hclust(manualTree)
              # for saving ggsave("membow_31_2tree.pdf", device=CairoPDF)
              #TEST FOR DIANA method # # # # # # 
              matdist_1 = matdist_
              colnames(matdist_1)<-tip.names;row.names(matdist_1)<-tip.names

              hc.manual_ = as.hclust( diana(as.dist(t(matdist_1))))

              ### EDIT FOR RF DISTRIBUTION 17th October 2018
              #this function calculates the RF distance between the reconstruction and the groun truth.
              #this should be comparable to previous measures of RF for the simulations.
              #RF =0 means random clustering
              #RF =1 means perfect reconstruction
              RF.list = normalized.RF.experiment(true.tree, manualTree , barcodes,posInfo)
              this.score = RF.list[[1]]
              alive.tree = RF.list[[2]]
              rand.dist = RF.list[[3]]

              #PLOTTING
              # # # # # # # #
               # # # # # # # #
              if(plot.all==1){
                  if(nCells_<4){
                    k=1
                  }

                  cex = 0.8
                  fviz_dend(hc.manual, k = k, cex = cex, horiz = TRUE,  k_colors = "jco",
                                      rect = F, rect_border = "jco", rect_fill = TRUE,xlab="cells",ylab="time",
                                    labels_track_height=2)

                  pdf.path = paste(plot.path,all.files$file.name[file.idx],".pdf",sep="")
                  scaling=0.7
                  ggsave(pdf.path, device=cairo_pdf,width = 9.32*scaling,height = 10.4*scaling,units="in")

                  pdf.path.truth = paste(plot.path,all.files$file.name[file.idx],"_groundTruth.pdf",sep="")
                  pdf(pdf.path.truth)
                    plot.phylo(true.tree,cex = cex)
                  dev.off()



                  #PDF WRITE for simplified ground.truth
                  pdf.path.truth.simple = paste(plot.path,all.files$file.name[file.idx],"_groundTruth_simple.pdf",sep="")
                  pdf(pdf.path.truth.simple)
                    plot.phylo(alive.tree,cex = cex)
                  dev.off()


                  #join the pdfs
                  system(paste("\'/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py\' -o ",plot.path,file.name,"_merged.pdf ",
                                  pdf.path," ",pdf.path.truth.simple," ",pdf.path.truth, sep=""))

                            #remove individual files:
                  system(paste("rm ",pdf.path,sep=""))
                  system(paste("rm ",pdf.path.truth,sep=""))
                  system(paste("rm ",pdf.path.truth.simple,sep=""))
              } #END PLOTTING
                # # # # # #
                 # # # # # #

              #we can make a vector of distances :
              colnames(matdist_)<-tip.names;row.names(matdist_)<-tip.names
              d1=this.score
              d2=RF.dist(as.phylo(as.hclust( diana(matdist_+t(matdist_) ) )),alive.tree,normalize=T)
            #  if(rand.dist>0){d2 = 1-d2/rand.dist}
              d2 = 1-d2
              #distance 2


              if(length(barcodes)>4){
                d3.list = normalized.RF.membow(true.tree,manualTree,barcodes,posInfo,global = F)
                d3 =d3.list[[1]]
                n.colony=d3.list[[2]]
              }else { d3 = d1;alive.tree_ = alive.tree;n.colony = length(barcodes)}

              # d3=RF.dist(as.phylo(upgma(matdist_1+t(matdist_1))), alive.tree  )
              # if(rand.dist>0){d3 = 1-d3/rand.dist}



              b=as.matrix(do.call(cbind,strsplit(substr(tip.names,4,14),"")))
              d4=sum(b=="1")/prod(dim(b)) #fraction of unedited sites (to filter bad colonies)
              #this.score is a colum vector of many clustering methods, based on the same distance matrix
              this.score = t( c (d1,d2,d3,d4))
              #just print the first score

              this.score = t(c(toString(all.files$file.name[file.idx]),d1,d3,d2,length(alive.tree$tip.label),d4))
            cat( paste(toString(file.idx),"Full:", toString(round(d1,digits =3)), "N=",toString(length(barcodes)) , "--- MEMbow:", toString(round(d3,digits = 3)),"N =",toString(n.colony),"\n",sep="\t" ))
              #WRITE FILE RF distance
                write.table(this.score,paste(file.path,"editRate/distribution_RF_perTree.txt",sep="",collapse="") ,
                          append=T,row.names = F,col.names = F)

              #WRITE FILE NEWICK reconstructed

              newick.out=paste(fileName,".nwk",sep="")
              write.tree(manualTree,file=newick.out)

            }



        #Alternative: use the hclust method from R

          if(clust.method ==2){
            matdist_names = matdist_
            colnames(matdist_names)<- paste(names(barcodes),barcodes,sep="_")
            dend<-hclust(as.dist(t(matdist_names)))
            fviz_dend(dend, k = k, cex = 0.8, horiz = TRUE,  k_colors = "jco",
                      rect = TRUE, rect_border = "jco", rect_fill = TRUE,xlab="time",ylab="cells")
          }

        }



    }# # # # # # # #
    # # # # # # # # #
    # # # # # # # # 
    #END MAIN FOR

}


#Clonal accuracy
# for a given tree it will return the genotypes that appeared in more than 2 cells along with their frequencies
# The returned object is a table (array with names)
get.all.clones<-function(file.idx){

  res<-inspect.tree(file.idx,plot.all = F,return.tree = T,global = F,globalG = 4,clust.method = "diana")

  if(length(res)>0){

      ground = res[[8]]
      reconstruction = res[[9]]

      genotypes = substr(ground$tip.label,4,14 )
      #table of frequencies for all genotypes
      a<-table(genotypes)
      #a[a>1]

    #      this_clone<-names(a[a>1])[7]
    #      this_size<-a[a>1][7]
    #      which_cells<-grep(this_clone , res[[8]]$tip.label)
    #Return the part of the table which includes the genotypes that appeared more than once
    #in this colony
      return(a[a>1])
    }else{
      return(NULL)
    }
}

# New method for clonal accuracy
# Aug 2019
# uses MRCA to calculate sub-trees and scores
clonal.score<-function(id,control=F){


      r = inspect.tree(id,plot.all = F, return.tree = T)
      # if ground truth is not good enough, then returns NULL so check for that
      if(length(r)>0){
        # tree is good quality then continue
        # inspect.tree returns a list where element 8 is the ground_truth tree
        ground_truth = r[[8]]

        # to calculate the random guess. negative CONTROL
        if(control)
          ground_truth$tip.label = sample(ground_truth$tip.label)




        dist_tips = distTips(ground_truth)
        dist_tips = as.matrix(dist_tips)
        barcodes = row.names(dist_tips)

        simple_barcodes= str_split(barcodes,"_",simplify=T)[,2]

        clones_here = get.all.clones(id)
        if(length(clones_here)>0){
          #for this barcodes, where are all the cells
          clone_score = c()
          best_Sb = c()
          best_Sa = c()
          clone_totaA = c()
          #for this BARCode
          for(i in 1:length(clones_here)){
              #i = 1
              # how many cells
              clone_size = length( which(simple_barcodes==names(clones_here[i])) )
              masked_barcodes = simple_barcodes
              masked_barcodes[ which(simple_barcodes==names(clones_here[i])) ] <- rep("A",clone_size)
              masked_barcodes[masked_barcodes != "A"] = "X"
              # > masked_barcodes
              # [1] "X" "X" "X" "X" "X" "X" "X" "A" "A" "A" "X" "X" "X"
              #rename ground truth
              ground_truth_masked = ground_truth
              ground_truth_masked$tip.label = masked_barcodes
              #get most common ancestor of this clade (includes ALL A)
              clone_mrca = MRCA( ground_truth ,which(simple_barcodes==names(clones_here[i])))
              # this is the largest tree that includes all A's
              clone_tree = extract.clade(ground_truth_masked, clone_mrca)
              # Here we need a function that iterates over all mrca partitions
              # Get all subtrees
              tree_partitions = subtrees(clone_tree)
              # Save score for all partitions
              partition_score = c()
              S_a_array = c()
              S_b_array = c()
              # Iterate over all partitions:
              Tot_A = clone_size # How many A's in the whole MRCA
              for(t in 1:length(tree_partitions)){
                  clone_subtree = tree_partitions[[t]] #this is a list
                  # how many A's in this partition
                    S_a = length(grep("A",clone_subtree$tip.label))
                  partition_score[t] = S_a / ( Tot_A + length(clone_subtree$tip.label) - S_a)
                   # S_a  /  Tot_A   +  cells_considered - S_a
                   # S_a  / Tot_A + S_b
                   # length(grep("A",clone_tree$tip.label)) / ( clone_size + length(clone_tree$tip.label) - length(grep("A",clone_tree$tip.label)) )
                   S_b = length(clone_subtree$tip.label) - S_a

                   # Save these numbers for false positive and false negative calculation
                   S_a_array[t] = S_a
                   S_b_array[t] = S_b


              }

              clone_score[i] = max(partition_score)
              # we can save the index of the best score and then get S_a, Tota_A, Sb
              best_score_index =which.max(partition_score)
              best_Sb[i] =  S_b_array[best_score_index]
              best_Sa[i] =  S_a_array[best_score_index]
              clone_totaA[i] = Tot_A
              #masked_barcodes = paste(as.character(1:length(masked_barcodes)),masked_barcodes,sep="" )
              #clone_score[i] = length(grep("A",clone_tree$tip.label)) / ( clone_size + length(clone_tree$tip.label) - length(grep("A",clone_tree$tip.label)) )
            }
          return(list(clone_totaA,clone_score, best_Sa, best_Sb))
        }else{
          return(NULL)
        }
      }else{
        return(NULL)
      }


}

#run for control and all data and make overlapping histogram
# Final plots (?)
make.clonal.plot<-function(){
  # Aug 11th 2019 works well
  #takes a minute
  a<-lapply(1:110,clonal.score)
  a_rand<-lapply(1:110,clonal.score,control=T)

  #make the plots
  hist(unlist(a),col=rgb(0,0,1,0.5)  ,xlim=c(0,1), ylim=c(0,250), main="Clonal accuracy for all clones", xlab="Score")
  hist(unlist(a_rand), col=rgb(1,0,0,0.5), add=T)
  legend("topright", c("All clones", "Random guess"), col=c("blue", "red"), lwd=2)

  #plot ECDF
  d<-ecdf(unlist(a))
  d_rand <-ecdf(unlist(a_rand))
  plot(seq(0,1,0.001) , d(seq(0,1,0.001)) ,type = "l",lwd = 2,
      main = "Clonal accuracy of all clones",xlab = "Score",
        ylab="Fraction of clones", col = "blue")
  lines(seq(0,1,0.001), d_rand(seq(0,1,0.001)),col ="red",lwd = 2)

}



 # # # # # #
# # # # # #
# OLD Methods for clonal accuracy

# clonal tracing accuracy  Jul 8th 2019
inspect.clusters<-function(file.idx,global = T,globalG = 4,clust.method = "complete",k= 4){
    res= inspect.tree(file.idx,plot.all=F,return.tree= T,global,globalG,clust.method)


    hc = res[[8]]
    hc_reconstructed = res[[9]]
    x11()
    p=fviz_dend(as.hclust(force.ultrametric(hc)), k = k, # Cut in four groups
              cex = 0.9, # label size
              k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
              color_labels_by_k = TRUE, # color labels by groups
              rect = TRUE, # Add rectangle around groups
              rect_border = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
              rect_fill = TRUE,
              main="reconstructed lineage (membow)",
              xlab="cells",ylab="time")

    p1=fviz_dend(as.hclust(force.ultrametric(hc_reconstructed)), k = k, # Cut in four groups
              cex = 0.9, # label size
              k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
              color_labels_by_k = TRUE, # color labels by groups
              rect = TRUE, # Add rectangle around groups
              rect_border = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
              rect_fill = TRUE,
              main="reconstructed lineage (membow)",
              xlab="cells",ylab="time")

    grid.arrange(p,p1,ncol = 2)
}
#   data.frame(RF=rate.dist,name=all.files$file.name[as.logical(all.files$FISH)])

clone_accuracy<-function(file.idx,k_ = 4,control =F){

    res<-inspect.tree(file.idx,plot.all = F,return.tree = T,global = F,globalG = 4,clust.method = "diana")



    ground = res[[8]]

    k = floor(length(ground$tip.label)/k_)
    #destroy lineage information as control
    #how much accuracy in clonal tracing you would get by random
    #because for small trees clonal tracing might be trivial
    if(control)
      ground$tip.label = sample(ground$tip.label)

    genotypes = substr(ground$tip.label,4,14)
    clones<-table(genotypes)
    clones<-clones[clones>1]

    if(length(clones)>0){
        #for each clone
        all_freqs = c()
        for(i in 1:length(names(clones))){
            cells_in_clone = names(clones)[i]
            clone_in_ground = grep( cells_in_clone ,ground$tip.label,value = T)

            ground_clusters = dendextend::cutree(force.ultrametric(ground),k =k)
            #get ground_cluster ID for these cells:
            clone_ID = ground_clusters[clone_in_ground]
            clone_freqs = table(clone_ID)/length(clone_in_ground)

            all_freqs[i] = clone_freqs
        }
        return(all_freqs)

    }else{
      return(NULL)
    }
}

compare_accuracy_control<-function(ks =c(4,3,2),large_colonies){

  #par(mfrow = c(ceiling(length(ks)/2),2))
  par(mfrow = c(2,2))
  for(k in ks){
    res_control<-lapply(large_colonies,clone_accuracy,k_ =k,control =T);
    res<-lapply(large_colonies,clone_accuracy,k_ =k,control =F);

    d_control = ecdf(unlist(res_control))
    d = ecdf(unlist(res))

    plot(seq(0,1,0.010),d(seq(0,1,0.010)),type = "l",lwd = 2,col ="blue",
    main ="Fraction of cells correctly classified (per clone)",
    xlab="Accuracy", ylab ="Fraction of clones")
    lines(seq(0,1,0.010),d_control(seq(0,1,0.010)),type = "l",lwd = 2,col ="red")

  }

}
# # # # # # # #
 # # # # # # #
#Additional plot schemes for paper FIGURES
# "npg" is the palette inspired by nature-like papers
plot.npg.tree<-function(){
  x11();fviz_dend(hc.manual, k = k, cex = 1.4, horiz = TRUE,  k_colors = "npg",
          rect = T, rect_border = "npg", rect_fill = T,xlab="cells",ylab="time",low_rect=100)
}

#FIGURES
#plot edit rates per site as a barplot with error bars WORKS

barplot.edit.rate<-function(rate.file="/Users/alejandrog/MEGA/Caltech/trees/GraceData/10mer_2019/editRate/All_rates.txt"){

  all.rates = as.matrix(read.table(rate.file))
  data<-data.frame(
    name= as.character(1:10), value=apply(all.rates,2,mean), sd=apply(all.rates,2,sd),order=1:10
  )
  x11()
  ggplot(data) + geom_bar( aes(x=reorder(name,order),y=value), stat="identity",fill="skyblue" )  +
    geom_errorbar( aes(x=name, ymin = value-sd,ymax = value+sd ),width=0.4,colour="orange",size=1.3 ) +
    theme(text = element_text(size=20)) + labs(x = "Element in array",y="Freq edited sites")
            #axis.text.x = element_text(angle=90, hjust=1))

}

compare.ecdf<-function(file.path){
  rate.dist.file = paste(file.path,"editRate/distribution_RF_perTree.txt",sep="",
  collapse="");
  rate.dist = read.table(rate.dist.file)
  n = dim(rate.dist)[1]

  for(i in 1:length(names(rate.dist))){
    a = sort(rate.dist[,i])
    if(i ==1){
      plot(a,(1:n)/n,type = "l",col = "red",ylab="Cumulative fraction",xlab="Norm RF score",
    main = "ECDF all MEMOIR colonies")
    }else {
      lines(a,(1:n)/n)

    }
  }

}

#paper plot
#plot the frequency of each edit per site as a stacked barplot
#this plots uses ggplot2
stackedplot.edit.rate<-function(rate.file="/Users/alejandrog/MEGA/Caltech/trees/GraceData/10mer_2019/editRate/All_rates.txt"){

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

histogram.plot.RF<-function(file.path=file.path,posInfo){

  rate.dist.file = paste(file.path,"editRate/distribution_RF_perTree.txt",sep="",collapse="")
  rate.dist = read.table(rate.dist.file)
  tree.analysis = data.frame(RF=rate.dist$V1,name=all.files$file.name[as.logical(all.files$FISH)])
  return(tree.analysis)
}


#plot 4
#read ALL files in the results/ folder and calculate a global rate of edits, including "w" & "R"
stacked.plot.allTrees<-function(which.sites=1:10,results.folder = "/Users/alejandrog/Downloads/results/"){

    all.results = list.files(results.folder)
    translate = T

    alphabet = c("u","r","x")
    N.sites=10

    all.mat.sites = matrix(, nrow=0,ncol=10)
    rate.per.site = c()
    for(i in 1:length(all.results)){
        posInfo=read.table(paste(results.folder,all.results[i],sep=""),header=F, colClasses = 'character')

        names(posInfo)<-c("cell","state")
        if(translate){
            posInfo$state = str_replace_all(posInfo$state, c("2" = "r", "1" = "u","0"="x"))
            posInfo$state = str_replace_all(posInfo$state, c("3" = "u"))
          }


        #edit Jun 10, 2019
      #  posInfo=fread(paste(results.folder,all.results[i],sep=""),header=F)
      #  posInfo = posInfo$V2 #skip the cell ID which is the first column in the file

        barcodes = posInfo$state[!is.na(posInfo$cell) & !posInfo$state==paste(rep("x",barcodeLength),collapse="")]
        names(barcodes) =posInfo$cell[!is.na(posInfo$cell)  & !posInfo$state==paste(rep("x",barcodeLength),collapse="")]
        tt=apply(as.matrix(barcodes),2,strsplit,"")
        mat.sites = t(do.call(cbind,do.call(cbind,tt)))
        #save the average edit rate for each tree:
        #choose which sites to include in the analysis (to look for biases in edit rates)
        #rate.per.site[i] = mean(apply(mat.sites[,which.sites]!="u",2,sum)/dim(mat.sites)[1])

        all.mat.sites= rbind(all.mat.sites,mat.sites)
    }

    library(RColorBrewer); coul = brewer.pal(3,"Pastel2")


    rates = apply(all.mat.sites,2,table)
    rates.pct = rates/apply(rates,2,sum)  * 100
    x11();
    barplot(as.matrix(rates.pct[c(2,3,1),]),col = coul[c(3,2,1)],border = "white",xlab = "barcode")

    return(list(rate.per.site,data,all.mat.sites))
}
#use the object returned above as input for these two functions:
plot.stacked.frequencies<-function(a){

  #making a stacked barplot with ggplot2
  ggplot() +
    geom_bar(aes(y=value,x=siten,fill=letter),data=a[[2]],stat="identity") +
    scale_fill_brewer(palette="Pastel2") +
  #
  # geom_bar(aes(y=value,x=siten,fill=letter),data=a[[2]],stat="identity") +
  #   scale_fill_brewer(palette="Pastel2") +
   labs(x="Site N",y="Frequency",title="Frequency of edits per site (48hrs)")
}
#histogram
plot.histogram.edit.rates<-function(a){

  qplot(a[[1]],geom="histogram",main=paste("Edit rate distribution, n=",toString(length(a[[1]])),sep=""),
  binwidth=0.12,xlab="Fraction of edited sites",ylab="N of colonies",
  fill=I("grey"),col=I("black"))

}

# #including all barcodes in the original data (includes xxxxx, and cells with no Movie.ID)
# x11()
# barplot(prop.table(table(posInfo$Summary)),las=2,ylab="Freq",main="Barcode dist (all cells)")
