# analysis of integrase data
# August 2019

#Load required libraries
library(data.table)
library(phangorn)
library(cluster)
library(stringr)
library(phytools)
library(dplyr)
# Goals:
# barcode.processing function that filters bad cells and bad colonies
# thus creating a data.frame with metadata for good trees to be included
# This is important because we will then simulate the lineages


# Set up the full path for the data
# Barcode data is in the folder FISH/
# Main folder is ../integrase_folder/10mer_2019/
if(length(grep("linux",read.table("../os.txt")$V1))){
    source("../integrase-data/RF.experiment.R")
    integrase_folder= "integrase-data/" # For Ubuntu
    file.path=paste("/home/agranado/MEGA/Caltech/trees/",integrase_folder,"10mer_2019/",sep="") # Ubuntu

}else{
    source("../GraceData/RF.experiment.R")
    integrase_folder = "GraceData/"     # For Mac
    file.path=paste("/users/alejandrog/MEGA/Caltech/trees/",integrase_folder,"10mer_2019/",sep="") # Mac
}

# One functionality is to apply sistematic filters to the data in order to find out the sources of error
# Thus the function needs to filter by cell or by tree and generate a new filtered_dataset

# Global parameters
# How many barcode files do we have? Each one represents a colony
existing.files = list.files(  paste("../",integrase_folder,"10mer_2019/FISH/",sep="") )
#read the xls with all the files info (some are not in the folder)
# This data.frame has names, ncells and ground truth:
all.files=fread(paste(file.path,"results.guide2019.tsv",sep=""))

# Index of files that have FISH data
files_to_process = which(all.files$file.name %in% str_split(existing.files,"\\.",simplify =T)[,1])
# Some colonies have N/A in the ground truth field. So let's remove that
files_to_process = files_to_process[files_to_process != grep("N/A",all.files$newick)]



# # # #
# # # #
# Main Functions
# Sep 4th

# High level functions
# Perform all analysis and call sub-functions







# GLOBAL PARAMETERS
params_global = estim.params.global(estimG = 4, fil = paste("../",integrase_folder,"10mer_2019/editRate/allBarcodes.txt",sep=""))
mu = params_global[[1]]
alpha = params_global[[2]]





# USAGE : process.all.files(files_to_process) %>% reconstruct.all.lineages() -> memoirData


# FIRST MAIN FUNCTION
process.all.files <- function(files_to_process,remove_largeDels = F) {

  aa = lapply(files_to_process,barcode.processing,remove_largeDels)
  aa=do.call(rbind,lapply(aa,unlist))
  colnames(aa)<-c("fileID","fileName","nCells","edit","dels","ground")
  aa<-as.data.frame(aa)

  #data frame is in factor format
  #to make our lives easier we can convert it to numeric
  indx<-sapply(aa,is.factor)
  indx[2] = FALSE #this is character
  indx[6] = FALSE # this is character
  aa[indx] <- lapply(aa[indx], function(x) as.numeric(as.character(x)))

  return(aa)
}

# SECOND main function
reconstruct.all.lineages <-function(dataTree,min.cells = 2,clust.method = "diana", recursive =F,second.clust.method = "ward.D2"){
  dataTree %>% dplyr::filter(nCells>min.cells) -> dataTree


  # dataTree is our first data sctructure
  aa = lapply(1:dim(dataTree)[1], reconstruct.lineages,dataTree, clust.method = clust.method,recursive.reconstruction = recursive, second.clust.method = second.clust.method )
  aa = do.call(rbind,aa)
  aa <- as.data.frame(aa)
  names(aa) <- c("fileName","nCells","edit","dels","RF","ground","rec")

  indx<-sapply(aa,is.factor)
  indx[1] = FALSE #this is character
  indx[6] = FALSE # this is character
  indx[7] = FALSE # this is character
  aa[indx] <- lapply(aa[indx], function(x) as.numeric(as.character(x)))

  return(aa)


}
# Global parameter for reconstructing ALL trees


# Works OK Sep 4th

# barcode.processing
# input: meta.data as .xls and barcodes as .txt
# outptu: tidy data frame with barcodes + ground truth after filtering.
# It has to create new barcode files and new ground.truth
# output folder for new barcode files and meta.data is ./filteredData/

barcode.processing<-function(file.idx,remove_largeDels = F){

        print( paste("processing file ", toString(file.idx), all.files$file.name[file.idx]))
        # Where to save the filtered barcode set
        barcode_path = paste(file.path,"filteredData/",sep="")

        file.name=all.files$file.name[file.idx]
        # file.path is global
        fileName=paste(file.path,"FISH/",file.name,sep="")

        # READ fish data from .txt
        posInfo = fread(paste(fileName,".txt",sep=""),sep="\t",colClasses = 'character')
        names(posInfo)<-c("cell","state")
        # some files have 3 as unknown state. We treat that as unedited:
        posInfo$state = str_replace_all(posInfo$state, c("3" = "1"))


        # This is the number of the cell
        posInfo$cell = as.character(as.numeric(substr(posInfo$cell,5,6)))

        barcodeLength =  nchar(as.character(posInfo$state[1]))

        # filter large deletions i.e. xxxxxxxxxx
        # state has the barcode readout
        # cell has the cell's name (number ID)
        keep_this_cells = !is.na(posInfo$cell) & !posInfo$state==paste(rep("0",barcodeLength),collapse="")

        # remove cells with large deletions: 4 or more deletions 0000
        if(remove_largeDels) keep_this_cells = keep_this_cells & !grepl("000",posInfo$state)
        # NOTE: I need to put some filter such that if we are left without cells after filtering, just skip this colony
        print( paste(toString(sum(keep_this_cells)), "cells included _  ", toString(sum(grepl("000",posInfo$state))) ," large deletions" ))
        #remove: IF (no cells ){ return (NULL ) } else { keep going}
        if(sum(keep_this_cells)>0){
                barcodes = posInfo$state[keep_this_cells]
                names(barcodes) =posInfo$cell[keep_this_cells]

                # WRITE the filtered barcodes to file.
                # This will match the new filtered ground truth for further reconstruction
                # WORKS
                write.table( data.frame(cell = names(barcodes), state=barcodes) , file=paste(file.path,'filteredData/', file.name,'.txt',sep=""), quote=FALSE, sep='\t',row.names = F)


                # Extract the ground truth from the file
                # Here the ground truth has ALL cells, including those with deletions, low-quality etc.
                ground.truth = all.files$newick[file.idx]

                true.tree=read.newick(text=toString(ground.truth))


                # barcodes is already a filtered list. We need to keep only those cells
                # so we filter everything else from the ground truth (dead cells, XXXXXX, etc. )
                alive.tree = filter.cells.groundTruth(barcodes,true.tree)

                #fraction of edited sites

                b=as.matrix(do.call(cbind,strsplit(barcodes,"")))
                d4=sum(b=="1")/prod(dim(b))
                # deletion rate:
                d0=sum(b=="0")/prod(dim(b))

                # FILE_ID   NCELLS    PR_EDIT   GROUND_TRUTH
                return(list(file.idx,file.name, length(barcodes), d4, d0,  write.tree(alive.tree)  ) )
                # WE can now save everything in a new data.frame
                # let.save individual components in
          }else{
            return(NULL)
          }

}



# Input: data frame with the names of all files and ground truth for each
# Output list of ground truth lineages: ready to be reconstructed
# Can be used to get ALL ground_truth lineages like :
# lapply(1:dim(dataTree)[1],reconstruct.lineages,dataTree=dataTree)
reconstruct.lineages <- function(file.idx,dataTree,control =F,clust.method = "diana",recursive.reconstruction = F,second.clust.method = "ward.D2"){


    print( paste( "processing", dataTree$fileName[file.idx] )  )
    posInfo = fread(paste(file.path,"filteredData/",dataTree$fileName[file.idx],".txt",sep=""),colClasses = "character")
    groundTruth = dataTree$ground[file.idx]
    ground_phylo = as.phylo(read.newick(text = toString(groundTruth)))

    # we keep just the cell number, i.e. remove the frame (216) which is no useful
    # For random control we just re-sample the tips
    if(!control){
      sample_tips<-str_split(ground_phylo$tip.label,"_",simplify = T)[,2]
    }else{
      sample_tips<-sample(str_split(ground_phylo$tip.label,"_",simplify = T)[,2])
    }
    # We need the labels to be in the same order as in the tree
    # We match the tips order to the barcode data
    # Then we take from the barcode data using the matching indexes
    ground_phylo$tip.label = paste(posInfo$cell,"_",posInfo$state[match(sample_tips,posInfo$cell)],sep="")

    # This function takes the barcodes fromt he groun truth tree and reconstructs an independent lineage
    manualTree = reconstructLineage(ground_phylo,mu,alpha,return_tree = T,clust.method = clust.method)
    #Calculation of different scoring metrics:
    #
    RF_score = 1- RF.dist(manualTree,ground_phylo,normalize = T)

    # NEw method for DIANA + hclust
    if(recursive.reconstruction)
      RF_score = 1- RF.dist(recursiveReconstruction(manualTree,mu,alpha, second.clust.method),ground_phylo,normalize = T)



    results_row = c(as.character(dataTree$fileName[file.idx]),toString(length(ground_phylo$tip.label)),
                  toString(dataTree$edit[file.idx]), toString(dataTree$dels[file.idx]),  toString(RF_score) ,
                      write.tree(ground_phylo), write.tree(manualTree) )

    return(results_row)



}

compute.distance <- function(index,dataTree){


}
# We have a list of already filtered barcodes and we need to remove everything else from the ground truth a
filter.cells.groundTruth<-function(barcodes,true.tree){

  # At this point barcodes is the list of cells that we want to keep for further analysis

  ground.labels = true.tree$tip.label
  # last recording frame
  # Cell ID has the last frame the sell was seen, so we keep the ones that were seen until the end of the experiment
  last.frame = 216 #we know this from the experiment
  lf.pattern = paste(toString(last.frame),"_",sep="")

  #last frame gives us alive cells while names(barcodes) gives us the cells that we filtered by diverse reasons
  keep.leaves = paste(last.frame,"_",names(barcodes),sep="")

  drop.leaves = ground.labels[-which( ground.labels %in% keep.leaves)]

  new.tree = drop.tip(true.tree,drop.leaves)

  return(new.tree)

}



# this function reconstructs a single lineage
# It takes the barcode from the ground truth tree

reconstructLineage<-function(ground_phylo,mu,alpha,return_tree = F,clust.method = "diana"){

  #get the barcode and cell id (cell id is not neccessarily continuous numbers)
  barcodes = str_split(ground_phylo$tip.label,"_",simplify=T)[,2]
  cell_ids = str_split(ground_phylo$tip.label,"_",simplify=T)[,1]

  # translate the barcode data to the old notation uxr
  barcodes_urx = str_replace_all(barcodes, c("2" = "r", "1" = "u","0"="x"))
  barcodes_urx = str_replace_all(barcodes_urx, c("3" = "u"))

  #recontruct the tree
  matdist_=manualDistML_2(as.character(barcodes_urx),mu = mu ,alpha = alpha,nGen = 4 )

  matdist_2 = cousinDistML(as.character(barcodes_urx),mu = mu ,alpha = alpha,nGen = 4 )

  row.names(matdist_)<- paste(cell_ids,barcodes,sep="_")
  colnames(matdist_)<- paste(cell_ids,barcodes,sep="_")

  hclust.tree = clusterDistMatrix(matdist_,clust.method)

  #here ground_phylo is "the ground truth" which is a simulation but this is what we want
  d = 1- RF.dist(hclust.tree,ground_phylo,normalize =T)
  if(return_tree){
    return(hclust.tree)
  }else{
    return(d)
  }

}

# In principle, we can use any clustering method to recover a dendrogram from the
# distance matrix. Divisive clustering (diana) show higher performance than other clustering method
# D2.ward & complete_linkage methods also perform well
clusterDistMatrix<-function(matdist_,clust.method = "diana"){

  if(clust.method =="diana"){
    hclust.tree = as.phylo(as.hclust( diana(as.dist(t(matdist_)))))
  }else {

    hclust.tree=as.phylo(hclust(as.dist(t(matdist_)),method = clust.method))
  }

    return(hclust.tree)
}








# Sep 9th 2019
# IDEA:
# use DIANA to reconstruct the early lineage
# cut the tree at 3 groups and then reconstruct each subgroup with hclust (which is best for sister and small groups)
library(dendextend)


recursiveReconstruction<-function(ground_phylo,mu,alpha, second.clust.method= "ward.D2",min.tree.size = 9){

  if(length(ground_phylo$tip.label)>min.tree.size){
      ngroups = 3
      #get the DIANA tree first
      diana_tree<-reconstructLineage(ground_phylo,mu,alpha,T,clust.method="diana")

      diana_labels = cutree(diana_tree,k=ngroups) #returns cluster labels for each cell
      # CONVERT  the whole tree to dendrogram class
      # ground_dendro = as.dendrogram.phylo(diana_tree)


       # FOR 1:3   (assuming 3 main clades)
       # reconstruct the first partition of DIANA tree using hclust complete (or Ward D2)
       # NOTE Only for the first patition:
       recursive_phylo = findCladeRecursive(diana_tree, diana_labels,mu,alpha,second.clust.method )


       # Then we have to find out, which is the other main partition
       # We have the labels, so we might need to explore

       ##xx$tip.label<-sample(xx$tip.label)

       # NOTE: do something to xx (i.e. reconstruct)

       # convert back to phylogram and replace
       ground_dendro[[1]]<-as.dendrogram.phylo(xx)
      # x11();plot(ground_dendro)

       return(recursive_phylo)

   }else{
       #do nothing
       return(ground_phylo)
   }

}


findCladeRecursive<-function(ground_phylo,diana_labels,mu,alpha,second.clust.method = "ward.D2"){

    # This function will use the drndrogram object which can be easily accessed
    ground_dendro = as.dendrogram.phylo(ground_phylo)

    # This function will only reconstruct the first level of the DIANA tree
    # this structure comes as a list so we could take the principal branches directly
    # ground_dendro[[1]] and ground_dendro[[2]] are the two main clades (this is a binary tree)
    main_clade1 = as.phylo(ground_dendro[[1]])
    main_clade2 = as.phylo(ground_dendro[[2]])


    # index_clade = c()
    #
    # for(i in 1:ngroups){
    #   index_clade[i]= length(which(main_clade2$tip.label %in% names(diana_labels)[diana_labels==i] ))>0
    # }


    #quick and dirty try for 2 groups
    for(j in 1:2){
      #test for leaf:
      if(length(ground_dendro[[j]])<2)
        next # don't processs this branch since it is a leaf

      main_clade =   as.phylo(ground_dendro[[j]])
      #how many leaves on this clade:
      # since we are going to re-reconstruct this clade, we need at least 3 leaves
      if(length(main_clade$tip.label)>2){
        main_clade_clust  = reconstructLineage(main_clade, mu,alpha, return_tree = T, clust.method = second.clust.method)
      }else{
        main_clade_clust = main_clade
      }

      ground_dendro[[j]]<- as.dendrogram.phylo(main_clade_clust)
    }
    #return phylo object
    # as.Node converts the tree to Node object and somehow fixes weird branching left by the clade replacement
    return(as.phylo(as.Node(ground_dendro)))

}



  #
  # for(i in 1:ngroups){
  #   subclade = extract.clade(diana_tree,MRCA(diana_tree,names(diana_labels[which(diana_labels==1)])))
  #
  #   # this command removes the tips we indicated and returns the rest of the tree
  #   outgroup = dendextend::prune(ground_dendro,leaves =  names(diana_labels[which(diana_labels==1)]))
  #
  #   # reconstruct subclade
  #
  #     node.leaves(phylogeny, node)
  # }

#}
