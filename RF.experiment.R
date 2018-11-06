    #function
    #distance between ground truth and reconstructed trees

    #this function takes the reconstructed tree and the ground truth to calculated;
    #the RF distance. Normalization comes from the distribution of RF distnaces
    #for random trees (with the same numner of leaves. )

rand.RF.dist<-function(new.tree){

  n.reps = 100
  rand.dist =c()
  for(i in 1:n.reps){
    random.tree = rcoal(length(new.tree$tip.label),tip.label = new.tree$tip.label)
    rand.dist[i] = RF.dist(new.tree, random.tree)
  }
  return(rand.dist)
}

normalized.RF.experiment<-function(true.tree, manualTree, barcodes,posInfo){
    library(ape)


    #Calculate the distance of all leaves to the root of the trees
    leaves.length= diag(vcv.phylo(true.tree))
    ground.labels = true.tree$tip.label

    #specific for pos40

    #barcodes = barcodes[-length(barcodes)]
    method = 2
    # keep only the cells that stayed until the end of the Movie
    #In theory, the max value of length should correspond to  the last movie frame
    #Some cells, though, have max-1 (or something like that)
    if(method==1){
      length.threshold = max(leaves.length)-1
      drop.leaves = leaves.length[ leaves.length < length.threshold]
      keep.leaves = leaves.length[ leaves.length >= length.threshold]
    #remove dead cells:
      #method 1 : more systematic
      #identify the lenght of the branches then remove dead cells
      new.tree = true.tree
      for(i in 1:length(drop.leaves)){
          new.tree = drop.tip(new.tree,names(drop.leaves)[i])
      }
    }else if(method==2){
      #method 2: more biologist-style
      #we know the lenght of the branches so just match the name (Grace put that!)
      last.frame = 216 #we know this from the experiment
      lf.pattern = paste(toString(last.frame),"_",sep="")
      #keep.index contains the cells that are alive  (i.e. have FISH readout)
      #however, still contains cells that have xxxx.. in their data. So we need to remove them
      keep.index = grep(lf.pattern,ground.labels)

      #remove XXXX
      good.cells = posInfo$cell[!is.na(posInfo$cell) & !posInfo$state==paste(rep("x",barcodeLength),collapse="")]
      bad.cells = posInfo$cell[!(!is.na(posInfo$cell) & !posInfo$state==paste(rep("x",barcodeLength),collapse=""))]

      keep.leaves = paste(lf.pattern,good.cells,sep="")

      keep.them = c();for(h in 1:length(keep.leaves)){keep.them[h] = grep(paste(keep.leaves[h],"$",sep=""),ground.labels)}
      drop.leaves = ground.labels[-keep.them]

      # if(length(bad.cells>0)){
      #   drop.leaves.x = paste(lf.pattern,bad.cells,sep="")
      #   drop.leaves = c(ground.labels[-grep(lf.pattern,ground.labels)],drop.leaves.x)
      # }else{
      #
      #     drop.leaves = ground.labels[-grep(lf.pattern,ground.labels)]
      # }


      #keep.leaves = ground.labels[grep(lf.pattern,ground.labels)]

      #EDIT: we can take the cells from the excel file directly since there is only FISH for cells that
      #make it to the end fo the movie, so the filtering for dead cells was already done by GRACE

      new.tree = true.tree
      for(i in 1:length(drop.leaves)){
          new.tree = drop.tip(new.tree,drop.leaves[i])#drop.leaves are the names already
      }
    }

    #THE NEW labels for the alive cells
    alive.cells = new.tree$tip.label
    new.alive.cells = array(0,length(alive.cells))

    id.n=names(barcodes)
    for(i in 1:length(barcodes)){
        this.cell= grep(paste(lf.pattern,id.n[i],"$",sep=""),alive.cells)
        new.alive.cells[this.cell] = paste(names(barcodes)[i],"_",  toString(barcodes[i]),sep="")
    }


    tip.names = new.alive.cells
    for(tt in 1:length(tip.names)){
      this.cell.n = str_extract(tip.names[tt],"\\d+")
      if(nchar(this.cell.n)==1){
        tip.names[tt]=paste("0",tip.names[tt],sep="")
      }
    }

    tip.names=gsub("r","2",tip.names)
    tip.names=gsub("x","1",tip.names)
    tip.names=gsub("u","0",tip.names)


    new.tree$tip.label<-tip.names

    #calculate the distribution of RF distances using random clustering based on the number of leaves:
    n.reps = 100
    rand.dist =c()
    for(i in 1:n.reps){
      random.tree = rcoal(length(new.tree$tip.label),tip.label = new.tree$tip.label)
      rand.dist[i] = RF.dist(new.tree, random.tree)
    }
    if(mean(rand.dist)==0){
      #The case for unrooted trees of size 3
      this.score = 1 #kind of cheating, but is true
    }else{
      this.score = 1-RF.dist(manualTree,new.tree)/mean(rand.dist)
    }

    if(this.score<0){this.score=0}
    return(list(this.score,new.tree,mean(rand.dist)))

  }
file.path.membow = file.path
file.path.memoir = "/Users/alejandrog/MEGA/Caltech/trees/GraceData/10mer/"

normalized.RF.membow<-function(true.tree,manualTree,barcodes,posInfo){


  # 1. Find clusters based on identical genotypes
      unique.genotypes = unique(substr(manualTree$tip.label,4,14))
      if(length(unique.genotypes)>2){
      # 2. Reconstruct based on the new clusters of identical cells
      #      {What is the right order for 1 & 2?}
          RF.list = normalized.RF.experiment(true.tree,manualTree,barcodes,posInfo)
          alive.tree =RF.list[[2]]
          rand.dist = RF.list[[3]]
          #simplify the reconstructed tree, by removing duplicated genotypes
          new.tree = manualTree
          new.alive.tree = alive.tree
          for(g in 1:length(unique.genotypes)){
              drop.leaves = manualTree$tip.label[grep(unique.genotypes[g],manualTree$tip.label)]
              #delete all except one

              if(length(drop.leaves)>1){
                  for(i in 1:(length(drop.leaves)-1)){
                      new.tree = drop.tip(new.tree,drop.leaves[i])
                      new.alive.tree = drop.tip(new.alive.tree,drop.leaves[i])
                  }
              }

          }

  # 3. Using the ground truth, remove cells that were missclassified, such that _groundTruth
  # and reconstruction are comparable.

  # 4. Somehow map the ground truth to the dendrogram by excluding cells that
  # were missclassified from both trees (make the clades identical)
          this.rand =mean(rand.RF.dist(new.tree))
          if(this.rand<2){ this.rand=2}

          this.score = 1-RF.dist(new.alive.tree,new.tree)/this.rand
    }else{
          this.score = 1

    }
  # 5.
  return(list(this.score,length(unique.genotypes)))
}
