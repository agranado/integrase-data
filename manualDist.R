manualDist <- function(barcodeLeaves,mu,alpha,nGen){
   alphabet = array()
   alphabet[1]= "u"
   alphabet[2]= "r"
   alphabet[3]= "x"
   #mu & alpha indicate an explicit probabilistic model
   #Probability of no mutation for nGen cell divisions (mu is rate of edit per cell division)

   #transition probabilities
   #beacuse transitions only happen from u ->  then this is a 1-D vector
   Tran.pr = c(1, alpha,1-alpha)


   #NULL MODEL: probability of observing sites as independent events:
   Pr = array()
   Pr[1] = (1-mu)^nGen
   #probability of nu mutation during nGen-1 devisions and then a mutation in generation nGen times Pr(alphabet[2])
   #then we use the choose to correct for all the order in which this could have happened, (we need still further correction to
   #to account for irreversibility)
   #Pr[2] = choose(nGen,nGen-1)*( 1-mu)^(nGen-1)*mu*alpha
   #corrected for irreversibility:
   Pr[2] = Pr_edit(nGen,mu,alpha)

   #same as before but using (1-alpha)
   #Pr[3] = choose(nGen,nGen-1)*(1-mu)^(nGen-1)*mu*(1-alpha)
   Pr[3] = Pr_edit(nGen,mu,1-alpha)

   PrMatrix  = array(0,dim=c(length(Pr),length(Pr)))
   #calcualte probabilistic model:
   #this just means the pr that those two sites have those characters by random. This is the expected pr for each site
   #it assummes independence but does not tell you how likely they are to come from a common ancestor
   for (p1 in 1:length(alphabet)){
     for (p2 in 1:length(alphabet)){
        PrMatrix[p1,p2] = Pr[p1] * Pr[p2]
     }
   }

   #weights for number of sustitutions
   equalU = 0 #barcodeLength #both cells with u IMPLY ancestry
   equalSust = 0 #both cells with
   oneSust = 1
   twoSust = 2

   nBarcodes = length(barcodeLeaves)
   barcodeLength=nchar(barcodeLeaves[1])
   distMat= array(0,dim =c(nBarcodes,nBarcodes))
   ratioMat = array(0,dim=c(nBarcodes,nBarcodes))
   productMat = array(0,dim=c(nBarcodes,nBarcodes))

    #go through all the elements in the barcode array
    for (i in 1:(nBarcodes-1)){
      for (j in (i+1):nBarcodes){
        barcodeArray1 =strsplit(barcodeLeaves[i],"")[[1]]
        barcodeArray2 =strsplit(barcodeLeaves[j],"")[[1]]
        distSum = 10
        ratio.sum =0
        ratio.product=1
        #for each pairwise comparison
        for (s in 1:barcodeLength){
            #Pr of observing these characters independtly arising Pr1 * Pr2 : assuming independence at nGen
            Pr.sust = PrMatrix[which(alphabet ==barcodeArray1[s]),which(alphabet ==barcodeArray2[s])] #this is just the product of both Pr
            Pr.sust.inv = 1/Pr.sust
            #both characters are the same (independently of their identity)
            if(barcodeArray1[s]==barcodeArray2[s]){
                #calcualte probability of both sites

                #Pr_sust
                #probabilities under assumption of sister cells
                if(barcodeArray1[s]=="u"){
                  distSum = distSum - equalU
                  #probability of u in the previous generations times pr(no sust) * pr(no sust)
                  Pr_sister=(1-mu)^(nGen-1) * (1-mu)^ 2
                }else if(barcodeArray1[s]=="r"){
                  #if both are r : probability that ancestor is r + pr_a (u) * Pr(sust) ^2
                  # Pr(r_{t-1}) + Pr(u_{t-1}) * Pr(u->r)^2 = Pr(u->r_{t},u->r_{t} | u_{t-1})
                  Pr_sister=Pr_edit(nGen-1,mu,alpha) + (1-mu)^(nGen-1) * (mu *alpha)^2

                  distSum = distSum - equalSust
                }else if(barcodeArray1[s]=="x"){
                  #if both are x : probability that ancestor is r + pr_a (u) * Pr(sust) ^2
                  Pr_sister=Pr_edit(nGen-1,mu,1-alpha) + (1-mu)^(nGen-1) * (mu * (1-alpha))^2

                  distSum = distSum - equalSust
                }
                  #ratio between sister probability and random probability
                  ratio.sum = ratio.sum + Pr_sister/Pr.sust
                  ratio.product = ratio.product* Pr_sister/Pr.sust
            }else{
              #characteres have different sites.
              #is any of the character a u
              b = grepl(alphabet[1],c(barcodeArray1[s],barcodeArray2[s]))

              #one of them is u
              if(length(which(b==FALSE))==1){
                   distSum = distSum + oneSust *Pr.sust
                   #are there any r
                   c = grepl(alphabet[2],c(barcodeArray1[s],barcodeArray2[s]))
                   #there is one r
                   if(length(which(c==TRUE))==1){
                      Pr_sister=(1-mu)^(nGen-1) * (1-mu) * mu * Tran.pr[2]

                   }else {
                     #it is an x:
                     # Pr_sister = Pr(u_{t-1}) * Pr(u->u_{t}) * Pr(u->r_{t})
                     Pr_sister=(1-mu)^(nGen-1) * (1-mu) * mu * Tran.pr[3]

                   }
                   ratio.sum = ratio.sum + Pr_sister/Pr.sust
              #none of them is u AND they are different
              }else if(length(which(b==FALSE))==2){
                  distSum = distSum + twoSust *Pr.sust

                  #Pr_sister= Pr(r,x | u_{t-1})
                  Pr_sister = (1-mu)^(nGen-1) * mu^2 * Tran.pr[2] * Tran.pr[3]
                  ratio.sum  = ratio.sum + Pr_sister/Pr.sust
              }
              #sister probabilities:
            }
            ratio.product = ratio.product* Pr_sister/Pr.sust
        } #end barcode Length loop (finished analysis these sisters)

    distMat[i,j]= distSum
    ratioMat[i,j]=1/ratio.sum  * distSum
    productMat[i,j] = 1/ratio.product
    }
  }
return(ratioMat)
}
