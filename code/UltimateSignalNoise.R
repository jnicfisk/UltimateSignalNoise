TreeCollapse.signal.noise<-function(freqs,subratevector,ratevector,thisTree=NULL,treePath=NULL, center=NULL){
  ##freqs in order T,C,A,G
  ##subratevector in order a,b,c,d,e,f
  require(ape) #ape is good for tree structures
  require(phylobase) #not sure if actually used now
  require(phytools)
  if(is.null(treePath)&is.null(thisTree)){
    print("Must provide an ape tree (thisTree) or a path to a .tre tree (treePath)")
    return(NULL)
  }
  if(!is.null(treePath)&!is.null(thisTree)){
    print("Both a tree and a path to a tree have been input. Defaulting to using the tree object (thisTree)...")
    treePath<-NULL
    #return(NULL)
  }
  
  if(!is.null(treePath)){
    if (!file.exists(treePath)) {
      print("Path to input tree does not exist")
      return(NULL)
    }
  }
  tree<-thisTree
  if(!is.null(treePath)){
    tree<-read.tree(treePath)
  }

  
  
  if(sum(freqs)!=1){
    print("Base Frequencies must sum to 1")
    return(NULL)
  }
  
  if(length(subratevector)!=6){
    print("subratevector must have length 6")
    print("Error in values a-f")
    return(NULL)
  }
  
  #change when packaging
  source("signal_noise_functions.R")
  
  #unload command args
  a <- subratevector[1]
  b <- subratevector[2]
  c <- subratevector[3]
  d <- subratevector[4]
  e <- subratevector[5]
  f <- subratevector[6]
  
  piT <- freqs[1]
  piC <- freqs[2]
  piA <- freqs[3]
  piG <- freqs[4]
  
  frequ<-c(piT,piC,piA,piG)
  ####ACTUAL PROGRAM STARTS HERE
  
  ####Get Shared EigenInformation
  
  eigenStuff <-
    doEigenMaths(calcTrueQ(
      calcQ(a, b, c, d, e, f, piT, piA, piG, piC),
      calcMu(a, b, c, d, e, f, piT, piA, piG, piC)
    ))
  
  btree<-as(tree,"phylo4")
  numTips<-nTips(btree)
  #allsibs<-getListOfSiblings(btree)
  allTipCombos<-combn(1:nTips(btree),2)
  #allTipCombos<-removeSibs(allTipCombos,allsibs)
  allNodesMRCA<-getAllNodeMatrix(tree,btree,allTipCombos)
  allNodesNext<-matrix(ncol=ncol(allNodesMRCA),nrow=nrow(allNodesMRCA))
  colnames(allNodesMRCA)<-tree$tip.label
  rownames(allNodesMRCA)<-tree$tip.label
  colnames(allNodesNext)<-tree$tip.label
  rownames(allNodesNext)<-tree$tip.label
  
  
  
  allInterPairs<-list()
  
  uniqueInternodes<-unique(as.vector(allNodes))
  #each unique internode is the larger of the two. So internode and the ancestor of the internode define the true internode
  uniqueInternodes<-uniqueInternodes[which(!is.na(uniqueInternodes))]
  uniqueRates<-unique(ratevector)
  masterMat<-matrix(nrow=length(uniqueRates),ncol=length(uniqueInternodes))
  colnames(masterMat)<-uniqueInternodes
  rownames(masterMat)<-uniqueRates
  for(lam in uniqueRates){
    clam<-as.character(lam)
    for(uI in uniqueInternodes){
      cuI<-as.character(uI)
    #split the tree 
    #collapse each half
    #rejoin the tree
    
    #split tree in 2
    internodePair<-c(ancestor(btree,uI),uI)
    internodeLen<-getEdgeBetweenTwoNodes(btree,internodePair)
    #not sure this is right splitting of trees (IT ISN'T)
    subtree2<-keep.tip(tree,tip = names(descendants(btree,internodePair[2])))
    subtree1<-keep.tip(tree,tree$tip.label[which(!tree$tip.label%in%subtree2$tip.label)])
    plot(subtree1)
    while(length(subtree1$edge.length) > 2) {
      subtree1<-collapseNextPair(subtree1,freqs,lam,eigenStuff)
      plot(subtree1)
    }
    plot(subtree2)
    while (length(subtree2$edge.length) > 2) {
      subtree2<-collapseNextPair(subtree2,freqs,lam,eigenStuff)
      plot(subtree2)
    }
    q1<-subtree1$edge.length
    q2<-subtree2$edge.length
    quartet<-read.newick(text = paste0("(Q1:",q1[1],",Q2:",q1[2],",(Q3:",q2[1],",Q4:",q2[2],"):",internodeLen,");"))
    Info_mes<-processQuartet(ratevector,eigenStuff,quartet)
    allInfoRes[[length(allInfoRes)+1]]<-Info_mes
    names(allInfoRes)[length(allInfoRes)]<-uI
    }
  }

  
  
  
  
   ###
  
  processQuartet<-function(ratevector,eigenStuff,quartet){
    #Initialize blanks
    eYsum <- 0
    eX1sum <- 0
    eX2sum <- 0
    eY2sum <- 0
    eX12sum <- 0
    eX22sum <- 0
    eX1Ysum <- 0
    eX2Ysum <- 0
    eX1X2sum <- 0
    
    for (lmbda in ratevector) {
      all <- evalLambda(lmbda,eigenStuff,quartet)
      y <- all[1]
      x1 <- all[2]
      x2 <- all[3]
      eYsum <- eYsum + y
      eX1sum <- eX1sum + x1
      eX2sum <- eX2sum + x2
      
      eY2sum <- eY2sum + (y ^ 2)
      eX12sum <- eX12sum + (x1 ^ 2)
      eX22sum <- eX22sum + (x2 ^ 2)
      
      eX1Ysum <- eX1Ysum + (x1 * y)
      eX2Ysum <- eX2Ysum + (x2 * y)
      eX1X2sum <- eX1X2sum + (x1 * x2)
    }
    
    Mu_1 <- eYsum - eX1sum
    Mu_2 <- eYsum - eX2sum
    
    
    Sigma_1 <-
      sqrt(eX1sum + eYsum - eX12sum - eY2sum + 2 * eX1Ysum)
    Sigma_2 <-
      sqrt(eX2sum + eYsum - eX22sum - eY2sum + 2 * eX2Ysum)
    Rho_ <-
      (-eX1X2sum + eX1Ysum + eX2Ysum + eYsum - eY2sum) / (Sigma_1 * Sigma_2)
    
    #Internal function for integration
    FofT <- function(t) {
      F1ofT = ((1 / Sigma_1) * dnorm((t - Mu_1) / Sigma_1) * pnorm(Rho_ * (t - Mu_1) /
                                                                     (Sigma_1 * sqrt(1 - Rho_ * Rho_)) - (t - Mu_2) / (Sigma_2 * sqrt(1 - Rho_ *
                                                                                                                                        Rho_))))
      F2ofT = ((1 / Sigma_2) * dnorm((t - Mu_2) / Sigma_2) *
                 pnorm(Rho_ * (t - Mu_2) / (Sigma_2 * sqrt(1 - Rho_ * Rho_)) - (t - Mu_1) /
                         (Sigma_1 * sqrt(1 - Rho_ * Rho_))))
      return(F1ofT + F2ofT)
    }
    
    princtree <- integrate(FofT,-Inf,-.5)
    prpolytomy = integrate(FofT,-.5, .5)
    prcortree  = integrate(FofT, .5, Inf)
    
    #print(paste0("Probablility Correct: ", prcortree$value))
    #print(paste0("Probability Incorrect: ", princtree$value))
    #print(paste0("Probability Polytomy: ", prpolytomy$value))
    return(c(prcortree$value,princtree$value,prpolytomy$value))
  }
  
  
  
  return(allInfoRes)
  
  #Internal function to evaluate lamda
  #should pass in the internode pair or subtrees here. 
  
}

#Independent test call
#TreeCollapse.signal.noise(freqs=c(0.34,0.16,0.32,0.18),subratevector=c(5.26,8.15,1,2.25,3.16,5.44),ratevector=c(rep(0.05,20),rep(0.1,10),rep(0.05,20)),treePath="~/GitHub/Ultimate_Signal_Noise/test.tre")

#Dependent tet calls
#TreeCollapse.signal.noise(freqs=c(0.34,0.16,0.32,0.18),subratevector=c(5.26,8.15,1,2.25,3.16,5.44),ratevector=new,treePath="test.tre")
# set.seed(66)
# for(i in 1:5){
#    print("i")
#    print(i)
#    write.tree(rtree(n=20*i,rooted = FALSE),file = "rtree.tre")
#    for(j in 1:5){
#      new<-runif(n=50*j, min=0, max=.7)
#      print("j")
#      print(j)
#      print((TreeCollapse.signal.noise(freqs=c(0.34,0.16,0.32,0.18),subratevector=c(5.26,8.15,1,2.25,3.16,5.44),ratevector=new,treePath="rtree.tre",specificNode=60)))
#    }
#  }
