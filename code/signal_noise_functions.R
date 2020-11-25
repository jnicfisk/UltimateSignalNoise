#signal_noise_functions
require(ape) #ape is good for tree structures
require(phylobase) #not sure if actually used now
require(phytools)
##Shared internal functions
calcMu <- function(a, b, c, d, e, f, piT, piA, piG, piC) {
  Mu_ <- (1 / 2) / (a * piT * piC + b * piT * piA + c * piT * piG +
                      d * piC * piA + e * piC * piG + f * piA * piG)
  return(Mu_)
}
#calc Q substitution rate matrix
calcQ <- function(a, b, c, d, e, f, piT, piA, piG, piC) {
  #Dim 4X4
  Q <- matrix(nrow = 4, ncol = 4)
  Q[1, 1] <- ((-a) * piC) - (b * piA) - (c * piG)
  Q[1, 2] <- a * piC
  Q[1, 3] <- b * piA
  Q[1, 4] <- c * piG
  Q[2, 1] <- a * piT
  Q[2, 2] <- ((-a) * piT) - (d * piA) - (e * piG)
  Q[2, 3] <- d * piA
  Q[2, 4] <- e * piG
  Q[3, 1] <- b * piT
  Q[3, 2] <- d * piC
  Q[3, 3] <- ((-b) * piT) - (d * piC) - (f * piG)
  Q[3, 4] <- f * piG
  Q[4, 1] <- c * piT
  Q[4, 2] <- e * piC
  Q[4, 3] <- f * piA
  Q[4, 4] <- ((-c) * piT) - (e * piC) - (f * piA)
  return(Q)
}
#calcTrueQ is an internal function, just multiplies Mu_*Q. Returns final Q matrix
calcTrueQ <- function(Q, Mu_) {
  return(Mu_ * Q)
}
#function to get the eigenvalues, vectors, and inverse of the vectors
doEigenMaths <- function(Q) {
  #Obtain the eigenvalues and vectors
  evects <- eigen(Q)
  evalues <- evects$values
  evectors <- evects$vectors
  
  #Reorder in ascending order, swap rows.
  evalues <- evalues[c(4, 3, 2, 1)] #same as mathematica
  evectors <- evectors[c(4, 3, 2, 1), c(4, 3, 2, 1)]
  tev <- evectors
  #Get inverse
  itev <- solve(tev)
  return(list(evalues, evectors, itev))
}
#### END SHARED INTERNAL FUNCTIONS

####Start Collapse Internal Functions

#gets probability transition matrix
P <- function(lamda, T_, eigenStuff) {
  evalues <- eigenStuff[[1]]
  tev <- eigenStuff[[2]]
  itev <- eigenStuff[[3]]
  return(tev %*% (diag(exp(
    evalues * lamda * T_
  )) %*% itev))
}
#get the joint probaility for the first conditional probability
JointProb3 <-
  function(freqs,
           root,
           char1,
           char2,
           T1,
           T2,
           lamda,
           eigenStuff) {
    return(freqs[root] * P(lamda, T1, eigenStuff)[root, char1] * P(lamda, T2, eigenStuff)[root, char2])
  }
#get the conditional probability for the first conditional entropy
CondProb3 <-
  function(freqs,
           root,
           char1,
           char2,
           T1,
           T2,
           lamda,
           eigenStuff) {
    output <- JointProb3(freqs, root, char1, char2, T1, T2, lamda, eigenStuff)
    outhold <- 0
    for (i in 1:length(freqs)) {
      outhold <-
        outhold + JointProb3(freqs, i, char1, char2, T1, T2, lamda, eigenStuff)
    }
    return(output / outhold)
  }
#Calculate the conditional entropy for 2 branches at average rate lamda
CondEntropy3 <- function(freqs, T1, T2, lamda, eigenStuff) {
  output <- 0
  len_freq <- 1:length(freqs)
  for (root in len_freq) {
    for (char1 in len_freq) {
      for (char2 in len_freq) {
        JP3 <- JointProb3(freqs, root, char1, char2, T1, T2, lamda, eigenStuff)
        logCP3 <-
          log(CondProb3(freqs, root, char1, char2, T1, T2, lamda, eigenStuff))
        output <- output + (JP3 * logCP3)
      }
    }
  }
  return(-1 * output)
}

#Compute the joint probability for the second conditional probability
JointProb2 <-
  function(freqs,
           rootprime,
           charprime,
           Tprime,
           lamda,
           eigenStuff) {
    return(freqs[rootprime] * P(lamda, Tprime, eigenStuff)[rootprime, charprime])
  }
#compute the conditional prob for the second conditional entropy
CondProb2 <-
  function(freqs,
           rootprime,
           charprime,
           Tprime,
           lamda,
           eigenStuff) {
    hold1 <- JointProb2(freqs, rootprime, charprime, Tprime, lamda, eigenStuff)
    tmp <- 0
    for (i in 1:length(freqs)) {
      tmp <- tmp + JointProb2(freqs, i, charprime, Tprime, lamda, eigenStuff)
    }
    return(hold1 / tmp)
  }
#get the second part of the conditional entropy
CondEntropy2 <- function(freqs, Tprime, lamda, eigenStuff) {
  output <- 0
  len_freq <- 1:length(freqs)
  for (root in len_freq) {
    for (charprime in len_freq) {
      cprob <- log(CondProb2(freqs, root, charprime, Tprime, lamda, eigenStuff))
      output <-
        output + JointProb2(freqs, root, charprime, Tprime, lamda, eigenStuff) *
        cprob
    }
  }
  return(-1 * output)
}
#get mutual information/entropy for collapse by finding root
InfoEquivalent <- function(freqs, T1, T2, lamda, eigenStuff) {
  cond_entropy3 <- CondEntropy3(freqs, T1, T2, lamda, eigenStuff)
  #internal function to minimize (pick Tprime such that
  # cond_entropy2-cond_entropy3 is 0)
  to_find_root <- function(Tprime) {
    cond_entropy2 <- CondEntropy2(freqs, Tprime, lamda, eigenStuff)
    return(cond_entropy2 - cond_entropy3)
  }
  #output<-uniroot(to_find_root,min(c(T1,T2)),maxiter = 500)
  #output<-uniroot(Vectorize(to_find_root),c(0.000000001,min(T1,T2)),maxiter = 10000)$root
  output <-
    uniroot((to_find_root), c(0.000000001, min(T1, T2)), maxiter = 10000)$root
  
  return(output)
}

#there isn't a good way to collapse the tree branches natively in ape
#this function helps get around that!
bind.tip <- function(tree,
                     tip.label,
                     edge.length = NULL,
                     where = NULL) {
  if (is.null(where))
    where <- length(tree$tip) + 1
  tip <- list(
    edge = matrix(c(2, 1), 1, 2),
    tip.label = tip.label,
    edge.length = edge.length,
    Nnode = 1
  )
  class(tip) <- "phylo"
  obj <- bind.tree(tree, tip, where = where)
  return(obj)
}
#using the current tree, collapse the outmost nested pair of branches
getInfoEquiv <- function(freqs, tree, lamda, eigenStuff) {
  #for now, use average+0.1+nodeabove
  #this needs to be the infocollapse
  
  #debug
  #print(tree$tip.label)
  ##end debug
  
  this <-
    InfoEquivalent(freqs,
                   tree$edge.length[1],
                   tree$edge.length[2],
                   lamda,
                   eigenStuff)
  #this<-mean(tree$edge.length[1]+tree$edge.length[2])+0.1
  nodeabove = tree$edge[1, 1] #do checks
  findCon = which(tree$edge[, 2] == nodeabove)
  addLen <- tree$edge.length[findCon]
  this <- this + addLen
  return(c(this, tree$edge[findCon, 1]))
}
#//////END INTERNAL FUNCTION DEFINITIONS////////

#functions for iterating over all internodes
getEdgeBetweenTwoNodes<-function(tree,nodePair){
  smNode<-min(nodePair)
  biNode<-max(nodePair)
  nodeLen<-edgeLength(btree)[paste0(smNode,"-",biNode)]
  if(is.na(nodeLen)){
    nodeLen<-edgeLength(btree)[paste0(biNode,"-",smNode)]
  }
  if(is.na(nodeLen)){
    print("invalid node pair")
    return("NA")
  }
  return(nodeLen)
}

getListOfSiblings<-function(btree){
  allsibs<-lapply(1:nTips(btree),function(x,btree){return(siblings(btree,x,include.self = T))},btree=btree)
  allsibs<-lapply(allsibs,function(x,btree){return(x[which(x<=nTips(btree))])},btree=btree)
  return(unique(allsibs[which(unlist(lapply(allsibs,length))==2)]))
}

removeSibs<-function(allTipCombos,allsibs){
  for(pair in allsibs){
    allTipCombos<-allTipCombos[,!apply(allTipCombos,2,identical,y=as.vector(pair))]
  }
  return(allTipCombos)
}

#might be able to save time by figuring out which sibling pairs are equivalent
getAllNodeMatrix<-function(ttree,btree,allTipCombos){
  allNodes<-matrix(ncol=nTips(btree),nrow=nTips(btree))
  for(i in 1:ncol(allTipCombos)){
    x<-allTipCombos[1,i]
    y<-allTipCombos[2,i]
    allNodes[x,y]<-MRCA(ttree,c(x,y))
    allNodes[y,x]<-allNodes[x,y]
  }
  return(allNodes)
}
collapseNextPair<-function(tree,freqs,lamda,eigenStuff){
  x<-tree
  x = reorder(x, "postorder") #sort
  #get branch 1 and 2
  first <- x$edge[1, 2]
  second = x$edge[2, 2]
  #return branch location info and values for collapse
  pair <- getInfoEquiv(freqs, x, lamda, eigenStuff)
  #using the value in pair, collapse the tree one level
  x2 <-
    bind.tip(
      x,
      tip.label = paste0(first, "_", second),
      edge.length = pair[1],
      where = pair[2]
    )
  x3a <- drop.tip(x2, tip = second)
  x3b <- drop.tip(x3a, tip = first)
  x <- x3b
  return(x)
}

evalLambda <- function(lamda,eigenStuff,quartet) {
  evalues <- eigenStuff[[1]]
  tev <- eigenStuff[[2]]
  itev <- eigenStuff[[3]]
  #x <- tree
  iter = 0 #to keep track of writing out files (unique names)
  
  internodeIndex <-
    which(quartet$edge[, 2] == 6) #definition of internode (node 5 to node 6)
  quartet_branches <- quartet$edge.length[-internodeIndex]
  internode <- c(quartet_branches, quartet$edge.length[internodeIndex])
  p <- list()
  p <- array(, dim = c(5, 4, 4))
  for (v in 1:length(internode)) {
    p[v, , ] <- (tev %*% (diag(exp(
      evalues * lamda * internode[v]
    )) %*% itev))
  }
  correct <- 0
  wrong1 <- 0
  wrong2 <- 0
  itera<-1:4
  for (original_character in itera) {
    for (internode_character in itera) {
      for (leaf_character_1 in itera) {
        for (leaf_character_2 in itera) {
          if (leaf_character_1 != leaf_character_2) {
            correct <- correct + (frequ[original_character] *
                                    p[5, original_character, internode_character] *
                                    p[1, original_character, leaf_character_1] *
                                    p[2, original_character, leaf_character_1] *
                                    p[3, internode_character, leaf_character_2] *
                                    p[4, internode_character, leaf_character_2])
            wrong1 <-
              wrong1 + (frequ[original_character] *
                          p[5, original_character, internode_character] *
                          p[1, original_character, leaf_character_1] *
                          p[2, original_character, leaf_character_2] *
                          p[3, internode_character, leaf_character_1] *
                          p[4, internode_character, leaf_character_2])
            wrong2 <-
              wrong2 + (frequ[original_character] *
                          p[5, original_character, internode_character] *
                          p[1, original_character, leaf_character_1] *
                          p[2, original_character, leaf_character_2] *
                          p[3, internode_character, leaf_character_2] *
                          p[4, internode_character, leaf_character_1])
          }
        }
      }
    }
  }
  all <- c(correct, wrong1, wrong2)
  return(all)
}