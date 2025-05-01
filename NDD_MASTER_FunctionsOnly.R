#Package dependencies
#Note:  I'm trying to add more of our own functions so that we can minimize 
#the number of dependencies
#
#
#These libraries help with tree structure and tree graphs
library(data.tree)
library(gtools)
library(igraph)
library(ggraph)
library(DirichletReg)
library(MonoPoly)
library(DiagrammeR)
library(compositions)

#g: igraph object
#r: a root or interior node within g
#function returns the terminal nodes nested under r or the children nodes
get.nodes <- function(g, r,children=F) {
  if(children==F){
  return(names(V(g))[is.finite(distances(g, r, mode = "out")) & degree(g) == 1])}
  if(children==T){return(NA)
  }
    }

#Plotting a tree frame object
#treeframe<- data frame in dataframeNetwork format. Two columns "from" and "to" indicating the root and interior nodes (from)
#mapping to the children nodes beneath
plot.treeframe<-function(treeframe){
  plot(FromDataFrameNetwork(treeframe))
}

#Random Number Generators
#Standard Dirichlet (code ripped straight from MCMC pack)
#Might be helpful to include mean/precision parameteriziation option
rDD<-function(n,alpha){
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x/as.vector(sm))
  
}

rNDD<-function(n,tree.frame){
  
  if(all(unique(tree.frame$from)=="Root")){
    dat<-rDD(n,tree.frame$alpha)
    colnames(dat)<-tree.frame$to
    return(dat)
  }
  
  tree.ig<- graph_from_data_frame( tree.frame )
  leafs<-get.nodes(tree.ig,V(tree.ig)[1])
  layers<-unique(tree.frame$from)
  
  subtrees<-vector(mode = "list", length = length(layers))
  names(subtrees)<-layers
  
  for(i in layers){
    index<-which(tree.frame$from==i)
    dat<-rDD(n,alpha=tree.frame$alpha[index])
    colnames(dat)<-tree.frame$to[index]
    subtrees[[i]]<-dat
  }
  
  leaf.vars<-setdiff(unique(tree.frame$to),unique(tree.frame$from))
  node.vars<-c("Root",setdiff(unique(tree.frame$to),leaf.vars))
  
  dat<-c()
  for(i in leaf.vars){
    index<-which(tree.frame$to==i)
    stop=tree.frame$from[index]
    myprod<-subtrees[[stop]][,colnames(subtrees[[stop]])==i]
    while(!stop=="Root"){
      index<-which(tree.frame$to==stop)
      nm<-stop
      stop=tree.frame$from[index]
      myprod<-myprod*subtrees[[stop]][,colnames(subtrees[[stop]])==nm]
    }
    dat<-cbind(dat,myprod)
  }
  colnames(dat)<-leaf.vars
  return(dat)
  
}

#############################################################################################

#Function to generate subtrees
#x<-compositional data set
#tree.x data frame in dataframeNetwork format. Two columns "from" and "to" indicating the root and interior nodes (from)
#mapping to the children nodes beneath
gen.subtree<-function(x,tree.x){
  
  if(all(unique(tree.x$from)=="Root")){
    subtrees<-list(Root=x)
    return(subtrees)
  }
  
  tree.ig<- graph_from_data_frame( tree.x )
  if(any(x<=0)) return("Values in x are not strictly positive.")
  if(!(any(names(tree.x)==c("from","to")) & ncol(tree.x)==2)) return("tree.x is not in Network format")
  x<-x/apply(x,1,sum)
  leafs<-get.nodes(tree.ig,V(tree.ig)[1])
  layers<-unique(tree.x$from)
  
  #Initializing list
  subtrees<-vector(mode = "list", length = length(layers))
  names(subtrees)<-layers
  
  #Initializing matrix to compute node sums
  node.mat<-c()
  #Computing node sums
  for(i in layers){
    node.mat<-cbind(node.mat,apply(x[,get.nodes(tree.ig,i)],1,sum))
  }
  colnames(node.mat)<-layers
  
  #Combining data with node sums
  x.full<-cbind(x,node.mat)
  
  #Creating list of subcompositions
  for(i in layers){
    child.nodes<-tree.x$to[which(tree.x$from==i)]
    subtrees[[i]]<-x.full[,child.nodes]/x.full[,i]
  }
  
  return(subtrees)
}

#More Quality of life functions
#Function sort the tree data frame from root node down to the lower terminal nodes (top to bottom, left to right)
tree.frame.sort<-function(tree){
  index<-order(tree$from)
  tree<-tree[index,]
  
  index2<-which(tree$from=="Root")
  tree<-rbind(tree[index2,],tree[-index2,])
  return(tree)
}

#Adjusted loglikelihood as defined in the tree finding paper
loglike.adj<-function(alpha,data){
  N<-dim(data)[1]
  logG<-apply(log(data),2,mean)
  return(N*lgamma(sum(alpha))+N*sum((alpha)*logG)-N*sum(lgamma(alpha)) )
}

#Adjusted loglike for a Nested Dirichlet distribution
#Adjusted just means the constant term is removed
#Tree frame with alpha estimates included
loglike.tree<-function(tree,data){
  N<-dim(data)[1]
  mytree<-tree.frame.sort(tree)
  mysub<-gen.subtree(data,mytree[,c("from","to")])
  layers<-names(mysub)
  logLike<-0
  for (i in layers){
    index<-which(mytree$from==i)
    logLike<-logLike+loglike.adj(alpha=mytree$alpha[index],data=mysub[[i]])
  }
  return(logLike)
}

#################################################################################################
#Function to estimate parameters of an NDD when the tree is given
#

#Given a composition data set and tree frame, fit the parameters. of the NDD
#Default method produces parameter estimates under the alpha parameters parms="alpha", means and precisions for each layer are provided if parms="mean"
#Default output="dataframe" appends estimate information to the tree data frame provided by user.
#If output="lists", then a list is produced with corresponding parameter estimates for each layer.
#Function always returns a list as the final loglikelihood value is provided.

tree.fit<-function(dataset,tree,parms="alphas",output="dataframe"){
  mytree<-tree.frame.sort(tree)
  mysub<-gen.subtree(dataset,mytree)
  mysub
  layers<-names(mysub)
  
  if(output=="dataframe"){
    alphas<-c()
    mu<-c()
    prec<-c()
    logLike<-0
    for(i in layers){
      newd<-DR_data(mysub[[i]])
      mod1<-DirichletReg::DirichReg(newd ~ 1 , model = "common")
      a<-exp(mod1$coefficients)
      names(a)<-colnames(newd)
      logLike<-logLike+loglike.adj(a,newd)
      alphas<-c(alphas,a)
      mu<-c(mu,a/sum(a))
      prec<-c(prec,rep(sum(a),length(a)))
    }
    if(parms=="alphas"){mytree$alpha<-alphas}
    if(parms=="mean"){
      mytree$mean<-mu
      mytree$precision<-prec
    }
    
    return(list(fit=mytree,logLike=logLike))
  }
  
  
  if(output=="list"){
    results<-list()
    logLike<-0
    for(i in layers){
      newd<-DR_data(mysub[[i]])
      mod1<-DirichletReg::DirichReg(newd ~ 1 , model = "common")
      alphas<-exp(mod1$coefficients)
      names(alphas)<-colnames(newd)
      logLike<-logLike+loglike.adj(alphas,newd)
      if(parms=="mean"){
        results[[i]]<-list(mean=alphas/sum(alphas),precision=sum(alphas))
      }
      if(parms=="alphas"){
        results[[i]]<-list(alphas=alphas)
      }
    }
    results$logLike<-logLike
    return(results)  
  }
}

####################################################################################################
#The following functions are used INSIDE the tree finding algorithm
#1. MLE estimation of standard dirichlet (we could also used DirichletReg package to do this like we do in other parts)
#2. Force.comp forces the compositional structure on the data set (X1,X2,...Xk)/sum(X's)
#3. split.tree function is used to determine if an internal node should be made for a given set of variables.

dirichletmlealphas<-function(X, ind, N)
{
  # Method: Using the Newton-Raphson algorithm (2.69)
  # Input: X: an m x n observed data matrix
  # ind=1: using (2.71) as the initial values
  # ind=2: using (2.73) as the initial values
  # N: the number of iterations required for
  # the Newton-Raphson algorithm
  m <- dim(X)[1]
  n <- dim(X)[2]
  one <- rep(1, n)
  xmean <- apply(X, 2, mean)
  G <- (apply(X, 2, prod))^(1/m)
  logG<-apply(log(X),2,sum)/m
  if(ind == 1) {
    a <- min(X) * one
  }
  if(ind == 2) {
    ga <- - digamma(1)
    b <- xmean
    de0 <- sum(b * log(G))
    aplus <- ((n - 1) * ga)/(sum(b * log(b)) - de0)
    a <- b * aplus
  }
  A <- matrix(0, N, n)
  for(tt in 1:N) {
    aplus <- sum(a)
    #g <- m * (digamma(aplus) - digamma(a) + log(G))
    g <- m * (digamma(aplus) - digamma(a) + logG)
    b <- m * trigamma(aplus)
    Binv <- - diag(1/trigamma(a))/m
    Hinv <- Binv - (Binv %*% one) %*% (t(one) %*%
                                         Binv)/(1/b + c(t(one) %*% Binv %*% one))
    a <- a - c(Hinv %*% g)
    A[tt, ] <- a
  }
  return(a)
}

force.comp<-function(data,index=1:ncol(data)){
  #Function that forces a composition in 1 on two ways
  #By default takes the entire matrix and divides by rowsums
  #If a subset of the columns are provided, it produces a 2 variable composition
  #where one variable is the sum of the variable index and the second variable is the 
  #sum of the remaining variables "-index"
  index<-index[index>0]
  if(length(index)==ncol(data)){newdata<-data/apply(data,1,sum)}
  if(length(index)<ncol(data)){newdata<-cbind(apply(data[,index,drop=F],1,sum),apply(data[,-index,drop=F],1,sum))}   
  return(newdata)
}

split.tree<-function(index,dataset,method="MaxLik"){
  #For a given subset of variables listed in "index" contained in "dataset",
  #transforms the subset into a composition, fits a standard dirichlet,
  #then computes all possible binary splits along with fit criterion for the NDD,
  #returns a matrix with two rows indicating the variables contained in the split
  #If the standard dirichlet is preferred, returns matrix of 0's
  dataset.dim<-dim(dataset)[2]
  index.vec<-index[!index==0]
  indicator.vec<-1:length(index.vec)
  
  #No splitting possible returning 0's
  if(length(index.vec)<3){return(matrix(rep(0,2*length(index)),nrow=2) )}
  
  #Splitting possible, finding best possible binary split
  if(length(index.vec)>2){
    datasettry<-force.comp(dataset[,index.vec]) #Normalizing data to composition
    n<-dim(datasettry)[[1]] #sample size
    n_var<-dim(datasettry)[[2]] #number of variables in composition
    j<-floor(n_var/2) #maximum number of variables to consider in a single split
    
    #Creating all possible combinations to create a binary split
    #via variable index
    labels<-matrix(rep(0,j),1,j) 
    for(l in 1:j){
      combos<-cbind(combinations(n_var,l),matrix(rep(0,(j-l)*dim(combinations(n_var,l))[1]),nrow=dim(combinations(n_var,l))[1]))
      if(l*2==n_var){combos<-combos[1:(nrow(combos)/2),]}
      labels<-rbind(labels,combos)
    }
    
    #Fitting standard dirichlet and computing fit criterion
    Criterion<-c()
    top.fitsDD<-dirichletmlealphas(datasettry,ind=1,N=40)
    top.fits<-c(0,0)
    Criterion[1]<- -2*loglike.adj(top.fitsDD,datasettry)
    
    #Evaluating all possible fits.  Storing the best while iterating through
    for(i in 2:nrow(labels)){
      alpha1<-dirichletmlealphas(force.comp(datasettry,index=labels[i,]),ind=1,N=40)
      alpha2<-dirichletmlealphas(force.comp(datasettry[,-labels[i,],drop=F]),ind=1,N=40)
      if(i<(n_var+2)){
        Criterion[i]<- -2*(loglike.adj(alpha1,force.comp(datasettry,index=labels[i,]))+loglike.adj(alpha2,force.comp(datasettry[,-labels[i,],drop=F])))  
      }
      else{
        alpha3<-dirichletmlealphas(force.comp(datasettry[,labels[i,],drop=F]),ind=1,N=40)
        Criterion[i]<- -2*(loglike.adj(alpha1,force.comp(datasettry,index=labels[i,]))+loglike.adj(alpha2,force.comp(datasettry[,-labels[i,],drop=F]))+loglike.adj(alpha3,force.comp(datasettry[,labels[i,],drop=F])))
      }
      top.fits<-rbind(top.fits,alpha1)
    }
  }
  k=c(n_var,rep(n_var+1,n_var),rep(n_var+2,nrow(labels)-n_var-1))
  AIC<- 2*k+Criterion
  AICc<- AIC+2*k*(k+1)/(n-k-1)
  BIC<- -2*Criterion+k*log(n)
  
  
  if(method=="MaxLik"){   loc<-which(Criterion==min(Criterion))[1]}
  if(method=="AIC"){   loc<-which(AIC==min(AIC))[1]}
  if(method=="AICc"){   loc<-which(AICc==min(AICc))[1]}
  if(method=="BIC"){   loc<-which(BIC==min(BIC))[1]}
  
  com.vec<-1:n_var
  subset1<-labels[loc,]
  subset2<-com.vec[-which(com.vec%in%subset1)]
  subset1.final<-index.vec[subset1]
  subset2.final<-index.vec[subset2]
  
  mymat<-rbind(c(subset1.final,rep(0,dataset.dim-length(subset1.final))), c(subset2.final,rep(0,dataset.dim-length(subset2.final))))
  #print(cbind(labels,-2*Criterion,AIC,AICc,BIC))
  
  return(mymat)
  
}

##################################################################################################################
#Tree finding algorithm
#Note: You must run the functions above first

#dataset must be compositional (need to add check and return error if not)
#method can be MaxLik, AIC,AICc,or BIC, options are fed to split.tree()

#Note: Currently the function plots the final resulting tree and also provides the tree.frame
#structure.  This means that you will have to fit the data set one more time with the tree
#that was determined to obtain the MLEs.  TODO:  The function should return the alpha estimates along
#with the tree frame.
LDM.alg<-function(dataset,meth="MaxLik"){
  
  data<-as.matrix(dataset)
  #direction.mat<-c()
  splits.mat<-c()
  splits.temp<-c()
  #direction.temp<-c()
  #dummy<-c()
  #j=2
  first.index<-1:dim(data)[2]
  
  
  splits.mat<-rbind(splits.mat,split.tree(first.index,data,method=meth))
  splits.temp<-splits.mat
  tree.frame<-data.frame(to=as.character(first.index),from=rep("Root",dim(data)[2])) 
  total<-sum(splits.mat)
  
  if(total==0){
    tree.fit<- FromDataFrameNetwork(tree.frame)
    x=plot(tree.fit) 
    x
    return(list(x,tree.frame[,2:1]))
  }
  
  while(total!=0){
    dummy<-c()
    for (i in 1:dim(splits.temp)[1]){
      dummy<-rbind(dummy,split.tree(splits.temp[i,],data,method=meth))
    }
    total<-sum(dummy)
    splits.temp<-dummy
    splits.mat<-rbind(splits.mat,splits.temp)
  }     
  
  splits.mat<-splits.mat[rowSums(splits.mat)!=0,]
  if(nrow(splits.mat)==0){
    tree.fit<- FromDataFrameNetwork(tree.frame)
    plot(tree.fit)
  }
  else{
    cnt=1
    for(j in 1:nrow(splits.mat)){
      subs<-splits.mat[j,]
      if(sum(subs>0)<2) {next}
      else{
        subs<-as.character(subs[subs>0])
        ind<-which(tree.frame$to %in% subs)
        nested.node<-tree.frame$from[ind[1]]
        tree.frame$from[ind]<-paste("N",cnt,sep="")
        tree.frame<-rbind(tree.frame,c(paste("N",cnt,sep=""),nested.node))
        cnt<-cnt+1
      }
    }
    tree.frame$to[first.index]<-colnames(data)
    tree.fit<- FromDataFrameNetwork(tree.frame)
    x=plot(tree.fit) 
    x
    #Playing around with other plots
    #plot(as.dendrogram(tree.fit),center=T)
    #plot(as.igraph(tree.fit, directed = TRUE))
    
    #useRtreeList <- ToListExplicit(tree.fit, unname = TRUE)
    #radialNetwork( useRtreeList)
  }
  return(list(x,tree.frame[,2:1]))
}
#Work In progress
################################################################
#Function to conduct LRT test statistic and p-value for 2 group
#comparison introduced in Turner et al. 2024.  
#
#Function inputs:
#   x,Y  two compositional data sets or X is a 

####################################################################
#LRT Test for the NDD model applied to Water Maze data
#
#


#Transforming to 3 sub-tree compositions based on the tree given 
#in Figure 6
#

LDM.test<-function(dataset,group,tree){
  mytree<-tree.frame.sort(tree)
  mysub<-gen.subtree(dataset,mytree)
  layers<-names(mysub)
  p<-ncol(dataset)
  LRT.vec<-c()   #LRT vector for each layer  
  for (i in layers){
    dat<-data.frame(mysub[[i]],Group=group)
    subp<-ncol(dat)-1
    dat$Y<-DR_data(dat[,1:subp])
    mod_alt <- DirichletReg::DirichReg(Y ~ Group | Group,data=dat, model="alternative")
    mod_null <- DirichletReg::DirichReg(Y ~ 1 | Group,data=dat, model = "alternative")
    LRT.vec<-c(LRT.vec, 2*(mod_alt$logLik-mod_null$logLik))
    
    }
  result1<-data.frame(Layer=layers,LRT.layer=LRT.vec)
  LRTtot<-sum(LRT.vec)
  result2<-data.frame(LRT.stat=LRTtot,p.value=pchisq(LRTtot,df=p-1,lower.tail=F))
  final<-list(result1,result2)
  names(final)<-c("LRT.Layers","Overall")
  return(final)
}
  

###############################################################################################################################
#These functions are used to produce pseudo residuals and marginal parameters for investigating an NDD model fit.
#At the end of this section the functions are utilized to make 1.qqplots,2. histogram with marginal fits,
# 3. Influential diagnostic, 4.  Aitchisons distance of each observation from the estimated mean
# For a package, it might be helpful to utilize these codes within the tree.fit function, so that everything is obtained in one 
#function call.  Then maybe a plot.ndd() function could be used to produce the plots, OR at a minimum, the objects would be available
#So information is easily extracted from lists for the user to do whatever they way.


######
#This function takes the paramater matrix from tree.fit and produces the appropriate Beta paramaters that completely parameterize the marginal's of each 
#component in the fitted NDD.
#Input:  tree.fit object
#Output:  list of matrices (1 for each variable in the composition). Matrices all have 2 cols (one for alpha and one for beta parameters).  The # rows depends on the component and is equal to howw many 
#beta distributions are used in the product to create the marginal.  This is equal to the number of branches starting at the root node leading down the path to a terminal node (component).
#See tree finding algorithm paper for technical details

beta.parms<-function(NDDfit){
  if(!is.data.frame(NDDfit)){NDDfit<-NDDfit$fit}
  leaf.vars<-setdiff(unique(NDDfit$to),unique(NDDfit$from))
  result<-list()
  node.vars<-c("Root",setdiff(unique(NDDfit$to),leaf.vars))
  
  #if tree fit provides alpha values must convert to mean/precision
  if(any(names(NDDfit)=="alpha")){
    NDDfit$mean<-0
    NDDfit$precision<-0
    for (i in unique(NDDfit$from)){
      index<-which(NDDfit$from==i)
      A<-sum(NDDfit$alpha[index])
      NDDfit$precision[index]<-A
    }
    NDDfit$mean<-NDDfit$alpha/NDDfit$precision
    NDDfit<-NDDfit[,c("from","to","mean","precision")]
  }
  
  for(i in leaf.vars){
    index<-which(NDDfit$to==i)
    stop=NDDfit$from[index]
    path<-index
    while(!stop=="Root"){
      index<-which(NDDfit$to==stop)
      stop=NDDfit$from[index]
      path<-c(path,index)
    }
    x<-NDDfit[path,3:4]
    x$alpha=x$mean*x$precision
    x$beta=x$precision-x$alpha
    result[[i]]<-x[,3:4]
  }
  
  return(result)
}

###################################################################################################################
#The following are all functions used to approximate the marginals
#Some of these functions should probably just be internal function (not for user use)
#The ones people will want to use are saddle.pdf and saddle.cdf, details on the function inputs
#are stated below in the first function.  The same options are consistent throughout the rest of the functions in this section. 
#saddle.pdf will be used for grahping density fits, saddle.cdf will be used for computing pseudoresiduals of an NDD fit.

#Cumulant generating function of the product of l indepdendent beta distributions Y=X_1*X_2*...*X_l.  
#The definition of k(s) is stated in the tree finding paper
#the inputs alpha and beta correspond to the l alpha and l beta parameters for the l Beta random variables.
#The int option will compute up to 3rd derivatives of k(s).  int=0 corresponds to k(s), int=1 corresponds to k'(s),
#int=2 corresponds to k''(s), and int=3 corresponds to k'''(s).
k.s<-function(s,alpha,beta,int=0){
  result<-c()
  if(int==0){
    for(i in 1:length(s)){
      result[i]<-sum(lgamma(alpha+s[i])+lgamma(alpha+beta)-lgamma(alpha)-lgamma(alpha+beta+s[i]))
    }
    return(result)
  }
  
  if(int==1){
    for(i in 1:length(s)){
      result[i]<-sum(digamma(alpha+s[i])-digamma(alpha+beta+s[i]))
    }
    return(result)
  }
  
  if(int==2){
    for(i in 1:length(s)){
      result[i]<-sum(trigamma(alpha+s[i])-trigamma(alpha+beta+s[i]))
    }
    return(result)
  }
  
  if(int==3){
    for(i in 1:length(s)){
      result[i]<- sum(psigamma(alpha+s,deriv=2)-psigamma(alpha+beta+s,deriv=2))
    }
    return(result)
  }
}


#Function to solve the saddle point equation k'(shat)=y  when approximating Y=X_1*X_2*...*X_l
#Function resutns a single value, shat.
s.hat<-function(x,alpha,beta){
  answer<-uniroot(function(s,a,b,y){k.s(s,a,b,int=1)-y},a=alpha,b=beta,y=x,interval=c(-min(alpha)+.00001,10000),extendInt="upX",tol=.Machine$double.eps^0.35)
  return(answer$root)
}


#Computes unnormalized saddlepoint pdf of Y=X_1*X_2*...*X_l
saddle.unnorm<-function(x,alpha,beta){
  result<-c()
  for(i in 1:length(x)){
    mys<-s.hat(x[i],alpha,beta)
    result[i]<-(2*pi*k.s(mys,alpha,beta,int=2))^(-.5)*exp(k.s(mys,alpha,beta,int=0)-mys*x[i])
  }
  return(result)
}


##Computes normalized saddlepoint pdf of Y=X_1*X_2*...*X_l
# The normalization occurs by numerically integrating saddle.unnorm to obtain a normalizing constant.  Function returns saddle.nnorm/normalizing constant so that it 
#is a proper density.  
#This function operates just like dnorm, dt, dexp and other density functions in R.  We could rename these functions to keep the same convention pbetaprod or pNDDmarg

saddle.pdf<-function(x,alpha,beta){
  result<-c()
  y<-x
  x<-log(x)
  constant<-integrate(saddle.unnorm,lower=-Inf,upper=-.0001,alpha=alpha,beta=beta)$value
  for(i in 1:length(x)){
    if(x[i]>=0){result[i]=0}
    else{
      mys<-s.hat(x[i],alpha,beta)
      result[i]<-(1/y[i])*constant^(-1)*(2*pi*k.s(mys,alpha,beta,int=2))^(-.5)*exp(k.s(mys,alpha,beta,int=0)-mys*x[i]) 
    }
  }
  return(result)
}


saddle.cdf.internal<-function(x,alpha,beta){
  result<-c()
  for(i in 1:length(x)){
    mys<-s.hat(x[i],alpha,beta)
    if(x[i]!=k.s(0,alpha,beta,int=1)){
      w.hat<-sign(mys)*sqrt(2*(mys*x[i]-k.s(mys,alpha,beta,int=0)))
      u.hat<-mys*sqrt(k.s(mys,alpha,beta,int=2))
    }
    if(x[i]==k.s(0,alpha,beta,int=1)) {result[i]<-0.5+k.s(0,alpha,beta,int=3)/6/sqrt(2*pi)/(k.s(0,alpha,beta,int=2))^(1.5)}
    else {result[i]<-pnorm(w.hat)+dnorm(w.hat)*(1/w.hat-1/u.hat)}
  }
  return(result)
}


#Computes saddlepoint cdf of Y=X_1*X_2*...*X_l, works like pnorm, pbeta, pexp, pgamma, etc.
saddle.cdf<-function(x,alpha,beta){
  probs<-c()
  x<-log(x)
  probs[x>=0]<-1
  if(all(!is.na(probs))) {return(probs)}
  
  true.mean<-k.s(0,alpha,beta,int=1)
  true.var<-k.s(0,alpha,beta,int=2)
  L<-true.mean-.1*sqrt(true.var)
  U<-true.mean+.1*sqrt(true.var)
  
  good.index<- (x<L) | (x>U & x<0)
  interp.index<- (x<=U) & (x>=L)
  if(any(good.index)){
    probs[good.index]<-saddle.cdf.internal(x[good.index],alpha,beta)
  }
  
  #only performing interpolation if necessary
  if(any(interp.index)){
    from=true.mean-1*sqrt(true.var)
    to=min(true.mean+1*sqrt(true.var),-.0001)
    myseq<-seq(from,to,length.out=1000)
    myseq<-myseq[!(myseq>=L & myseq<=U)]
    myseq<-c(myseq,true.mean)
    fhats<-saddle.cdf.internal(myseq,alpha,beta)
    mod<-monpol(fhats~myseq,degree=9,a=from,b=to)
    probs[interp.index]<-predict(mod,newdata=data.frame(myseq=x[interp.index]))  
  }
  
  return(probs)
  
}
########
#Function to compute residuals and diagnostics
#This function is designed to be the go to user function rather than tree.fit.  
#This function will basically operate just like tree.fit, but now additional results will be stored and returned

#The example after this function will be two full blown analysis of data sets as well as graphics that can be added to the paper
#In terms of an R package, we should probably create a unique plot function for a fitted result.  Probably only appropriate for like
#up to 9 components (3x3) graphics.
NDD.fit<-function(dataset,tree,residuals=F,diagnostics=F){
  mytree<-tree.frame.sort(tree)
  mysub<-gen.subtree(dataset,mytree)
  mysub
  layers<-names(mysub)
  n=nrow(dataset)
  p=ncol(dataset)
 
    alphas<-c()
    mu<-c()
    prec<-c()
    logLike<-0
    for(i in layers){
      newd<-DR_data(mysub[[i]])
      mod1<-DirichletReg::DirichReg(newd ~ 1 , model = "common")
      a<-exp(mod1$coefficients)
      names(a)<-colnames(newd)
      logLike<-logLike+loglike.adj(a,newd)
      alphas<-c(alphas,a)
      mu<-c(mu,a/sum(a))
      prec<-c(prec,rep(sum(a),length(a)))
    }
      mytree$alpha<-alphas
      mytree$mean<-mu
      mytree$precision<-prec
    results<-list(fit=mytree,logLike=logLike)
    marginfo<-beta.parms(mytree)
    #return(results)
    #Diagnostics
    if (residuals==T){
    
    cnames<-names(marginfo)
    res<-c()
    results$parmsmarg<-marginfo
    #return(results)
    for(i in 1:p){
      cname<-cnames[i]
      res<-cbind(res,qnorm(saddle.cdf(dataset[,cname],alpha=marginfo[[cname]]$alpha,beta=marginfo[[cname]]$beta)))
    }
    res<-data.frame(res)
    names(res)<-cnames
  
    results$resids<-res
    }
    
    if(diagnostics==T){
      d<-c()
      ll<-c()
      lsample<-log(dataset)
      lmean<-apply(lsample,1,mean)
      pihat<-unlist(lapply(marginfo,function(x){prod(x$alpha/(x$alpha+x$beta))}))
      constant<-log(pihat)/mean(log(pihat))
      d<-apply((t(lsample-lmean)-constant)^2,2,sum)
      for (i in 1:n){
        #ll[i]=tree.fit(dataset[-i,],tree=tree)$logLike
        ll[i]=loglike.tree(tree.fit(dataset[-i,],tree=tree)$fit,data=dataset)
        }
      results$D<-d
      results$LD<- 2*(logLike-ll)  
    }
    
    
    return(results)  
  
}

################################
#A summary function for the NDD.
#input is a tree frame with from,to, and alpha colums to specify the distribution
#outputs the mean, sd, and covariance matrix
#Note: use cov2cor() to obtain correlation matrix
#Note: summary can be used for theoretical calculations OR for sample estimates 
#      by feeding in the treeframe object returned by NDD.fit
summary.tree<-function(tree.frame){
  if(!all(c("from","to","alpha") %in% names(tree.frame))){stop("Tree frame is not complete.  Must have 3 colums: 'from','to',and 'alpha' ")}
  tree.frame<-tree.frame.sort(tree.frame)
  tree.ig<- graph_from_data_frame( tree.frame )
  leafs<-get.nodes(tree.ig,V(tree.ig)[1])
  layers<-unique(tree.frame$from)
  p<-length(leafs)
  
  #obtain 1st and 2nd moment and precision of branches
  tree.frame$x<-0
  tree.frame$phi<-0
  tree.frame$x2<-0
  layer.phis<-c()
  for (i in layers){
    index<-which(tree.frame$from==i)
    a<-tree.frame$alpha[index]
    phi<-sum(a)
    tree.frame$x[index]<-a/phi
    tree.frame$phi[index]<-phi
    tree.frame$x2[index]<-(a*(a+1))/(phi*(phi+1))
    layer.phis<-c(layer.phis,phi)
  }
  names(layer.phis)<-layers
  
  
  #computing mean vector
  #Adding branch means to tree object
  mypaths<-shortest_paths(tree.ig,from="Root",to=get.nodes(tree.ig,V(tree.ig)[1]))$vpath
  E(tree.ig)$weight<-tree.frame$x
  names(mypaths)<-leafs
  
  #computing mean by taking product of branches along the paths from root to terminal variable
  mean.mat<-rep(0,p)
  for (i in 1:length(mypaths)){
    mean.mat[i]=prod(E(tree.ig,path=mypaths[[i]])$weight)
  }
  
  #computing variances
  var.mat<-rep(0,p)
  E(tree.ig)$weight<-tree.frame$x2
  for (i in 1:length(mypaths)){
    var.mat[i]=prod(E(tree.ig,path=mypaths[[i]])$weight)-mean.mat[i]^2
  }
  
  
  #Computing E(X_i,X,j) for COV computations
  mypaths.nam<-lapply(mypaths,as_ids)
  
  #Initializing
  cov.mat<-matrix(rep(0,p*p),p,p)
  for(i in 1:p){
    for(j in i:p){
      if(i==j){cov.mat[i,j]<-var.mat[i]}
      else{
        com.ver<-intersect(mypaths.nam[[i]],mypaths.nam[[j]])
        l=length(com.ver)
        final.phi<-layer.phis[com.ver[l]]
        const<-final.phi/(final.phi+1)
        if( l>1){
          index<-which(tree.frame$to %in% com.ver)
          alphas<-tree.frame$alpha[index]
          phis<-tree.frame$phi[index]
          const<-const*prod((alphas+1)*phis/alphas/(phis+1))
        }
        cov.mat[i,j]<-(const-1)*mean.mat[i]*mean.mat[j]
        cov.mat[j,i]<-cov.mat[i,j]
      }
    }
  }
  mu=mean.mat
  names(mu)<-names(mypaths)
  stdev=sqrt(var.mat)
  names(stdev)<-names(mypaths)
  covar=cov.mat
  colnames(covar)<-names(mypaths)
  rownames(covar)<-names(mypaths)
  results<-list(mu,stdev,covar)
  names(results)<-c("mu","stdev","covar")
  return(results)
  
  
}
