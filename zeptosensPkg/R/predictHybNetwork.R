#' Choose the optimal regulization parameter and scale paramter for prior information adjusted network construction.
#' 
#' @param data input expression data frame. Gene in coloumns and samples in row. Wih colnames as gene tags and rownames as sample tags.
#' @param prior Prior information data frame ,with colnames and rownames as gene tags.With colnames and rownames as gene tags. Can be inferred from predictBioNetwork or any network resources.
#' @param cut_off Manually set up cut off value for strength of edge. Default at 0.1.
#' @param maxDist Maximun distance for the network. Default at 1.
#' @param proteomicResponses RPPA data tested for drug pertubation.
#' @param nProt number of Proteins contained in the data.
#' @return "parameters" as the Parameter list of regulization parameter decided by the prior information and the algorithmn lowest BIC.Including regularize parameter(L1 norm parameter) as "rho", scale parameter(decided how much prior information contribute) as "kappa", and regulization matrix for the expression data as "rho_m". 
#' @return "bic"as the Model's BIC matrix for differnet regularization parameters.
#' @return "wk" as the predicted network. 
#' @return "wks" TBA
#' @return "dist_ind" TBA
#' @return "inter" TBA
#' @return "edgelist" as the edgelist for predicted network.
#' @return "nedges" as the number of edges of the predicted network.
#' @examples predictHybNetwork(data=GeneExpresssion,prior=Priorinformation,proteomicResponses=proteomicResponses,nProt=nProt)
#' @concept zeptosensPkg
#' @export
predictHybNetwork<-function(data,prior=NULL, cut_off=0.1,proteomicResponses,nProt,maxDist=1,antibodyMapFile){
  if(is.null(prior)){
    wk=zeptosensPkg:::predictBioNetwork(nProt = 304 ,proteomicResponses = proteomicResponses,maxDist = 1,antibodyMapFile =antibodyMapFile)$wk
    prior=wk
  }  
  
    #HybNetwork
    
    index=colnames(prior[,which(colnames(prior)%in%colnames(data))])#match the data
    data=data[,index]
    prior1=prior[index,index]
    
    #prior information extration
    prior1=ifelse(prior1!=0,1,0)#information matrix of prior
    prior2=prior1               #symmetrical prior information
    for (i in 1:nrow(prior1)){
        for (j in 1:ncol(prior1)){
            if(prior1[i,j]!=0){
                prior2[i,j]=prior1[i,j]
                prior2[j,i]=prior1[i,j]
            }
        }
    }
    prior2=ifelse(prior2!=0,1,0)
    
    # getting the best tuning parameter from BIC minimization
    Covmatrix=cov(data)
    rho <- seq(0.01,1,length=100)
    bic <- matrix(NA,100,100)
    kappa<-rho
    rho_m=c()
    g.result=c()
    U=matrix(1,nrow(prior2),ncol(prior2))
    p_off_d=c()
  for(i in 1:100){
    for(j in 1:i){
        rho_m=rho[i]*U-kappa[j]*prior2
        g.result  <- glasso(Covmatrix,rho_m)
        p_off_d <- sum(g.result$wi!=0 & col(Covmatrix)<row(Covmatrix))
        bic[i,j]  <- -2*(g.result$loglik) + p_off_d*log(nrow(data))
        bic=as.data.frame(bic)
        rownames(bic)=rho
        colnames(bic)=kappa
    }
   }
    pos=which(bic==min(bic,na.rm = TRUE),arr.ind=T)
    rho=rho[pos[1]]
    kappa=kappa[pos[2]]
    rho_m=rho*U-kappa*prior2
    parameters=list(rho_m,rho,kappa)
   
     # Estimated inverse covariance (precision)
    tmp=glasso(Covmatrix,rho = rho_m)
    sigma.matrix=tmp$wi
    niter=tmp$niter
    print(niter)# if niter = 10,000
    if(niter==10000){
      stop("ERROR: Algorithmn does not convergence!")
    }
    
    pcor.matrix=matrix(0,nrow=ncol(data),ncol=ncol(data))
    for(i in 1:ncol(data)){
      for(j in 1:ncol(data)){
        pcor.matrix[i,j]=-sigma.matrix[i,j]/sqrt(sigma.matrix[i,i]*sigma.matrix[j,j])
      }
    }
    
    t.edges=pcor.matrix
    
    #cut off = cut off value
    t.net <- as.data.frame(ifelse( abs(t.edges)>=cut_off &row(t.edges)!=col(t.edges),t.edges,0))
    colnames(t.net)=colnames(data)
    rownames(t.net)=colnames(data) 
    
    #directed network 
    #get direction for the network as the bigger covariance estimated indicated the upper stream gene;
    t.net.d<-matrix(0,nrow = nrow(t.net),ncol = ncol(t.net))
    for(i in 1:ncol(t.net)){
      for(j in 1:nrow(t.net)){
        if(prior1[i,j]!=0&prior1[j,i]==0){
          t.net.d[i,j]=t.net.d[i,j]
        }
        if(prior1[i,j]==0&prior1[j,i]==0){
          t.net.d[i,j]=t.net[i,j]
          t.net.d[j,i]=t.net[j,i]
        }
        if(prior1[i,j]!=0&prior1[j,i]!=0){
          t.net.d[i,j]=t.net[i,j]
          t.net.d[j,i]=t.net[j,i]
        }
      }}
    colnames(t.net.d)=colnames(data)
    rownames(t.net.d)=colnames(data)
    
    #Inferred missing proteins from bio-network
    network=t.net.d
    index2=which(colnames(prior)%in%index)
    prior1.extra=prior[-index2,-index2]
    index.extra=colnames(prior1.extra)
    
    #with the edgevalue(prior information have the value of the median(abd(netowrk)))
    prior3=as.data.frame(prior*median(abs(c(na.omit(c(ifelse(network!=0,network,NA)))))))
    colnames(prior3)=colnames(prior)
    rownames(prior3)=rownames(prior)
    
    part1=prior3[which(rownames(prior3)%in%index.extra),which(colnames(prior)%in%index)]
    part2=prior3[which(rownames(prior3)%in%index),colnames(prior)%in%index.extra]
    part3=prior3[which(rownames(prior3)%in%index.extra),which(colnames(prior)%in%index.extra)]
    
    network.total=cbind(rbind(network,part1),rbind(part2,part3))
    index3=match(colnames(prior),colnames(network.total))
    network.total=network.total[index3,index3]
    
    edgelist.total=zeptosensPkg:::createSifFromMatrix(t.net = network.total,genelist = colnames(network.total))
    
    #number of edges
    nedges=sum(network.total!=0)
     
    #
    wk=network.total
    networks <- zeptosensPkg:::network2(wk=wk,nProt=nProt, 
                         proteomicResponses=proteomicResponses,
                         maxDist=maxDist)
    wk <- networks$wk
    wks <- networks$wks
    dist_ind <- networks$dist_ind
    inter <- networks$inter
    
    #return result
    result=list(parameters=parameters,nedges=nedges, inter=inter,
                wk=wk,wks=wks, dist_ind=dist_ind, edgelist=edgelist.total,
                bic=bic)
    return(result)
}