#' predict data-driven only network.
#' 
#' @param data input expression data frame. Gene in coloumns and samples in row. Wih colnames as gene tags and rownames as sample tags.
#' @param cut_off Manually Setup cut off value for the strength of edge. Default at 0.1.
#' @return Parameter of regulization decided lowest BIC.Including regularize parameter(L1 norm parameter) as rho.
#' @examples optimizeParameter(data=GeneExpresssion,prior=Priorindormation)
#' @concept zeptosensPkg
#' @export
predictDatNetwork<-function(data,cut_off=0.1){
    Covmatrix=cov(data)
    
    #optimize penalty parameter rho
    rho <- seq(0.01,1,length=100)
    bic <- rho
    g.result=c()
    p_off_d=c()
    for(i in 1:100){
            g.result  <- glasso(Covmatrix,rho[i])
            p_off_d <- sum(g.result$wi!=0 & col(Covmatrix)<row(Covmatrix))
            bic[i]  <- -2*(g.result$loglik) + p_off_d*log(nrow(data))
    }
    parameter=rho[which.min(bic)]
    
    # Estimated inverse covariance (precision matrix)
    tmp=glasso(Covmatrix,rho = parameter)
    sigma.matrix=tmp$wi
    niter=tmp$niter
    print(niter) # if niter = 10,000
    if(niter == 10000){
        stop("ERROR: Algorithmn does not convergence!")
    }
    pcor.matrix=matrix(0,nrow=ncol(data),ncol=ncol(data))
    for(i in 1:ncol(data)){
        for(j in 1:ncol(data)){
            pcor.matrix[i,j]=-sigma.matrix[i,j]/sqrt(sigma.matrix[i,i]*sigma.matrix[j,j])
        }
    }
    
    t.edges=pcor.matrix
    
    #cut off small edge value edges
    t.net <- as.data.frame(ifelse( abs(t.edges)>=cut_off &row(t.edges)!=col(t.edges),t.edges,0))
    colnames(t.net)=colnames(data)
    rownames(t.net)=colnames(data) 
    
    #sum of edges 
    nedges=sum(t.net!=0)
    
    #Network to edgelist
    edgelist=zeptosensPkg:::createSifFromMatrix(t.net=t.net,genelist = colnames(t.net))
    
    result=list(rho=rho,      nedges=nedges,
                t.net=t.net,  edgelist=edgelist,
                bic=bic)
    return(result)
}
