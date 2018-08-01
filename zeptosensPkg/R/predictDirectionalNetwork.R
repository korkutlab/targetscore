#' Predicted Directional GeneNetwork
#' 
#' @param data input expression data. Coloumns as the gene, rows as the sample.With colnames as the gene tags, rownames as the sample tags.
#' @param prior prior information matrix in form of p*p gene interaction matrix.Colnames and rownames as the gene tags.Coloumns and rows correspond to the data matrix gene tags.
#' @param rho regulization parameter
#' @param kappa scaler parameter
#' @return estimated.network include list of estimated directional partial correlation gene network rho as the regulization parametwe and kappa as the scaler parameter.
#' @concept zeptosensPkg
#' @export

predictDirectionalNetwork=function(data,prior,rho,kappa,cut.off){
    U=matrix(1,ncol(data),ncol(data))
    rho_m=rho*U-kappa*prior
    pc=cov(data) 
    #Network construction with directional prior information
    sigma.matrix<-glasso(pc,rho=rho_m)$wi
    pcor.matrix=matrix(0,nrow=ncol(data),ncol=ncol(data))
    for(i in 1:ncol(data)){
        for(j in 1:ncol(data)){
            pcor.matrix[i,j]=-sigma.matrix[i,j]/sqrt(sigma.matrix[i,i]*sigma.matrix[j,j])
        }
    }
    
    t.edges.rhoadjusted=pcor.matrix
    
    #get direction for the network as the bigger covariance estimated indicated the upper stream gene;
    t.edges.rhoadjusted.d<-matrix(0,nrow = nrow(t.edges.rhoadjusted),ncol = ncol(t.edges.rhoadjusted))
    for(j in 1:ncol(t.edges.rhoadjusted)){
        for(i in 1:nrow(t.edges.rhoadjusted)){
            if(abs(t.edges.rhoadjusted[i,j])>abs(t.edges.rhoadjusted[j,i])){
                t.edges.rhoadjusted.d[i,j]=t.edges.rhoadjusted[i,j]
            }
            if(abs(t.edges.rhoadjusted[i,j])<abs(t.edges.rhoadjusted[j,i])){
                t.edges.rhoadjusted.d[j,i]=t.edges.rhoadjusted[j,i]
            }
            if(abs(t.edges.rhoadjusted[i,j])==abs(t.edges.rhoadjusted[j,i])){
                t.edges.rhoadjusted.d[i,j]=t.edges.rhoadjusted[i,j]
                t.edges.rhoadjusted.d[j,i]=t.edges.rhoadjusted[j,i]
            }
        }
    }
    #cutoff at cutoff point
    t.edges.rhoadjusted.d=ifelse(abs(t.edges.rhoadjusted.d)>cut.off,t.edges.rhoadjusted.d,0)
    estimated.network=list(t.edges.rhoadjusted.d,rho,kappa)
    return(estimated.network)
}