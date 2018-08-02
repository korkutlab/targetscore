#' Choose the optimal regulization parameter and scale paramter.
#' 
#' @param data input expression data. Gene in coloumns and samples in row. Wih colnames as gene tags and rownames as sample tags.
#' @param prior Prior information matrix,with colnames and rownames as gene tags.
#' @return Parameter list of regulization parameter decided by the prior information and the algorithmn lowest BIC.Including regularize parameter as rho, scale parameter as kappa, and regulization matrix decided.
#' @concept zeptosensPkg
#' @export
optimizeParameter<-function(data,prior){
    index=colnames(prior[,which(colnames(prior)%in%colnames(data))])#match the data

    data=data[,index]
    prior=prior[index,index]
    prior=ifelse(prior!=0,1,0)#information matrix of prior
    
         Covmatrix=cov(data)
         rho <- seq(0.01,1,length=100)
         bic <- matrix(NA,100,100)
         kappa<-rho
         rho_m=c()
        g.result=c()
         U=matrix(1,nrow(prior),ncol(prior))
         p_off_d=c()
for(i in 1:100){
    for(j in 1:i){
        rho_m=rho[i]*U-kappa[j]*prior
        g.result  <- glasso(Covmatrix,rho_m)
        p_off_d <- sum(g.result$wi!=0 & col(Covmatrix)<row(Covmatrix))
        bic[i,j]  <- -2*(g.result$loglik) + p_off_d*log(nrow(data))
        bic=as.data.frame(bic)
        rownames(bic)=rho
        colnames(bic)=kappa
    }
}
for(i in 1:100){
    for(j in 1:100){
    if(bic[i,j]==which.min(bic)){
        rho=rho[i]
        kappa=kappa[j]
    }
}}
 rho_m=rho*U-kappa*prior
 parameters=list(rho_m,rho,kappa)
return(parameters)
}