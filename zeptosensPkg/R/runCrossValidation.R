#' Run 10-fold Cross validation for data through the current algorythm.
#' 
#' @param data input expression data. Coloumns as the gene, rows as the sample.With colnames as the gene tags, rownames as the sample tags.
#' @param prior prior information matrix with colnames and rownames as the gene tags.
#' @param boot.time Bootstrap time mannually set.
#' @param rho regulization parameter
#' @param kappa scaler parameter
#' @return result lists of positive predictive values, negative predictive values,sensitivity, specificity,for random network and the network predicted with the algorithm.
#' @concept zeptosensPkg
#' @export
runCrossValidation=function(data,prior,boot.time,rho,kappa,cut.off){
    index=colnames(prior[,which(colnames(prior)%in%colnames(data))])#match the data
    
    data=data[,index]
    prior=prior[index,index]
    prior=ifelse(prior!=0,1,0)#information matrix of prior
    
    score1.tp=array(0,dim = c(boot.time,1))
    score1.fp=array(0,dim = c(boot.time,1))
    score1.fn=array(0,dim = c(boot.time,1))
    score1.tn=array(0,dim = c(boot.time,1))

    scorer.tp=array(0,dim = c(boot.time,1))
    scorer.fp=array(0,dim = c(boot.time,1))
    scorer.fn=array(0,dim = c(boot.time,1))
    scorer.tn=array(0,dim = c(boot.time,1))
    
    for(r in 1:boot.time){
        
        #Split the data into two parts
        valid_n=sample(1:nrow(data),0.1*nrow(data),replace=F)
        valid_data=data[valid_n,]
        train_data=data[-valid_n,]
        
        #set up the parameters for regulization
        rho=rho
        kappa=kappa
        U=matrix(1,nrow(prior),ncol(prior))
        rho_m=rho*U-kappa*prior
        pc=cov(train_data) 
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
        #cutoff at cut.off point
        t.edges.rhoadjusted.d=ifelse(abs(t.edges.rhoadjusted.d)>cut.off,t.edges.rhoadjusted.d,0)
        #t.net in signaling way(0,1,-1)
        t.net.rhoadjusted.d=ifelse(col(t.edges.rhoadjusted.d)!=row(t.edges.rhoadjusted.d)& t.edges.rhoadjusted.d>0,1,
                                   ifelse(row(t.edges.rhoadjusted.d)!=col(t.edges.rhoadjusted.d)&t.edges.rhoadjusted.d<0,-1,0))
        #random network 
        random.pcor <- ggm.simulate.pcor(ncol(data),etaA = sum(t.net.rhoadjusted.d!=0)/(ncol(data)^2))
        
        t.net.r=ifelse( col(random.pcor)!=row(random.pcor)&random.pcor>0,1,
                        ifelse(row(random.pcor)!=col(random.pcor)&random.pcor<0,-1,0))
        
        # validation data correlation
        valid.cor=cor(valid_data) #pearson correlation
        
        ##cutoff point as the mean sdv
        mean.v=median(cor(train_data))
        std.v=sd(cor(train_data))
        
        valid.pcor <-ifelse(col(valid.cor)!=row(valid.cor)&valid.cor>=mean.v+std.v,1,
                            ifelse(col(valid.cor)!=row(valid.cor)&valid.cor<=mean.v-std.v,-1,0))
        
        #10 fold cross validation as comparisn of the trainning network with the signal of the valid data  
        score1.1=0
        score1.2=0
        score1.3=0
        score1.4=0
        
        score_r.1=0
        score_r.2=0
        score_r.3=0
        score_r.4=0
        
        #True Positive#
        for(t in 1:nrow(valid.pcor)){
            for(p in 1:ncol(valid.pcor)){
                if(valid.pcor[t,p]==t.net.rhoadjusted.d[t,p]&valid.pcor[t,p]!=0){
                    score1.1=score1.1+1
                }
                if(valid.pcor[t,p]==t.net.r[t,p]&valid.pcor[t,p]!=0){
                    score_r.1=score_r.1+1
                }
            }}
        
        #False Positive
        for(t in 1:nrow(valid.pcor)){
            for(p in 1:ncol(valid.pcor)){
                if(valid.pcor[t,p]!=t.net.rhoadjusted.d[t,p]&t.net.rhoadjusted.d[t,p]!=0){
                    score1.2=score1.2+1
                }
        
                if(valid.pcor[t,p]!=t.net.r[t,p]&t.net.r[t,p]!=0){
                    score_r.2=score_r.2+1
                }
            }}
        
        #False Negative#
        for(t in 1:nrow(valid.pcor)){
            for(p in 1:ncol(valid.pcor)){
                if(valid.pcor[t,p]!=t.net.rhoadjusted.d[t,p]&t.net.rhoadjusted.d[t,p]==0){
                    score1.3=score1.3+1
                }
            
                if(valid.pcor[t,p]!=t.net.r[t,p]&t.net.r[t,p]==0){
                    score_r.3=score_r.3+1
                }
            }}
        
        #True Negative#
        for(t in 1:nrow(valid.pcor)){
            for(p in 1:ncol(valid.pcor)){
                if(valid.pcor[t,p]==t.net.rhoadjusted.d[t,p]&valid.pcor[t,p]==0){
                    score1.4=score1.4+1
                }
            
                if(valid.pcor[t,p]==t.net.r[t,p]&valid.pcor[t,p]==0){
                    score_r.4=score_r.4+1
                }
            }}
        #score in table###and the explaination##3
        
        score1.tp[r]=score1.1
        score1.fp[r]=score1.2
        score1.fn[r]=score1.3
        score1.tn[r]=score1.4
        
        scorer.tp[r]=score_r.1
        scorer.fp[r]=score_r.2
        scorer.fn[r]=score_r.3
        scorer.tn[r]=score_r.4
    }
    
    sensitivity1=mean(score1.tp/(score1.tp+score1.fn))#estimated directional network
    specificity1=mean(score1.tn/(score1.fp+score1.tn))
    ppv1=mean(score1.tp/(score1.tp+score1.fp))
    npv1=mean(score1.tn/(score1.tn+score1.fn))
    fdr1=mean(score1.fp/(score1.tp+score1.fp))
    
    sensitivityr=mean(scorer.tp/(scorer.tp+scorer.fn))#random network
    specificityr=mean(scorer.tn/(scorer.fp+scorer.tn))
    ppvr=mean(scorer.tp/(scorer.tp+scorer.fp))
    npvr=mean(scorer.tn/(scorer.tn+scorer.fn))
    fdrr=mean(scorer.fp/(scorer.tp+scorer.fp))
    
    result=list(sensitivity=sensitivity,specificity=specificity,ppv=ppv1,npv=npv1,fdr=fdr1,
                sensitivityr=sensitivityr,specificityr=specificityr,ppvr=ppvr,npvr=npvr,fdrr=fdrr)
    return(result)
}