#' Run 10-fold Cross validation for data through the current algorythm.
#' 
#' @param data input expression data. Coloumns as the gene, rows as the observations. With colnames as the gene tags, rownames as the sample barcodes.
#' @param prior prior information matrix with colnames and rownames as the gene tags.
#' @param boot.time Bootstrap time mannually set.
#' @param rho regulization parameter
#' @param kappa scaler parameter
#' @param cut.off the cut.off poingt for edge.value.Default value at 0.05.
#' @return result lists of bootstrap cross validation for random network(scorer) and the network predicted with the algorithm.As in the list(score1).
#' @concept zeptosensPkg
#' @export
runCrossValidation=function(data,prior,boot.time,rho,kappa,cut.off=0.05){
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
        
        # validation data signal
        valid.mean=colMedians(as.matrix(valid_data))
        valid.std=colSds(as.matrix(valid_data))
        
        valid.sig <-ifelse(valid_data>=valid.mean+valid.std,1,
                           ifelse(valid_data<=valid.mean-valid.std,-1,0))
        
        #10 fold cross validation as comparisn of the trainning network with the signal of the valid data  
        score1.1=array(0,dim = c(nrow(valid_data),1))
        score1.2=array(0,dim = c(nrow(valid_data),1))
        
        score_r.1=array(0,dim = c(nrow(valid_data),1))
        score_r.2=array(0,dim = c(nrow(valid_data),1))
       
        #valid_data signal
         for(a in 1:nrow(valid_data)){
            valid.pcor.d=array(0,dim = c(ncol(valid_data),ncol(valid_data)))
            valid.pcor=array(0,dim = c(ncol(valid_data),ncol(valid_data)))
            for(i in 1:ncol(valid_data)){
                for(j in 1:ncol(valid_data)){
                    if(valid.sig[a,i]==1&valid.sig[a,j]==1)
                    {valid.pcor[i,j]=1  
                    valid.pcor.d[i,j]=1}
                    if(valid.sig[a,i]==-1&valid.sig[a,j]==-1)
                    {valid.pcor[i,j]=1
                    valid.pcor.d[i,j]=1}
                    if(valid.sig[a,i]==1&valid.sig[a,j]==-1)
                    {valid.pcor[i,j]=-1
                    valid.pcor.d[i,j]=-1}
                    if(valid.sig[a,i]==-1&valid.sig[a,j]==1)
                    {valid.pcor[i,j]=-1
                    valid.pcor.d[i,j]=-1}
                    if(valid.sig[a,i]==1&valid.sig[a,j]==0)
                    {valid.pcor[i,j]=0
                    valid.pcor.d[i,j]=0}
                    if(valid.sig[a,i]==-1&valid.sig[a,j]==0)
                    {valid.pcor[i,j]=0
                    valid.pcor.d[i,j]=0}
                    if(valid.sig[a,i]==0&valid.sig[a,j]==1)
                    {valid.pcor[i,j]=0
                    valid.pcor.d[i,j]=NA}
                    if(valid.sig[a,i]==0&valid.sig[a,j]==-1)
                    {valid.pcor[i,j]=0
                    valid.pcor.d[i,j]=NA}
                    if(valid.sig[a,i]==0&valid.sig[a,j]==0)
                    {valid.pcor[i,j]=NA
                    valid.pcor.d[i,j]=NA}
                }
            }   
            valid.pcor.d=ifelse( col(valid.pcor.d)==row(valid.pcor.d),NA,valid.pcor.d)
            valid.pcor=ifelse( col(valid.pcor)==row(valid.pcor),NA,valid.pcor)
        
            #True Positive#(#1)
            for(t in 1:ncol(valid_data)){
                for(p in 1:ncol(valid_data)){
                    if(valid.pcor.d[t,p]==t.net.rhoadjusted.d[t,p]&is.na(valid.pcor.d[t,p])==F&t.net.rhoadjusted.d[t,p]!=0){
                        score1.1[a]=score1.1[a]+1
                    }
                    if(valid.pcor[t,p]==t.net.r[t,p]&is.na(valid.pcor[t,p])==F&t.net.r[t,p]!=0){
                        score_r.1[a]=score_r.1[a]+1
                    }
                }}
            
            #False Positive(#0)
            for(t in 1:nrow(valid.pcor)){
                for(p in 1:ncol(valid.pcor)){
                    if(valid.pcor.d[t,p]!=t.net.rhoadjusted.d[t,p]&is.na(valid.pcor.d[t,p])==F&t.net.rhoadjusted.d[t,p]!=0){
                        score1.2[a]=score1.2[a]+1
                    }
                    if(valid.pcor[t,p]!=t.net.r[t,p]&is.na(valid.pcor[t,p])==F&t.net.r[t,p]!=0){
                        score_r.2[a]=score_r.2[a]+1
                    }
                }}
         }  
        #Score for each boot.time
        
        score1.tp[r]=mean(score1.1)
        score1.fp[r]=mean(score1.2)
        
        scorer.tp[r]=mean(score_r.1)
        scorer.fp[r]=mean(score_r.2)
    }
    
    score1=mean(score1.tp/(score1.tp+score1.fp))#rho with directional prior
    scorer=mean(scorer.tp/(scorer.tp+scorer.fp))#rho without prior
    
    result=list(score1=score1,scorer=scorer)
    return(result)
}