#' Run 10-fold Cross validation for data through the current algorythm.
#' 
#' @param data input expression data. Coloumns as the gene, rows as the sample.With colnames as the gene names, rownames as the sample tags.
#' @param boot.time Bootstrap time mannually set.
#' @param rho regulization parameter
#' @param kappa scaler parameter
#' @return result lists of positive predictive values, negative predictive values,sensitivity, specificity,false discovery rate
#' @concept zeptosensPkg
#' @export
runCrossValidation=function(data,boot.time,rho,kappa){
    score1.tp=array(0,dim = c(boot.time,1))
    score1.fp=array(0,dim = c(boot.time,1))
    score1.fn=array(0,dim = c(boot.time,1))
    score1.tn=array(0,dim = c(boot.time,1))
    
    score2.tp=array(0,dim = c(boot.time,1))
    score2.fp=array(0,dim = c(boot.time,1))
    score2.fn=array(0,dim = c(boot.time,1))
    score2.tn=array(0,dim = c(boot.time,1))
    
    score3.tp=array(0,dim = c(boot.time,1))
    score3.fp=array(0,dim = c(boot.time,1))
    score3.fn=array(0,dim = c(boot.time,1))
    score3.tn=array(0,dim = c(boot.time,1))
    
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
        U=matrix(1,nrow(prior1),ncol(prior1))
        rho_m=rho*U-kappa*prior1
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
        #cutoff at 0.05
        t.edges.rhoadjusted.d=ifelse(abs(t.edges.rhoadjusted.d)>0.05,t.edges.rhoadjusted.d,0)
        
        #network with symmetric prior information
        rho=rho
        kappa=kappa
        U=matrix(1,nrow(prior2),ncol(prior2))
        rho_m2=rho*U-kappa*prior2
        
        sigma.matrix<-glasso(pc,rho=rho_m)$wi
        pcor.matrix=matrix(0,nrow=ncol(data),ncol=ncol(data))
        for(i in 1:ncol(data)){
            for(j in 1:ncol(data)){
                pcor.matrix[i,j]=-sigma.matrix[i,j]/sqrt(sigma.matrix[i,i]*sigma.matrix[j,j])
            }
        }
        
        t.edges.rhoadjusted2=pcor.matrix
        t.edges.rhoadjusted2=ifelse(abs(t.edges.rhoadjusted2)>0.05,t.edges.rhoadjusted2,0)
        
        
        #network with no prior information
        sigma.matrix<-glasso(pc,rho=0.04)$wi
        pcor.matrix=matrix(0,nrow=ncol(data),ncol=ncol(data))
        for(i in 1:ncol(data)){
            for(j in 1:ncol(data)){
                pcor.matrix[i,j]=-sigma.matrix[i,j]/sqrt(sigma.matrix[i,i]*sigma.matrix[j,j])
            }
        }
        
        t.edges.rho=pcor.matrix
        t.edges.rho=ifelse(abs(t.edges.rho)>0.05,t.edges.rho,0)
        
        
        
        #t.net in signaling way(0,1)
        t.net.rhoadjusted.d=ifelse(col(t.edges.rhoadjusted.d)!=row(t.edges.rhoadjusted.d)& t.edges.rhoadjusted.d>0,1,
                                   ifelse(row(t.edges.rhoadjusted.d)!=col(t.edges.rhoadjusted.d)&t.edges.rhoadjusted.d<0,-1,0))
        
        t.net.rhoadjusted2=ifelse(col(t.edges.rhoadjusted2)!=row(t.edges.rhoadjusted2)& t.edges.rhoadjusted2>0,1,
                                  ifelse(row(t.edges.rhoadjusted2)!=col(t.edges.rhoadjusted2)&t.edges.rhoadjusted2<0,-1,0))
        
        t.net.rho=ifelse(row(t.edges.rho)!=col(t.edges.rho)& t.edges.rho>0,1,
                         ifelse(row(t.edges.rho)!=col(t.edges.rho)&t.edges.rho<0,-1,0))
        #network random
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
        
        score2.1=0
        score2.2=0
        score2.3=0
        score2.4=0
        
        score3.1=0
        score3.2=0
        score3.3=0
        score3.4=0
        
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
                if(valid.pcor[t,p]==t.net.rhoadjusted2[t,p]&valid.pcor[t,p]!=0){
                    score2.1=score2.1+1
                }
                if(valid.pcor[t,p]==t.net.rho[t,p]&valid.pcor[t,p]!=0){
                    score3.1=score3.1+1
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
                if(valid.pcor[t,p]!=t.net.rhoadjusted2[t,p]&t.net.rhoadjusted2[t,p]!=0){
                    score2.2=score2.2+1
                }
                if(valid.pcor[t,p]!=t.net.rho[t,p]&t.net.rho[t,p]!=0){
                    score3.2=score3.2+1
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
                if(valid.pcor[t,p]!=t.net.rhoadjusted2[t,p]&t.net.rhoadjusted2[t,p]==0){
                    score2.3=score2.3+1
                }
                if(valid.pcor[t,p]!=t.net.rho[t,p]&t.net.rho[t,p]==0){
                    score3.3=score3.3+1
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
                if(valid.pcor[t,p]==t.net.rhoadjusted2[t,p]&valid.pcor[t,p]==0){
                    score2.4=score2.4+1
                }
                if(valid.pcor[t,p]==t.net.rho[t,p]&valid.pcor[t,p]==0){
                    score3.4=score3.4+1
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
        
        score2.tp[r]=score2.1
        score2.fp[r]=score2.2
        score2.fn[r]=score2.3
        score2.tn[r]=score2.4
        
        score3.tp[r]=score3.1
        score3.fp[r]=score3.2
        score3.fn[r]=score3.3
        score3.tn[r]=score3.4
        
        scorer.tp[r]=score_r.1
        scorer.fp[r]=score_r.2
        scorer.fn[r]=score_r.3
        scorer.tn[r]=score_r.4
    }
    
    sensitivity1=mean(score1.tp/(score1.tp+score1.fn))#rho with directional prior
    specificity1=mean(score1.tn/(score1.fp+score1.tn))
    ppv1=mean(score1.tp/(score1.tp+score1.fp))
    npv1=mean(score1.tn/(score1.tn+score1.fn))
    fdr1=mean(score1.fp/(score1.tp+score1.fp))
    
    sensitivity2=mean(score2.tp/(score2.tp+score2.fn))#rho with directional prior
    specificity2=mean(score2.tn/(score2.fp+score2.tn))
    ppv2=mean(score2.tp/(score2.tp+score2.fp))
    npv2=mean(score2.tn/(score2.tn+score2.fn))
    fdr2=mean(score2.fp/(score2.tp+score2.fp))
    
    sensitivity3=mean(score3.tp/(score3.tp+score3.fn))#rho with directional prior
    specificity3=mean(score3.tn/(score3.fp+score3.tn))
    ppv3=mean(score3.tp/(score3.tp+score3.fp))
    npv3=mean(score3.tn/(score3.tn+score3.fn))
    fdr3=mean(score3.fp/(score3.tp+score3.fp))
    
    sensitivityr=mean(scorer.tp/(scorer.tp+scorer.fn))#rho with directional prior
    specificityr=mean(scorer.tn/(scorer.fp+scorer.tn))
    ppvr=mean(scorer.tp/(scorer.tp+scorer.fp))
    npvr=mean(scorer.tn/(scorer.tn+scorer.fn))
    fdrr=mean(scorer.fp/(scorer.tp+scorer.fp))
    
    # p.score1=t.test(n_score1,n_scorer)$p.value
    # p.score2=t.test(n_score2,n_scorer)$p.value
    # p.score3=t.test(n_score3,n_scorer)$p.value
    
    result=list(sensitivity1=sensitivity1,specificity1=specificity1,ppv1=ppv1,npv1=npv1,fdr1=fdr1,
                sensitivity2=sensitivity2,specificity2=specificity2,ppv2=ppv2,npv2=npv2,fdr2=fdr2,
                sensitivity3=sensitivity3,specificity3=specificity3,ppv3=ppv3,npv3=npv3,fdr3=fdr3,
                sensitivityr=sensitivityr,specificityr=specificityr,ppvr=ppvr,npvr=npvr,fdrr=fdrr)
    return(result)
}