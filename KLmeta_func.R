
##===================================================================================================================================
## "KLmeta.func" is a function that is used to detect genetic interactions via KL-meta method, and we only consider the dominant model
## 
##
##  Input: data (e.g., ID, Y, Z, SNP1,...,SNPL), L (# of markers), K (# of populations/studies), pre.alpha (tuning parameter), alpha
## Output: SNP1 (id for the 1st interacting SNP), SNP2 (id for the 2nd interacting SNP), bet3.overall.est, beta3.overall.se, pvalue   
##====================================================================================================================================

KLmeta.func = function(mydata, L, K, pre.alpha, alpha){
    
    ##************************************
    ## 1.1: Do the subgroup analysis
    ##*************************************
    eta1.est = eta1.se = matrix(NA, nrow=K ,ncol=L)
    for(k in 1:K){
        data.temp0=subset(mydata,subset=Z==k) ##ID,Y,Z,SNP1,...SNPL
        for(l in 1:L){
            fit = summary(glm(data.temp0[,2]~1+data.temp0[,(3+l)],data=data.temp0,family=binomial(logit)))
            eta1.est[k,l] = as.numeric(fit$coefficients[2,1])
            eta1.se[k,l]= as.numeric(fit$coefficients[2,2])
        }
    }
    
    ##*********************************************************
    ## 1.2: Test for K0: eta1^{overall}=0 for each marker
    ##**********************************************************
    eta1.overall.se = eta1.overall.est = pvalue0 = L.star = rep(NA,L)
    
    for(l in 1:L){
       eta1.overall.se[l] = sqrt(1/sum(1/(eta1.se[,l]^2)))
       eta1.overall.est[l] = sum((1/eta1.se[,l]^2)*eta1.est[,l])/sum(1/(eta1.se[,l]^2))
       pvalue0[l] = 2*pnorm(-abs(eta1.overall.est[l]/eta1.overall.se[l])) ## conduct the Wald test
       
       if(pvalue0[l]> pre.alpha){L.star[l] = 0}
       if(pvalue0[l]<= pre.alpha){L.star[l] = 1}
    }
    
    ## Collect the significant genes and count the length of each subset
    L.star.subset = which(L.star==1); nL.star=length(L.star.subset)
    alpha.sig = alpha/(nL.star*(nL.star-1)/2) ##multiple testing correction in stage2
    
    ##*****************************************
    ## 1.3: Test for H0: beta3^{overall}=0
    ##*****************************************
    beta3.est = beta3.se = matrix(NA, nrow=K ,ncol=1)
    output = NULL ##collect info about those SNPs with interaction effects
    
    if(nL.star>=2){
        for(l1 in 1:(nL.star-1)){
            for(l2 in (l1+1):nL.star){
                    data.temp1 = mydata[,c(2,3+L.star.subset[l1],3+L.star.subset[l2],3)] #ID, SNP1, SNP2, Z
                    for(k in 1:K){
                        data.temp=subset(data.temp1,subset=Z==k) 
                        fit.int = summary(glm(data.temp[,1]~1+data.temp[,2]+data.temp[,3]+ data.temp[,2]*data.temp[,3],data=data.temp,family=binomial(logit)))
                        beta3.est[k] = as.numeric(fit.int$coefficients[4,1])
                        beta3.se[k] = as.numeric(fit.int$coefficients[4,2])
                    }
                    
                    ## calculate estiamte of overall beta3
                    beta3.overall.se = sqrt(1/sum(1/(beta3.se^2)))
                    beta3.overall.est = sum((1/(beta3.se^2))*beta3.est)/sum(1/(beta3.se^2))
                    
                    ## calculate p-value based on a Wald test
                    pvalue = 2*pnorm(-abs(beta3.overall.est/beta3.overall.se)) 
                    if(pvalue<= alpha.sig){
                    output = rbind(output,c(L.star.subset[l1],L.star.subset[l2],beta3.overall.se,beta3.overall.est,pvalue))             
                   }
            
              }
       }
    }
    return(list(SNP1=output[,1], SNP2=output[,2], beta3.overall.est=output[,3], beta3.overall.se=output[,4], pvalue=output[,5]))
}


