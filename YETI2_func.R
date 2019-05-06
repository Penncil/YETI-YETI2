
##===================================================================================================================================
## "YETI2.func" is a function that is used to detect genetic interactions via YETI2 method,  and we only consider the dominant model
##
##  Input: data (e.g., ID, Y, Z, SNP1,...,SNPL), L (# of markers), K (# of populations/studies), pre.alpha (tuning parameter), alpha
## Output: SNP1 (id for the 1st interacting SNP), SNP2 (id for the 2nd interacting SNP), bet3.overall.est, beta3.overall.se, pvalue   
##====================================================================================================================================

YETI2.func = function(mydata, L, K, pre.alpha, alpha){
    
    ##************************************
    ## 1.1: Do the subgroup analysis
    ##************************************
    eta1.est = eta1.se = matrix(NA, nrow=K ,ncol=L)
    for(k in 1:K){
        data.temp0 = subset(mydata,subset=Z==k) #data: ID, Y, Z, SNP1,..., SNP2
        for(l in 1:L){
            fit = summary(glm(data.temp0[,2]~1+data.temp0[,(3+l)],data=data.temp0,family=binomial(logit)))
            eta1.est[k,l] = as.numeric(fit$coefficients[2,1])
            eta1.se[k,l]= as.numeric(fit$coefficients[2,2])
        }
    }
    
    ##********************************************************************************************************
    ## 1.2: Test for marginal effects and heterogeneity together, ie, H0:eta1^(1)=eta1^(2)=....=eta1^(k) =0
    ##*********************************************************************************************************
    Test.stat = L.star.star = rep(NA,L)
    
    ## threshold for chi-square distribution with df= K
    cutoff.value = qchisq(df = K, 1-pre.alpha)
    
    ## Fixed-effect model
    for(l in 1:L){
    	    Test.stat[l] = sum((eta1.est[,l]/eta1.se[,l])^2)
        if(Test.stat[l] >= cutoff.value){L.star.star[l]=1}
        if(Test.stat[l] < cutoff.value){L.star.star[l]=0}
    }
    ## Collect the significant genes and count the length of each subset
    L.star.star.subset = which(L.star.star==1); nL.star.star=length(L.star.star.subset)
    alpha.sig = alpha/(nL.star.star*(L-1)) ## multiple testing correction for stage2
    
    
    ##**********************************************************
    ## 1.3: Test for H0: beta3^{overall}=0 via meta-analysis
    ##**********************************************************
    beta3.est = beta3.se = matrix(NA, nrow=K ,ncol=1)
    output = NULL ##collect info about those SNPs with interaction effects
        
    if(nL.star.star>0){
        for(l1 in 1:nL.star.star){
            for(l2 in 1:L){
                if((2+L.star.star.subset[l1])!=(2+l2)){
                    data.temp0 = mydata[,c(2,3+L.star.star.subset[l1],3+l2,3)] #Y,SNP1, SNP2, Z
                    for(k in 1:K){
                        data.temp = subset(data.temp0,subset=Z==k)
                        fit.int = summary(glm(data.temp[,1]~1+data.temp[,2]+data.temp[,3]+ data.temp[,2]*data.temp[,3],data=data.temp,family=binomial(logit)))
                        beta3.est[k] = as.numeric(fit.int$coefficients[4,1])
                        beta3.se[k] = as.numeric(fit.int$coefficients[4,2])
                    }
                    
                    ## calcualte the overall estimate of beta3
                    beta3.overall.se = sqrt(1/sum(1/(beta3.se^2)))
                    beta3.overall.est = sum((1/(beta3.se^2))*beta3.est)/sum(1/(beta3.se^2))
                    
                    ## calculate the p-value based on a Wald test
                    pvalue = 2*pnorm(-abs(beta3.overall.est/beta3.overall.se))
                    if(pvalue<= alpha.sig){
                    output = rbind(output,c(L.star.star.subset[l1],l2,beta3.overall.se,beta3.overall.est,pvalue))
                    }
              }
         }
      }
    }
    return(list(SNP1=output[,1], SNP2=output[,2], beta3.overall.est=output[,3], beta3.overall.se=output[,4], pvalue=output[,5]))
}


