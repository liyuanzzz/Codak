gen_model_glmnet=function(piargs=list(),alphaargs=list(),betaargs=list()){
  
  pifun_init=function(q,indReveal,data,...){
    args=list(...)

    # revealed 0; masked negative/left 1; masked positive/right 2
    temph=rep(0,length(q))
    temph[!indReveal & q<=0] <- 1
    temph[!indReveal & q>0] <- 2
    
    # give more weight to revealed ones
    weighttt=ifelse(indReveal,1/(1-2*0.3),1)
    weighttt=weighttt/(min(weighttt)+max(weighttt))
    censorweights=pminmax2(weighttt)
    temph1=ifelse(temph==0,1,0)
    param11=safe_glmnet(data,temph1,family="binomial",weights=censorweights)$param
    temph1=ifelse(temph==0,1,2-temph)
    param12=safe_glmnet(data,temph1,family="binomial",weights=ifelse(temph==0,0,censorweights))$param
    param1=cbind(param11,(1-param11)*param12,(1-param11)*(1-param12))
    return(param1)
  }
  
  alphafun_init=function(q,indReveal,data,...){
    args=list(...)
    
    tempy=ifelse(indReveal,q,sign(q)*pmax(abs(q),1-abs(q)))
    tempy1=-log(1/2+pminmax(tempy,1e-15-1,1-1e-15)/2)

    weightt=ifelse(q<0, 1, 0)
    param2=safe_glmnet(data,tempy1,family="Gamma",weights=weightt)$param
    param2=pmin(1/param2,1)
    return(param2)
  }
  betafun_init=function(q,indReveal,data,...){
    args=list(...)
    
    tempy=ifelse(indReveal,q,sign(q)*pmax(abs(q),1-abs(q)))
    tempy2=-log(1/2-pminmax(tempy,1e-15-1,1-1e-15)/2)

    weightt=ifelse(q>0, 1, 0)
    param3=safe_glmnet(data,tempy2,family="Gamma",weights=weightt)$param
    param3=pmin(1/param3,1)
    return(param3)
  }
  pifun=function(data=NULL, weights=NULL, pre=NULL, init=NULL, ...){
    y=data[,1]
    x=as.matrix(data[,-1])
    param1=safe_glmnet(x,y,family="binomial",weights=weights)$param
    return(list(param=param1))
  }
  alphabetafun=function(data=NULL, weights=NULL, pre=NULL, init=NULL, ...){
    y=data[,1]
    x=as.matrix(data[,-1])
    param2=safe_glmnet(x,y,family="Gamma",weights=weights)$param
    param2=pmin(1/param2,1)
    return(list(param=param2))
  }
  pipre=function(data=NULL, weights=NULL, ...){
    return(NULL)
  }
  alphabetapre=function(data=NULL, weights=NULL, ...){
    return(NULL)
  }

  piargs=c(list(type="glmnet"),piargs)
  piargs_init=piargs
  alphaargs=c(list(type="glmnet"),alphaargs)
  alphaargs_init=alphaargs
  betaargs=c(list(type="glmnet"),betaargs)
  betaargs_init=betaargs
  gen_Codak_model(pifun,alphabetafun,alphabetafun,pifun_init,alphafun_init,betafun_init,pipre,alphabetapre,alphabetapre,
                  piargs,alphaargs,betaargs,piargs_init,alphaargs_init,betaargs_init,name="glmnet")
}

#'Covariate and direction adaptive knockoff procedure with group structure
#'
#'Codak is used to solve multiple hypothesis testing problems accompanied by generic covariate information.
#'    \code{codak.glmnet} fits parameters by \code{\link[glmnet]{glmnet}} from \code{glmnet} package.
#'    This function can not only utilize numerical covaraites but also work under situations where null hypotheses are divided into multiple groups and accompanied by different sets of covariates in different groups.
#'
#'@param p a vector of the p-values.
#'@param signy a vector of the signs of test statistics.
#'@param x a matrix of the accompanying covariates. Each column corresponds to a covariate.
#'@param group a numerical vector of group labels. It ranges from 1 to \eqn{g}{g}, where \eqn{g}{g} is the number of groups. It is set to NULL by default, and by doing so this function will not consider the group strcture.
#'@param groupinfo a logical matrix with number of rows equal to the number of groups and number of columns equal to the number of column in \eqn{x}. A TRUE value in the ith row and jth column indicates that the jth covariate will be considered in the model for the ith group. It is set to NULL by default, which means that each covariate will be considered in all of the groups.
#'@param initparam A matrix. The initial parameter of the EM-algorithm. It is generated from the algorithm described in Appendx B.3 in the paper by default.
#'@param nfits An integer. The maximum number of model-fitting. Supposing \eqn{n}{n} is the number of null hypotheses, the model will be refitted every time after accepting \eqn{\left \lceil{n/nfits}\right \rceil}{[n/nfits]} null hypotheses.
#'@param gamma A number concatenating the second shape parameter for the beta densities. Default to 1.
#'
#'@return A vector containing the q-value of each null hypothesis. The idea of q-values is from Storey (2002).
#'  The rejection rule of \code{codak_group} is equivalent to thresholding on this value.
#'  Setting the nominal FDR level as \eqn{\alpha}{\alpha} will reject and only reject all the null hypotheses with
#'  q-values smaller than or equal to \eqn{\alpha}{\alpha}.
#'
#'@examples
#'m <- 4000
#'mu1 <- -1.5
#'mu2 <- 3
#'omega=0.25
#'alpha=0.1
#' 
#'x1 <- 2*runif(m)-1
#'x2 <- 2*runif(m)-1
#'pi0=0.6+0.3*(1-omega)*(x1^2+x2^2-2/3)-0.375*omega*x1*x2
#'pi1=(1-pi0)*(1+x1*x2)/2
#'pi2=1-pi0-pi1
#'temp=rep(0,m)
#'for(j in 1:m){
#'   temp[j]=sample(c(0,1,2),1,replace=TRUE,prob=c(pi0[j],pi1[j],pi2[j]))
#'}
#'y=rnorm(m)+mu1*(temp==1)+mu2*(temp==2)
#' 
#'x=cbind(x1,x2)
#'pval=2-2*pnorm(abs(y))
#'out=codak.glmnet(pval,sign(y),x)
#'rejresult=(out<=alpha)
#'
#'fdr=sum((rejresult)&(temp==0))/max(sum(rejresult),1)
#'power=sum((rejresult)&(temp!=0))/max(sum(temp!=0),1)
#'
#'@export

codak.glmnet=function(p,signy,x,piargs=list(),alphaargs=list(),betaargs=list(),group=NULL,groupinfo=NULL,initparam=NULL,nfits=20,gamma=1){
  if(!is.null(group)){
    if(is.null(groupinfo)){
      groupinfo=matrix(TRUE,length(unique(group)),ncol(x))
    }else{
      if((!is.logical(groupinfo))||(nrow(groupinfo)!=length(unique(group)))||(ncol(groupinfo)!=ncol(x))) stop("wrong format of groupinfo")
    }
  }
  if(is.data.frame(x)){
    xx=as.matrix(x)
  }else{
    xx=x
  }
  if(is.null(group)){
    model=gen_model_glmnet(piargs,alphaargs,betaargs)
  }else{
    model=gen_model_glmnet_group(piargs,alphaargs,betaargs,group,groupinfo)
  }
  
  codak(p=p,signy=signy,x=xx,model=model,param0=initparam,nfits=nfits,gamma=gamma)
}
