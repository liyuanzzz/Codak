
gen_model_xgboost=function(piargs=list(),alphaargs=list(),betaargs=list(),xgparam){
  
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
    
    # null prob, predicting indReveal? Treat revealed as null, OK
    piargs <- complete_args(x=data, response=ifelse(temph==0,1,0), args=args, weights=censorweights)
    piargs <- c(list(family="binomial"), piargs)
    param11 <- do.call(safe_xgboost, piargs)$param
    
    # 0 to 1, 1 to 1, 2 to 0. Masked left as alpha.
    piargs <- complete_args(x=data, response=ifelse(temph==0,1,2-temph), args=args, weights=ifelse(temph==0,0,censorweights))
    piargs <- c(list(family="binomial"), piargs)
    param12 <- do.call(safe_xgboost, piargs)$param
    
    # prob being null, alpha and beta components, respectively
    param1=cbind(param11,(1-param11)*param12,(1-param11)*(1-param12))
    return(param1)
  }
  alphafun_init=function(q,indReveal,data,...){
    args=list(...)
    
    tempy=ifelse(indReveal,q,sign(q)*pmax(abs(q),1-abs(q)))
    tempy1=-log(1/2+pminmax(tempy,1e-15-1,1-1e-15)/2)
    
    alpha_args <- complete_args(x=data, response=tempy1, args=args, weights=ifelse(q<0, 1, 0))
    alpha_args <- c(list(family="Gamma"), alpha_args)
    param2 <- do.call(safe_xgboost, alpha_args)$param
    param2 <- pmin(1/param2,1)
    
    return(param2)
  }
  betafun_init=function(q,indReveal,data,...){
    args=list(...)
    
    tempy=ifelse(indReveal,q,sign(q)*pmax(abs(q),1-abs(q)))
    tempy2=-log(1/2-pminmax(tempy,1e-15-1,1-1e-15)/2)
    
    beta_args <- complete_args(x=data, response=tempy2, args=args, weights=ifelse(q>0, 1, 0))
    beta_args <- c(list(family="Gamma"), beta_args)
    param3 <- do.call(safe_xgboost, beta_args)$param
    param3 <- pmin(1/param3,1)
    
    return(param3)
  }
  
  pifun <- function(data=NULL, weights=NULL, eta=NULL, max_depth=NULL, init=NULL,...){
    safe_xgboost(data, family = "binomial", weights = weights, eta=eta, max_depth=max_depth, xgb_model=init,...)
  }
  
  alphabetafun <- function(data=NULL, weights=NULL, eta=NULL, max_depth=NULL, init=NULL,...){
    res <- safe_xgboost(data, family = "Gamma", weights = weights, eta=eta, max_depth=max_depth, xgb_model=init,...)
    res$param <- pmin(1/res$param,1)
    return(res)
  }
  
  pipre=function(data=NULL, weights=NULL, ...){
    return(NULL)
  }
  alphabetapre=function(data=NULL, weights=NULL, ...){
    return(NULL)
  }
  
  piargs=c(list(type="xgboost",eta=xgparam[1],max_depth=xgparam[2]),piargs)
  piargs_init=piargs
  alphaargs=c(list(type="xgboost",eta=xgparam[1],max_depth=xgparam[2]),alphaargs)
  alphaargs_init=alphaargs
  betaargs=c(list(type="xgboost",eta=xgparam[1],max_depth=xgparam[2]),betaargs)
  betaargs_init=betaargs
  gen_Codak_model(pifun,alphabetafun,alphabetafun,pifun_init,alphafun_init,betafun_init,pipre,alphabetapre,alphabetapre,
                  piargs,alphaargs,betaargs,piargs_init,alphaargs_init,betaargs_init,name="xgboost")
}


#'Covariate and direction adaptive knockoff procedure with group structure
#'
#'Codak is used to solve multiple hypothesis testing problems accompanied by generic covariate information.
#'    \code{codak.xgboost} fits parameters by \code{\link[xgboost]{xgboost}} from \code{xgboost} package.
#'    This function can not only utilize numerical covaraites but also work under situations where null hypotheses are divided into multiple groups and accompanied by different sets of covariates in different groups.
#'
#'@param p a vector of the p-values.
#'@param signy a vector of the signs of test statistics.
#'@param x a data frame of the accompanying covariates.
#'@param group a numerical vector of group labels. It ranges from 1 to \eqn{g}{g}, where \eqn{g}{g} is the number of groups. It is set to NULL by default.
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
#'x3 <- sample(c(0,1),m,replace=TRUE,prob=c(0.7,0.3))
#'pi0=0.6+0.3*(1-omega)*(x1^2+x2^2-2/3)-0.375*omega*x1*x2
#'pi1=(1-pi0)*(1+x1*x2)/2
#'pi2=1-pi0-pi1
#'temp=rep(0,m)
#'for(j in 1:m){
#'  temp[j]=sample(c(0,1,2),1,replace=TRUE,prob=c(pi0[j],pi1[j],pi2[j]))
#'}
#'y=rnorm(m)+mu1*(temp==1)+mu2*(temp==2)
#'
#'x=data.frame(x1=x1,x2=x2)
#'pval=2-2*pnorm(abs(y))
#'out=codak.xgboost(pval,sign(y),x)
#'rejresult=(out<=alpha)
#'
#'fdr=sum((rejresult)&(temp==0))/max(sum(rejresult),1)
#'power=sum((rejresult)&(temp!=0))/max(sum(temp!=0),1)
#'
#'@export

#'@param eta A number or a set of candidate numbers for model selection. The learning rate used in \code{xgboost} package. Default to 0.3. 
#'@param max_depth A number or a set of candidate numbers for model selection. The maximum depth of the trees used in \code{xgboost} package. Default to 6. 


codak.xgboost=function(p,signy,x,piargs=list(),alphaargs=list(),betaargs=list(),group=NULL,initparam=NULL,nfits=20,gamma=1){
  # if(is.null(group)&&max(c(length(piformula),length(alphaformula),length(betaformula)))>1) stop("Please only give one formula with no group setting.")
  # if(!is.null(group)&&(length(unique(group))!=length(piformula))) stop("The number of formulas does not match the number of groups.")
   if (!requireNamespace("xgboost", quietly = TRUE)){
    stop("'xgboost' package is required for 'adapt_xgboost'. Please intall.")
  }
  if(is.matrix(x)||is.vector(x)){
    xx=data.frame(x)
  }else{
    xx=x
  }
  # if(is.null(group)){
  #   model=gen_model_gam(piformula,alphaformula,betaformula,piargs,alphaargs,betaargs)
  # }else{
  #   model=gen_model_gam_group(piformula,alphaformula,betaformula,piargs,alphaargs,betaargs,group)
  # }
  # 
  eta=seq(0.001,1,0.001);max_depth=c(1:6)
  if(length(eta)>1 || length(max_depth)>1){
    param_grid <-  expand.grid(eta,max_depth)
    colnames(param_grid) <- c("eta","max_depth")
    model <- apply(param_grid, 1, gen_model_xgboost,piargs=piargs,alphaargs=alphaargs,betaargs=betaargs)
    codak(p=p,signy=signy,x=xx,model=model,param0=initparam,nfits=nfits,gamma=gamma,ms=TRUE)
  }else{
    model=gen_model_xgboost(piargs=piargs,alphaargs=alphaargs,betaargs=betaargs,xgparam=c(eta,max_depth))
    codak(p=p,signy=signy,x=xx,model=model,param0=initparam,nfits=nfits,gamma=gamma)
  }
}
