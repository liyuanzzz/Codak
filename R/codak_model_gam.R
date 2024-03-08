
gen_model_gam=function(piformula,alphaformula,betaformula,piargs=list(),alphaargs=list(),betaargs=list()){

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
    piargs <- c(list(family=quasibinomial()), piargs)
    param11 <- do.call(safe_gam, piargs)$param

    # 0 to 1, 1 to 1, 2 to 0. Masked left as alpha.
    piargs <- complete_args(x=data, response=ifelse(temph==0,1,2-temph), args=args, weights=ifelse(temph==0,0,censorweights))
    piargs <- c(list(family=quasibinomial()), piargs)
    param12 <- do.call(safe_gam, piargs)$param

    # prob being null, alpha and beta components, respectively
    param1=cbind(param11,(1-param11)*param12,(1-param11)*(1-param12))
    return(param1)
  }
  alphafun_init=function(q,indReveal,data,...){
    args=list(...)

    tempy=ifelse(indReveal,q,sign(q)*pmax(abs(q),1-abs(q)))
    tempy1=-log(1/2+pminmax(tempy,1e-15-1,1-1e-15)/2)

    alpha_args <- complete_args(x=data, response=tempy1, args=args, weights=ifelse(q<0, 1, 0))
    alpha_args <- c(list(family=Gamma(link="log")), alpha_args)
    param2 <- do.call(safe_gam, alpha_args)$param
    param2 <- pmin(1/param2,1)

    return(param2)
  }
  betafun_init=function(q,indReveal,data,...){
    args=list(...)

    tempy=ifelse(indReveal,q,sign(q)*pmax(abs(q),1-abs(q)))
    tempy2=-log(1/2-pminmax(tempy,1e-15-1,1-1e-15)/2)

    beta_args <- complete_args(x=data, response=tempy2, args=args, weights=ifelse(q>0, 1, 0))
    beta_args <- c(list(family=Gamma(link="log")), beta_args)
    param3 <- do.call(safe_gam, beta_args)$param
    param3 <- pmin(1/param3,1)

    return(param3)
  }
  
  
  # in case of only G is given, default other args are null 
  
  pipre <- function(formula=NULL, data=NULL, weights=NULL, ...){
    safe_gam(formula=formula, family = quasibinomial(), data=data, weights = weights, fit=FALSE, ...)
  }

  pifun <- function(formula=NULL, data=NULL, weights=NULL, pre=NULL, init=NULL, ...){
    if(!is.null(pre)){
      pre$w <- weights
      pre$mf$`(weights)` <- weights
    }
    sp=NULL
    if(!is.null(init)){
      sp=init$sp
    }
    safe_gam(formula=formula, family = quasibinomial(), data=data, weights = weights, G=pre, sp=sp, ...)
  }
  
  alphabetapre <- function(formula=NULL, data=NULL, weights=NULL, ...){
    safe_gam(formula=formula, family = Gamma(link="log"), data=data, weights = weights, fit=FALSE, ...)
  }
  
  
  alphabetafun <- function(formula=NULL, data=NULL, weights=NULL, pre=NULL, init=NULL, ...){
    if(!is.null(pre)){
      pre$w <- weights
      pre$mf$`(weights)` <- weights
    }
    sp=NULL
    if(!is.null(init)){
      sp=init$sp
    }
    res <- safe_gam(formula=formula, family = Gamma(link="log"), data=data, weights = weights, G=pre, sp=sp, ...)
    res$param <- pmin(1/res$param,1)
    return(res)
  }

  piargs=c(list(formula=piformula,type="gam"),piargs)
  piargs_init=piargs
  alphaargs=c(list(formula=alphaformula,type="gam"),alphaargs)
  alphaargs_init=alphaargs
  betaargs=c(list(formula=betaformula,type="gam"),betaargs)
  betaargs_init=betaargs
  gen_Codak_model(pifun,alphabetafun,alphabetafun,pifun_init,alphafun_init,betafun_init,pipre,alphabetapre,alphabetapre,
                  piargs,alphaargs,betaargs,piargs_init,alphaargs_init,betaargs_init,name="gam")
}

#'Covariate and direction adaptive knockoff procedure with group structure
#'
#'Codak is used to solve multiple hypothesis testing problems accompanied by generic covariate information.
#'    \code{codak.gam} fits parameters by \code{\link[mgcv]{gam}} from \code{mgcv} package.
#'    This function can not only utilize numerical covaraites but also work under situations where null hypotheses are divided into multiple groups and accompanied by different sets of covariates in different groups.
#'
#'@param p a vector of the p-values.
#'@param signy a vector of the signs of test statistics.
#'@param x a data frame of the accompanying covariates.
#'@param piformula a vector of strings. The formulas for fitting \eqn{\pi_0(x)}{\pi_0(x)} and \eqn{\lambda(x)}{\lambda(x)} by \code{gam} in package \code{mgcv}. Each formula matches one group.
#'@param alphaformula a vector of strings. The formulas for fitting \eqn{\alpha(x)}{\alpha(x)} by \code{gam} in package \code{mgcv}. Each formula matches one group.
#'@param betaformula a vector of strings. The formulas for fitting \eqn{\beta(x)}{\beta(x)} by \code{gam} in package \code{mgcv}. Each formula matches one group.
#'@param group a numerical vector of group labels. It ranges from 1 to \eqn{g}{g}, where \eqn{g}{g} is the number of groups. It is set to NULL by default.
#'@param initparam A matrix. The initial parameter of the EM-algorithm. It is generated from the algorithm described in Appendx B.3 in the paper by default.
#'@param nfits An integer. The maximum number of model-fitting. Supposing \eqn{n}{n} is the number of null hypotheses, the model will be refitted every time after accepting \eqn{\left \lceil{n/nfits}\right \rceil}{[n/nfits]} null hypotheses.
#'@param init logical; if TRUE (default), the smoothing parameters are obtained from the first EM iteration and set as fixed for the following EM iterations. Otherwise, the smoothing parameters are estimated separately in each iteration.
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
#'out=codak.gam(pval,sign(y),x,c("s(x1,x2)","s(x1)"),group=1+x3)
#'rejresult=(out<=alpha)
#'
#'fdr=sum((rejresult)&(temp==0))/max(sum(rejresult),1)
#'power=sum((rejresult)&(temp!=0))/max(sum(temp!=0),1)
#'
#'@export

codak.gam=function(p,signy,x,piformula,alphaformula=piformula,betaformula=alphaformula,piargs=list(),alphaargs=list(),betaargs=list(),group=NULL,initparam=NULL,nfits=20,init=TRUE,gamma=1){
  if(is.null(group)&&max(c(length(piformula),length(alphaformula),length(betaformula)))>1) stop("Please only give one formula with no group setting.")
  if(!is.null(group)&&(length(unique(group))!=length(piformula))) stop("The number of formulas does not match the number of groups.")
  if(is.matrix(x)||is.vector(x)){
    xx=data.frame(x)
  }else{
    xx=x
  }
  if(is.null(group)){
    model=gen_model_gam(piformula,alphaformula,betaformula,piargs,alphaargs,betaargs)
  }else{
    model=gen_model_gam_group(piformula,alphaformula,betaformula,piargs,alphaargs,betaargs,group)
  }

  codak(p=p,signy=signy,x=xx,model=model,param0=initparam,nfits=nfits,initial=init,gamma=gamma)
}
