
gen_model_gam_group=function(piformula,alphaformula,betaformula,piargs=list(),alphaargs=list(),betaargs=list(),group){
  
  pifun_init=function(q,indReveal,data,...){
    args=list(...)
    
    if((!"formula" %in% names(args))&&(!"piformula" %in% names(args))){
      stop("Formula is not found.")
    }
    if(!"piformula" %in% names(args)){
      piformula=args$formula
    }else{
      piformula=args$piformula
    }
    if(!"group" %in% names(args)){
      stop("Group is not found.")
    }
    group<-args$group
    # revealed 0; masked negative/left 1; masked positive/right 2
    temph=rep(0,length(q))
    temph[!indReveal & q<=0] <- 1
    temph[!indReveal & q>0] <- 2
    
    # give more weight to revealed ones
    weighttt=ifelse(indReveal,1/(1-2*0.3),1)
    weighttt=weighttt/(min(weighttt)+max(weighttt))
    censorweights=pminmax2(weighttt)
    temph1=ifelse(temph==0,1,0)
    

    param11=rep(0,length(group))
    # fit by group
    uniGroup <- unique(group)
    for(ijk in uniGroup){
      args$formula <- piformula[ijk]
      piargs <- complete_args(x=data[group==ijk,], response=temph1[group==ijk], args=args, weights=censorweights[group==ijk])
      piargs <- c(list(family=quasibinomial()), piargs)
      param11[group==ijk] <- do.call(safe_gam, piargs)$param
    }
  
    temph1=ifelse(temph==0,1,2-temph)
    censorweights=ifelse(temph==0,0,censorweights)

    param12=rep(0,length(group))
    for(ijk in uniGroup){
      args$formula <- piformula[ijk]
      piargs <- complete_args(x=data[group==ijk,], response=temph1[group==ijk], args=args, weights=censorweights[group==ijk])
      piargs <- c(list(family=quasibinomial()), piargs)
      param12[group==ijk] <- do.call(safe_gam, piargs)$param
    }
    
    param1=cbind(param11,(1-param11)*param12,(1-param11)*(1-param12))
    return(param1)
  }
  alphafun_init=function(q,indReveal,data,...){
    args=list(...)

    if((!"formula" %in% names(args))&&(!"alphaformula" %in% names(args))&&(!"muformula" %in% names(args))){
      stop("Formula is not found.")
    }
    if("alphaformula" %in% names(args)){
      alphaformula=args$alphaformula
    }else if("muformula" %in% names(args)){
      alphaformula=args$muformula
    }else{
      alphaformula=args$formula
    }
    if(!"group" %in% names(args)){
      stop("Group is not found.")
    }
    group<-args$group

    tempy=ifelse(indReveal,q,sign(q)*pmax(abs(q),1-abs(q)))
    tempy1=-log(1/2+pminmax(tempy,1e-15-1,1-1e-15)/2)
    weightt=ifelse(q<0, 1, 0)

    uniGroup <- unique(group)
    param2=rep(0,length(group))
    for(ijk in uniGroup){
      args$formula <- alphaformula[ijk]
      alpha_args <- complete_args(x=data[group==ijk,], response=tempy1[group==ijk], args=args, weights=weightt[group==ijk])
      alpha_args <- c(list(family=Gamma(link="log")), alpha_args)
      param2[group==ijk] <- do.call(safe_gam, alpha_args)$param
    }
    
    param2=pmin(1/param2,1)
    return(param2)
  }
  betafun_init=function(q,indReveal,data,...){
    args=list(...)

    if((!"formula" %in% names(args))&&(!"betaformula" %in% names(args))&&(!"muformula" %in% names(args))){
      stop("Formula is not found.")
    }
    if("betaformula" %in% names(args)){
      betaformula=args$betaformula
    }else if("muformula" %in% names(args)){
      betaformula=args$muformula
    }else{
      betaformula=args$formula
    }
    if(!"group" %in% names(args)){
      stop("Group is not found.")
    }
    group<-args$group

    tempy=ifelse(indReveal,q,sign(q)*pmax(abs(q),1-abs(q)))
    tempy2=-log(1/2-pminmax(tempy,1e-15-1,1-1e-15)/2)
    weightt=ifelse(q>0, 1, 0)

    uniGroup <- unique(group)
    param3=rep(0,length(group))
    for(ijk in uniGroup){
      args$formula <- betaformula[ijk]
      beta_args <- complete_args(x=data[group==ijk,], response=tempy2[group==ijk], args=args, weights=weightt[group==ijk])
      beta_args <- c(list(family=Gamma(link="log")), beta_args)
      param3[group==ijk] <- do.call(safe_gam, beta_args)$param
    }
    param3=pmin(1/param3,1)
    return(param3)
  }
  
  pipre <- function(formula=NULL, data=NULL, weights=NULL, ...){
    args=list(...)
    if(!"group" %in% names(args)){
      stop("Group is not found.")
    }
    group<-args$group

    Ginfo <- vector(mode="list",length=length(unique(group)))
    # fit by group
    uniGroup <- unique(group)
    for(ijk in uniGroup){
      indGroup <- rep(group==ijk, times=2)
      res <- safe_gam(formula=formula[ijk], family = quasibinomial(), data=data[indGroup,], weights = weights[indGroup], fit=FALSE, ...)
      Ginfo[[ijk]] <- res
    }
    return(Ginfo)
  }
  
  pifun=function(formula=NULL, data=NULL, weights=NULL, pre=NULL, init=NULL, ...){
    args=list(...)
    if(!"group" %in% names(args)){
      stop("Group is not found.")
    }
    group<-args$group
    
    sp=NULL
    if(!is.null(init)){
      sp=init
    }
    
    param1=rep(0,2*length(group))
    info <- vector(mode="list",length=length(unique(group)))
    # fit by group
    uniGroup <- unique(group)
    for(ijk in uniGroup){
      indGroup <- rep(group==ijk, times=2)
      if(!is.null(pre[[ijk]])){
        Gijk <- pre[[ijk]]
        Gijk$w <- weights[indGroup]
        Gijk$mf$`(weights)` <- weights[indGroup]
      }
      spijk=sp[[ijk]]$sp
      res <- safe_gam(formula=formula[ijk], family = quasibinomial(), data=data[indGroup,], weights = weights[indGroup], G=Gijk, sp=spijk, ...)
      param1[indGroup] <- res$param
      info[[ijk]] <- res$info
    }
    return(list(param=param1,info=info))
  }
  
  alphabetapre <- function(formula=NULL, data=NULL, weights=NULL, ...){
    args=list(...)
    if(!"group" %in% names(args)){
      stop("Group is not found.")
    }
    group<-args$group
    
    Ginfo <- vector(mode="list",length=length(unique(group)))
    # fit by group
    uniGroup <- unique(group)
    for(ijk in uniGroup){
      indGroup <- rep(group==ijk, times=2)
      res <- safe_gam(formula=formula[ijk], family = Gamma(link="log"), data=data[indGroup,], weights = weights[indGroup], fit=FALSE, ...)
      Ginfo[[ijk]] <- res
    }
    return(Ginfo)
  }
  
  
  alphabetafun <- function(formula=NULL, data=NULL, weights=NULL, pre=NULL, init=NULL, ...){
    args=list(...)
    if(!"group" %in% names(args)){
      stop("Group is not found.")
    }
    group<-args$group
    
    sp=NULL
    if(!is.null(init)){
      sp=init
    }
    
    param2=rep(0,2*length(group))
    info <- vector(mode="list",length=length(unique(group)))
    # fit by group
    uniGroup <- unique(group)
    for(ijk in uniGroup){
      indGroup <- rep(group==ijk, times=2)
      if(!is.null(pre[[ijk]])){
        Gijk <- pre[[ijk]]
        Gijk$w <- weights[indGroup]
        Gijk$mf$`(weights)` <- weights[indGroup]
      }
      spijk=sp[[ijk]]$sp
      res<- safe_gam(formula=formula[ijk], family = Gamma(link="log"), data=data[indGroup,], weights = weights[indGroup], G=Gijk, sp=spijk, ...)
      param2[indGroup] <- res$param
      info[[ijk]] <- res$info
    }
    param2 <- pmin(1/param2,1)
    return(list(param=param2,info=info))
  }


  piargs=c(list(piformula=piformula,formula=piformula,group=group,type="gam"),piargs)
  piargs_init=piargs
  alphaargs=c(list(alphaformula=alphaformula,formula=alphaformula,group=group,type="gam"),alphaargs)
  alphaargs_init=alphaargs
  betaargs=c(list(betaformula=betaformula,formula=betaformula,group=group,type="gam"),betaargs)
  betaargs_init=betaargs
  gen_Codak_model(pifun,alphabetafun,alphabetafun,pifun_init,alphafun_init,betafun_init,pipre,alphabetapre,alphabetapre,
                  piargs,alphaargs,betaargs,piargs_init,alphaargs_init,betaargs_init,name="gam")
}
