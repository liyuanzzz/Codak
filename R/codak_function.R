emiter=function(param,x,q,indReveal,model,preinfo=NULL,initial=NULL,gamma){
  pi0=param[,1]
  pi1=param[,2]
  pi2=param[,3]
  alph=param[,4]
  bet=param[,5]
  n <- length(q)
  mas=(1:n)[!indReveal]
  aa=rep(0,n)
  bb=rep(0,n)
  cc=rep(0,n)
  dd=rep(0,n)
  tildq=ifelse(q>0,1-q,-1-q)
  # alpha1 <- alph * (1/2+q/2)^(alph-1) # alpha component density
  # alpha2 <- alph * (1/2+tildq/2)^(alph-1) # alpha component density for q_tilda, the mirror
  # beta1 <- bet * (1/2-q/2)^(bet-1) # beta component density
  # beta2 <- bet * (1/2-tildq/2)^(bet-1) # beta component density for q_tilda, the mirror
  
  alpha1 <- alph * (1/2+q/2)^(alph-1) * (1/2-q/2)^(gamma-1) # alpha component density
  alpha2 <- alph * (1/2+tildq/2)^(alph-1) * (1/2-tildq/2)^(gamma-1) # alpha component density for q_tilda, the mirror
  beta1 <- bet * (1/2-q/2)^(bet-1) * (1/2+q/2)^(gamma-1) # beta component density
  beta2 <- bet * (1/2-tildq/2)^(bet-1) * (1/2+tildq/2)^(gamma-1) # beta component density for q_tilda, the mirror
  
  # reveal indices
  reve_true_indices <- which(indReveal)

  # revealed
  denominator1 <- pi0[reve_true_indices] + pi1[reve_true_indices] * alpha1[reve_true_indices] + pi2[reve_true_indices] * beta1[reve_true_indices]
  aa[reve_true_indices] <- pi1[reve_true_indices] * alpha1[reve_true_indices] / denominator1 # posterior prob being alpha
  bb[reve_true_indices] <- pi2[reve_true_indices] * beta1[reve_true_indices] / denominator1 #posterior prob being beta
  cc[reve_true_indices] <- aa[reve_true_indices]
  dd[reve_true_indices] <- bb[reve_true_indices]

  # masked
  reve_false_indices <- which(!indReveal)
  denominator2 <- 2 * pi0[reve_false_indices] + pi1[reve_false_indices] * (alpha1[reve_false_indices] + alpha2[reve_false_indices]) + pi2[reve_false_indices] * (beta1[reve_false_indices] + beta2[reve_false_indices])
  aa[reve_false_indices] <- pi1[reve_false_indices] * (alpha1[reve_false_indices] + alpha2[reve_false_indices]) / denominator2 # posterior prob being alpha
  bb[reve_false_indices] <- pi2[reve_false_indices] * (beta1[reve_false_indices] + beta2[reve_false_indices]) / denominator2
  cc[reve_false_indices] <- pi1[reve_false_indices] * alpha1[reve_false_indices] / denominator2 # posterior prob being alpha and q (instead of tildq)
  dd[reve_false_indices] <- pi2[reve_false_indices] * beta1[reve_false_indices] / denominator2

  # fit pifun twice, can them be combined? Yes, gam has multinom family, not sure if faster.
  h_aug <- rep(c(1,0),rep(n,2))
  x_aug <- rbind(x,x)
  

  # null vs alt
  weights=c(pminmax2(1-aa-bb),pminmax2(aa+bb))
  piargs <- complete_args(x=x_aug, response=h_aug, args=model$args$piargs, weights=weights)
  piargs <- c(list(pre=preinfo$pipre,init=preinfo$init[[1]]),piargs)
  out11=do.call(model$algo$pifun,piargs)
  param11 <- out11$param

  
  # alpha vs beta
  tempweightt=c(pminmax2(aa),pminmax2(bb))
  weights=tempweightt/sum(tempweightt)*n
  piargs2 <- complete_args(x=x_aug, response=h_aug, args=model$args$piargs, weights=weights)
  piargs2 <- c(list(pre=preinfo$pipre,init=preinfo$init[[2]]),piargs2)
  out12=do.call(model$algo$pifun,piargs2)
  param12 <- out12$param

  # probability for null, alpha, beta
  param1=cbind(param11,(1-param11)*param12,(1-param11)*(1-param12))

  # -log(p) transformation
  y=c(-log(1/2+q/2),-log(1/2+tildq/2))
  z=c(-log(1/2-q/2),-log(1/2-tildq/2))
  
  weights=c(pminmax2(cc),pminmax2(aa-cc))
  alphaargs <- complete_args(x=x_aug, response=y, args=model$args$alphaargs, weights=weights)
  alphaargs <- c(list(pre=preinfo$alphapre,init=preinfo$init[[3]]),alphaargs)
  out2=do.call(model$algo$alphafun,alphaargs)
  param2 <- out2$param

  weights=c(pminmax2(dd),pminmax2(bb-dd))
  betaargs <- complete_args(x=x_aug, response=z, args=model$args$betaargs, weights=weights)
  betaargs <- c(list(pre=preinfo$betapre,init=preinfo$init[[4]]),betaargs)
  out3=do.call(model$algo$betafun,betaargs)
  param3 <- out3$param

  
  param_new=cbind(param1[1:n,],param2[1:n],param3[1:n])
  
  # if(isTRUE(initial)||ms){
  #   init11=out11$info
  #   init12=out12$info
  #   init2=out2$info
  #   init3=out3$info
  #   init=list(init11,init12,init2,init3)
  #   return(list(param=param_new,init=init))
  # }
  # return(param_new)
  init11=out11$info
  init12=out12$info
  init2=out2$info
  init3=out3$info
  init=list(init11,init12,init2,init3)
  return(list(param=param_new,init=init))
}

logl=function(param,x,q,indReveal,gamma){
  pi0=param[,1]
  pi1=param[,2]
  pi2=param[,3]
  alph=param[,4]
  bet=param[,5]
  revv=which(indReveal)
  mas=setdiff(1:length(q),revv)
  tildq=ifelse(q>0,1-q,-1-q)
  temp1=exp(log(alph)-log(1/2+q/2)*(1-alph)-log(1/2-q/2)*(1-gamma))
  temp2=exp(log(bet)-log(1/2-q/2)*(1-bet)-log(1/2+q/2)*(1-gamma))
  temp3=exp(log(alph)-log(1/2+tildq/2)*(1-alph)-log(1/2-tildq/2)*(1-gamma))
  temp4=exp(log(bet)-log(1/2-tildq/2)*(1-bet)-log(1/2+tildq/2)*(1-gamma))
  sum(log(1e-15+pi0[revv]+pi1[revv]*temp1[revv]+pi2[revv]*temp2[revv]))+sum(log(1e-15+2*pi0[mas]+pi1[mas]*(temp1[mas]+temp3[mas])+pi2[mas]*(temp2[mas]+temp4[mas])))
}

modelfitting=function(initialparam,x,q,indReveal,model,initial,gamma,ms=FALSE){
  trialparam=initialparam
  #i=0
  newlogl=logl(trialparam,x,q,indReveal,gamma)
  oldlogl=length(q)*log(2)
  if(ms==FALSE){
    cat("Start: old", oldlogl, "; new", newlogl, "\n")
  }
  
  # todo: initial fit to obtain G's
  # args <- model$emInit(trialparam,x,q,indReveal,model)
  
  #preinfo pi, alpha, beta
  n <- length(q)
  h_aug <- rep(c(1,0),rep(n,2))
  x_aug <- rbind(x,x)
  piargs <- complete_args(x=x_aug, response=h_aug, args=model$args$piargs, weights=runif(length(h_aug)))
  pipre=do.call(model$algo$pipre,piargs)

  tildq=ifelse(q>0,1-q,-1-q)
  y=c(-log(1/2+q/2),-log(1/2+tildq/2))
  z=c(-log(1/2-q/2),-log(1/2-tildq/2))

  alphaargs <- complete_args(x=x_aug, response=y, args=model$args$alphaargs, weights=runif(length(h_aug)))
  alphapre=do.call(model$algo$alphapre,alphaargs)


  betaargs <- complete_args(x=x_aug, response=z, args=model$args$betaargs, weights=runif(length(h_aug)))
  betapre=do.call(model$algo$betapre,betaargs)
  
  pre=list()
  pre <- c(list(pipre=pipre),list(alphapre=alphapre),list(betapre=betapre))
  
  
  i=0
  
  # for the init emiter, we record the info and pass to the following emiter, in gam model
  
  # if(isTRUE(initial)){
  #   oldlogl=newlogl
  #   trialparam_old <- trialparam
  #   out=emiter(trialparam,x,q,indReveal,model,preinfo=pre,initial=initial,gamma=gamma)
  #   trialparam=out$param
  #   init=out$init
  #   pre <- c(pre,list(init=init))
  #   newlogl=logl(trialparam,x,q,indReveal,gamma)
  #   i=1
  #   cat(i, ": old", oldlogl, "; new", newlogl, "; max param", apply(abs(trialparam_old-trialparam), 2, max), "\n")
  #   if(newlogl < oldlogl){
  #     cat("log-likelihood decrease, revert to the previous pamameters.\n")
  #     i=10
  #     return(trialparam_old)
  #   }
  # }

  init<-NULL
  preinfo <- c(pre,list(init=init))
  while(abs(newlogl-oldlogl)>0.001&&i<10){
    oldlogl=newlogl
    trialparam_old <- trialparam
    out=emiter(trialparam,x,q,indReveal,model,preinfo=preinfo,initial=initial,gamma=gamma)
    trialparam=out$param
    init=out$init
    preinfo <- c(pre,list(init=init))
    newlogl=logl(trialparam,x,q,indReveal,gamma)
    i=i+1
    if(ms==FALSE){
      cat(i, ": old", oldlogl, "; new", newlogl, "; max param", apply(abs(trialparam_old-trialparam), 2, max), "\n")
    }
    if(newlogl < oldlogl){
      if(ms){
        return(init)
      }
      cat("log-likelihood decrease, revert to the previous pamameters.\n")
      return(trialparam_old)
    }
  }
  if(ms){
    return(init)
  }
  return(trialparam)
}

calcLfdr=function(param,q,gamma){
  pi0=param[,1]
  pi1=param[,2]
  pi2=param[,3]
  alph=param[,4]
  alph=pmin(alph,1-0.01)
  bet=param[,5]
  bet=pmin(bet,1-0.01)
  tildq=ifelse(q>0,1-q,-1-q)
  # temp1=exp(log(alph)-log(1/2+q/2)*(1-alph))
  # temp2=exp(log(bet)-log(1/2-q/2)*(1-bet))
  # temp3=exp(log(alph)-log(1/2+tildq/2)*(1-alph))
  # temp4=exp(log(bet)-log(1/2-tildq/2)*(1-bet))
  alpha1 <- alph * (1/2+q/2)^(alph-1) * (1/2-q/2)^(gamma-1)
  alpha2 <- alph * (1/2+tildq/2)^(alph-1) * (1/2-tildq/2)^(gamma-1)
  beta1 <- bet * (1/2-q/2)^(bet-1) * (1/2+q/2)^(gamma-1)
  beta2 <- bet * (1/2-tildq/2)^(bet-1) * (1/2+tildq/2)^(gamma-1)
  
  # raw estimate
  # 2*pi0/(2*pi0+pi1*(alpha1+alpha2)+pi2*(beta1+beta2))
  
  # over estimate
  #(2*pi0+2*pi1*alph+2*pi2*bet)/(2*pi0+pi1*(alpha1+alpha2)+pi2*(beta1+beta2))
  (2*pi0+2*pi1*alph*1^(alph-1)*0^(gamma-1)+2*pi2*bet*1^(bet-1)*0^(gamma-1))/(2*pi0+pi1*(alpha1+alpha2)+pi2*(beta1+beta2))
  
}

hatFDP=function(A, R){
  (1+A)/pmax(R, 1)
}

# to do: move to model specific function, allow tailored initialization, such as G in mgcv
initialparameter=function(q,indReveal,x,model){
  piargs_init=c(list(q=q,indReveal=indReveal,data=x),model$args$piargs_init)
  param1=do.call(model$algo$pifun_init,piargs_init)
  alphaargs_init=c(list(q=q,indReveal=indReveal,data=x),model$args$alphaargs_init)
  param2=do.call(model$algo$alphafun_init,alphaargs_init)
  betaargs_init=c(list(q=q,indReveal=indReveal,data=x),model$args$betaargs_init)
  param3=do.call(model$algo$betafun_init,betaargs_init)
  cbind(param1,param2,param3)
}

## Main steps
## 1. Initial reveal and parameter setup
## 2. Batch mode: reveal and update parameters

codak <- function(p,signy,x,model,nfits=20, s0=0.3, param0=NULL,initial=NULL,gamma,ms=FALSE){

  # Compute signed p-values. If signy > 0, it multiplies (1-p) by 1, else by -1.
  signp=ifelse(signy>0,1,-1)*(1-p)

  # If x is a vector, convert it to a dataframe
  if(is.vector(x)) x=data.frame(x)

  n <- length(p)
  # Clip the signed p-values to be within the range [-0.99999, 0.99999], to change?
  q=pmax(pmin(signp,1-1e-10),1e-10-1)

  # Initialization
  indReveal <- (abs(q)<=1-s0) & (abs(q)>=s0)
  nInitialReveal <- sum(indReveal)
  revealOrder <- rep(-1,n)
  # initial reveal order by p-value, does not matter as they will not be rejected in the end
  indexReveal <- which(indReveal)
  revealOrder[1:nInitialReveal] <- indexReveal[order(p[indexReveal])]
  # compute initial FDP
  indMask <- !indReveal
  initialA <- sum(indMask & (abs(q) < s0))
  initialR <- sum(indMask & (abs(q) > 1-s0))
  initialFDP <- hatFDP(initialA, initialR)

  # # Compute initial FDR estimate
  # minfdp=hatFDP(q,indReveal)

  # If all values in indReveal are TRUE, set qval to 1 for all entries
  if(all(indReveal)){
    cat("All p-values are revealed at begining, no fitting needed.\n")
    return(list(qval=rep(1,length(q))))
  }
  
  # Model selection if needed for coarse tuning 
  if(ms){
    ms1 <- seq(1,5901,100)
    models <- model
    model1 <- models[ms1]
    ms_res <- model_selection_cv(x,q,indReveal,model1,gamma)
    model <- ms_res$model
    # extract models with top 3 loglikelihood for further fine tuning
    modellogl <- order(ms_res$loglinfo,decreasing=T)[1:3]
    modellogl <- unique(c(sapply(ms1[modellogl], function(x){c((x-50):(x+50))})))
    modellogl <- modellogl[modellogl>0]
    cat("logl_res", ms_res$loglinfo, "\n")
  }

  # If initial parameters aren't provided, compute them
  # todo: compute initial parameters and G
  if(is.null(param0)){
    param=initialparameter(q,indReveal,x,model)
  }else{
    param=param0
  }

  # plus the initial fit, total nfits
  nBatch <- nfits
  batchSize <- floor((n-nInitialReveal)/nBatch)

  # data.frame for holding batch starting and ending position
  batchDF <- data.frame(start=rep(0, nBatch), end=nInitialReveal+(1:nBatch)*batchSize)
  batchDF$start <- batchDF$end-batchSize+1
  batchDF$end[nBatch] <- n

  for(i in 1:nBatch){
    cat("\nBatch", i, date(), "\n")

    Lfdr <- calcLfdr(param, q, gamma)
    # reveal the masked ones with the largest Lfdr
    Lfdr[indReveal] <- -1
    newReveal <- order(Lfdr, decreasing = TRUE)[1:(batchDF$end[i]-batchDF$start[i]+1)]
    revealOrder[batchDF$start[i]:batchDF$end[i]] <- newReveal
    indReveal[newReveal] <- TRUE
    if(ms){
      indMask <- !indReveal
      batchFDP <- hatFDP(sum(indMask & (abs(q) < s0)), sum(indMask & (abs(q) > 1-s0)))
      cat("\nBatchFDP", batchFDP, "\n")
      # model selection for fine tuning
      if(batchFDP<s0/2){
        model2 <- models[modellogl]
        ms_res <- model_selection_cv(x,q,indReveal,model2,gamma)
        model <- ms_res$model
        cat("logl_res", ms_res$loglinfo, "\n")
        ms=FALSE
      }
    }
    # refit the model
    if(i < nBatch) param <- modelfitting(param,x,q,indReveal,model,initial,gamma)
    # todo: pass the coefficient from the last batch to the next
    # if(i < nBatch){
    #   fit <- modelfitting(fit,x,q,indReveal,model)
    #   fit$param; fit$coefficients
    # }
  }

  # indicator whether q is in potential rejection region
  indRejRegion <- abs(q)>0.5
  # whether q is in rejection region along the reveal order
  indRej_RO <- indRejRegion[revealOrder[(nInitialReveal+1):n]]
  # Adec <- 1-indRej_RO
  # A <- (abs(q[revealOrder[(nInitialReveal+1):(n-1)]]) < 0.5)
  # R <- 1-A

  # FDP along the reveal order
  cumR <- cumsum(indRej_RO)
  cumA <- 1:length(cumR)-cumR
  FDP_RO <- cummin(c(initialFDP, hatFDP(initialA-cumA, initialR-cumR)))
  # remove the last one
  FDP_RO <- FDP_RO[-length(FDP_RO)]

  qval <- rep(1, n)
  # q-value for revealed and in rejection region
  qval[revealOrder[(nInitialReveal+1):n][indRej_RO]] <- FDP_RO[indRej_RO]
  return(qval)

}
