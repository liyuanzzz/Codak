gen_model_glmnet_group=function(piargs=list(),alphaargs=list(),betaargs=list(),group,groupinfo=NULL){
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
    
    param111=list()
    for(ijk in unique(group)){
      if(is.null(groupinfo)){
        param111[[ijk]]=safe_glmnet(data[group==ijk,],temph1[group==ijk],family="binomial",weights=censorweights[group==ijk])$param
      }else{
        param111[[ijk]]=safe_glmnet(data[group==ijk,groupinfo[ijk,]],temph1[group==ijk],family="binomial",weights=censorweights[group==ijk])$param
      }
    }
    param11=rep(0,length(group))
    for(ijk in unique(group)){
      param11[group==ijk]=param111[[ijk]]
    }
    
    temph1=ifelse(temph==0,1,2-temph)
    
    param121=list()
    for(ijk in unique(group)){
      if(is.null(groupinfo)){
        param121[[ijk]]=safe_glmnet(data[group==ijk,],temph1[group==ijk],family="binomial",weights=ifelse(temph[group==ijk]==0,0,censorweights[group==ijk]))$param
      }else{
        param121[[ijk]]=safe_glmnet(data[group==ijk,groupinfo[ijk,]],temph1[group==ijk],family="binomial",weights=ifelse(temph[group==ijk]==0,0,censorweights[group==ijk]))$param
      }
    }
    param12=rep(0,length(group))
    for(ijk in unique(group)){
      param12[group==ijk]=param121[[ijk]]
    }
    
    param1=cbind(param11,(1-param11)*param12,(1-param11)*(1-param12))
    return(param1)
  }
  alphafun_init=function(q,indReveal,data,...){
    args=list(...)
    
    tempy=ifelse(indReveal,q,sign(q)*pmax(abs(q),1-abs(q)))
    tempy1=-log(1/2+pminmax(tempy,1e-15-1,1-1e-15)/2)
    
    weightt=ifelse(q<0, 1, 0)
    
    param21=list()
    for(ijk in unique(group)){
      if(is.null(groupinfo)){
        param21[[ijk]]=safe_glmnet(data[group==ijk,],tempy1[group==ijk],family="Gamma",weights=weightt[group==ijk])$param
      }else{
        param21[[ijk]]=safe_glmnet(data[group==ijk,groupinfo[ijk,]],tempy1[group==ijk],family="Gamma",weights=weightt[group==ijk])$param
      }
    }
    param2=rep(0,length(group))
    for(ijk in unique(group)){
      param2[group==ijk]=param21[[ijk]]
    }
    
    param2=pmin(1/param2,1)
    return(param2)
  }
  betafun_init=function(q,indReveal,data,...){
    args=list(...)
    
    tempy=ifelse(indReveal,q,sign(q)*pmax(abs(q),1-abs(q)))
    tempy2=-log(1/2-pminmax(tempy,1e-15-1,1-1e-15)/2)
    
    weightt=ifelse(q>0, 1, 0)
    
    param31=list()
    for(ijk in unique(group)){
      if(is.null(groupinfo)){
        param31[[ijk]]=safe_glmnet(data[group==ijk,],tempy2[group==ijk],family="Gamma",weights=weightt[group==ijk])$param
      }else{
        param31[[ijk]]=safe_glmnet(data[group==ijk,groupinfo[ijk,]],tempy2[group==ijk],family="Gamma",weights=weightt[group==ijk])$param
      }
    }
    param3=rep(0,length(group))
    for(ijk in unique(group)){
      param3[group==ijk]=param31[[ijk]]
    }
    
    param3=pmin(1/param3,1)
    return(param3)
  }
  pifun=function(data=NULL, weights=NULL, pre=NULL, init=NULL, ...){
    y=data[,1]
    x=as.matrix(data[,-1])
    param111=list()
    for(ijk in unique(group)){
      if(is.null(groupinfo)){
        param111[[ijk]]=safe_glmnet(x[c(group==ijk,group==ijk),],y[c(group==ijk,group==ijk)],family="binomial",weights=weights[c(group==ijk,group==ijk)])$param
      }else{
        param111[[ijk]]=safe_glmnet(x[c(group==ijk,group==ijk),groupinfo[ijk,]],y[c(group==ijk,group==ijk)],family="binomial",weights=weights[c(group==ijk,group==ijk)])$param
      }
    }
    param11=rep(0,length(group))
    for(ijk in unique(group)){
      param11[group==ijk]=param111[[ijk]][1:sum(group==ijk)]
    }
    return(list(param=param11))
  }
  alphabetafun=function(data=NULL, weights=NULL, pre=NULL, init=NULL, ...){
    y=data[,1]
    x=as.matrix(data[,-1])
    
    param21=list()
    for(ijk in unique(group)){
      if(is.null(groupinfo)){
        param21[[ijk]]=safe_glmnet(x[c(group==ijk,group==ijk),],y[c(group==ijk,group==ijk)],family="Gamma",weights=weights[c(group==ijk,group==ijk)])$param
      }else{
        param21[[ijk]]=safe_glmnet(x[c(group==ijk,group==ijk),groupinfo[ijk,]],y[c(group==ijk,group==ijk)],family="Gamma",weights=weights[c(group==ijk,group==ijk)])$param
      }
    }
    param2=rep(0,length(group))
    for(ijk in unique(group)){
      param2[group==ijk]=param21[[ijk]][1:sum(group==ijk)]
    }
    param2=pmin(1/param2,1)
    return(list(param=param2))
  }
  
  pipre=function(data=NULL, weights=NULL, ...){
    return(NULL)
  }
  alphabetapre=function(data=NULL, weights=NULL, ...){
    return(NULL)
  }
  

  piargs=c(list(group=group,groupinfo=groupinfo,type="glmnet"),piargs)
  piargs_init=piargs
  alphaargs=c(list(group=group,groupinfo=groupinfo,type="glmnet"),alphaargs)
  alphaargs_init=alphaargs
  betaargs=c(list(group=group,groupinfo=groupinfo,type="glmnet"),betaargs)
  betaargs_init=betaargs
  gen_Codak_model(pifun,alphabetafun,alphabetafun,pifun_init,alphafun_init,betafun_init,pipre,alphabetapre,alphabetapre,
                  piargs,alphaargs,betaargs,piargs_init,alphaargs_init,betaargs_init,name="glmnet")
}
