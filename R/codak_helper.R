pminmax=function(x,low=1e-15,up=1-1e-15) pmin(pmax(x,low),up)

pminmax2=function(x,low=1e-15,up=1-1e-15) ifelse(x<low,0,ifelse(x>up,1,x))

gen_Codak_model=function(pifun,alphafun,betafun,pifun_init,alphafun_init,betafun_init,pipre,alphapre,betapre,piargs=list(),alphaargs=list(),betaargs=list(),piargs_init=list(),alphaargs_init=list(),betaargs_init=list(),name=""){
  args=list(piargs=piargs,alphaargs=alphaargs,betaargs=betaargs,piargs_init=piargs_init,alphaargs_init=alphaargs_init,betaargs_init=betaargs_init)
  if(!is.function(pifun)){
    stop("pifun must be a function.")
  }
  if(!is.function(alphafun)){
    stop("alphafun must be a function.")
  }
  if(!is.function(betafun)){
    stop("betafun must be a function.")
  }
  if(!is.function(pifun_init)){
    stop("pifun_init must be a function.")
  }
  if(!is.function(alphafun_init)){
    stop("alphafun_init must be a function.")
  }
  if(!is.function(betafun_init)){
    stop("betafun_init must be a function.")
  }
  algo=list(pifun=pifun,alphafun=alphafun,betafun=betafun,pifun_init=pifun_init,alphafun_init=alphafun_init,betafun_init=betafun_init,pipre=pipre,alphapre=alphapre,betapre=betapre)
  model <- structure(
    list(name = name, algo = algo, args = args),
    class = "Codak_model"
  )
  return(model)
}
