
logl_ms=function(fitinfo,x,q,indReveal,gamma){

  x = as.matrix(x)
  param11 <- predict(fitinfo[[1]], newdata = x)
  param12 <- predict(fitinfo[[2]], newdata = x)
  param1=cbind(param11,(1-param11)*param12,(1-param11)*(1-param12))

  param2 <- predict(fitinfo[[3]], newdata = x)
  param2 <- pmin(1/param2,1)

  param3 <- predict(fitinfo[[4]], newdata = x)
  param3 <- pmin(1/param3,1)

  param_new=cbind(param1,param2,param3)

  newlogl=logl(param_new,x,q,indReveal,gamma)

  return(newlogl)
}


