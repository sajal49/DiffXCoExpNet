test_dis_coexpnet = function(nvar, sm){
  
  p = matrix(0, ncol=nvar, nrow=sm)
  for(i in 1:nvar){
    p[,i] = runif(sm, min=0, max=100);
  }
  
  c = matrix(0, ncol=nvar, nrow=sm)
  for(i in 1:nvar){
    p[,i] = runif(sm, min=0, max=100);
  }
  
  pnames = paste0("P",c(1:nvar))
  cnames = paste0("C",c(1:nvar))
  
  res = Discrete_Coexpnet(p, c, pnames, cnames, 10)
  
  return(res)
}