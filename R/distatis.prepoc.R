distatis.preproc <- function(Dist.mat, method = "distance", masses = NULL){ # start of function DiStatis.preproc
  # (1) preprocessing as MDS
  ## masses
  if (is.null(masses)){
    M_vec <- rep(1/nrow(Dist.mat), nrow(Dist.mat)) # sum of M_vec = 1
  }else{
    M_vec <- masses
  }
  ## centering matrix
  Cm <- diag(rep(1,dim(Dist.mat)[1])) - (as.matrix(rep(1,dim(Dist.mat)[1])) %*% t(M_vec))
  ## double-centering
  if (method == "covariance"){
    S <- 0.5*(Cm %*% Dist.mat %*% t(Cm))
  }else{
    S <- -0.5*(Cm %*% Dist.mat %*% t(Cm))
  }
  # (2) devided by the first eigen value
  Lambda <- eigen(S)$values[1]
  S.hat <- S/Lambda
  # output
  return(S.hat)
} # end of function DiStatis.preproc
