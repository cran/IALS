aPCAR<-function(W1,W2,k1,alpha){
  M=alpha*W1+W2
  return(svds(M,k1,k1,k1)$u)
}

aPCAC<-function(S1,S2,k2,alpha){
  M=alpha*S1+S2
  return(svds(M,k2,k2,k2)$u)
}

Orthogonalize=function(Z){
  gs=gramSchmidt(Z)
  Q <- gs$Q
  return(Q)
}

Distance=function(Z1,Z2){
  r1=ncol(Z1)
  r2=ncol(Z2)
  r=max(r1,r2)
  Q1=Orthogonalize(Z1)
  Q2=Orthogonalize(Z2)
  d=sqrt(1-1/r*sum(diag(Q1%*%t(Q1)%*%Q2%*%t(Q2))))
  return(d)
}
