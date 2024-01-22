IALS <- function(X,W1=NULL,W2=NULL,m1,m2,max_iter=100,ep=1e-6){
  T=dim(X)[[1]]
  p1=dim(X)[[2]]
  p2=dim(X)[[3]]
  
  if(class(W1)[1]=="NULL"&class(W2)[1]=="NULL"){
    if (m1==1 & m2==1){
      W1=W2=matrix(1,p1,1)
    }
    else if (m1==1 & m2!=1){
      W1=matrix(1,p1,1)
      W2=matrix(1,p2,m2)
      for(i in 2:m2){
        w2i=c(rep(1,i-1),rep(-1,i-1))
        w2i_length=length(w2i)
        d1=p2%/%w2i_length
        d2=p2%%w2i_length
        if(d2==0){
          W2[,i]=rep(w2i,d1)
        } else{
          W2[,i]=c(rep(w2i,d1),w2i[1:d2])
        }
      }
    }
    else if (m1!=1 & m2==1){
      W1=matrix(1,p1,m1)
      W2=matrix(1,p2,1)
      for(i in 2:m1){
        w1i=c(rep(1,i-1),rep(-1,i-1))
        w1i_length=length(w1i)
        d1=p1%/%w1i_length
        d2=p1%%w1i_length
        if(d2==0){
          W1[,i]=rep(w1i,d1)
        } else{
          W1[,i]=c(rep(w1i,d1),w1i[1:d2])
        }
      }
    }else{
      W1=matrix(1,p1,m1)
      W2=matrix(1,p2,m2)
      for(i in 2:m1){
        w1i=c(rep(1,i-1),rep(-1,i-1))
        w1i_length=length(w1i)
        d1=p1%/%w1i_length
        d2=p1%%w1i_length
        if(d2==0){
          W1[,i]=rep(w1i,d1)
        } else{
          W1[,i]=c(rep(w1i,d1),w1i[1:d2])
        }
      }
      for(i in 2:m2){
        w2i=c(rep(1,i-1),rep(-1,i-1))
        w2i_length=length(w2i)
        d1=p2%/%w2i_length
        d2=p2%%w2i_length
        if(d2==0){
          W2[,i]=rep(w2i,d1)
        } else{
          W2[,i]=c(rep(w2i,d1),w2i[1:d2])
        }
      }
    }
  }else{
    W1=W1;W2=W2
  }
  
  Fhat=array(0,c(T,m1,m2))
  for(t in 1:T){
    Fhat[t,,]=t(W1)%*%X[t,,]%*%W2/(p1*p2)
  }
  
  C_old=W2
  R_old=W1
  Fhat_old=Fhat
  
  for(i in 1:max_iter){
    sumR=matrix(0,p1,m1)
    for(t in 1:T){
      sumR=sumR+X[t,,]%*%C_old%*%t(Fhat_old[t,,])
    }
    
    E <- eigen(t(sumR)%*%sumR)
    A = matrix(0,m1,m1)
    diag(A) <- 1/sqrt(E$values)
    R_new=sqrt(p1)*sumR%*%E$vector%*%A%*%t(E$vector)
    
    sumC=matrix(0,p2,m2)
    for(t in 1:T){
      sumC=sumC+t(X[t,,])%*%R_new%*%Fhat_old[t,,]
    }
    
    e <- eigen(t(sumC)%*%sumC)
    B = matrix(0,m2,m2)
    diag(B) <- 1/sqrt(e$values)
    C_new=sqrt(p2)*sumC%*%e$vector%*%B%*%t(e$vector)
    
    Fhat_new=array(0,c(T,m1,m2))
    for(t in 1:T){
      Fhat_new[t,,]=t(R_new)%*%X[t,,]%*%C_new/(p1*p2)
    }
    
    CC_old=array(0,c(T,p1,p2))
    CC_new=array(0,c(T,p1,p2))
    CCdiff=c()
    for(t in 1:T){
      CC_new[t,,]=R_new%*%Fhat_new[t,,]%*%t(C_new)
      CC_old[t,,]=R_old%*%Fhat_old[t,,]%*%t(C_old)
      CCdiff[t]=norm(CC_new[t,,]-CC_old[t,,],"F")
    }
    
    if(sum(CCdiff) <= ep*T*p1*p2){
      return(list(R=R_new,C=C_new,F=Fhat_new,iter=i))
    }else{
      C_old=C_new
      R_old=R_new
      Fhat_old=Fhat_new
    }
  }
  return(list(R=R_new,C=C_new,F=Fhat_new,iter=max_iter))
}

KIALS <- function(X,W1=NULL,W2=NULL,kmax,max_iter=100,ep=1e-6){
  fit=IALS(X,W1,W2,kmax,kmax,max_iter,ep)
  Ftilde=fit$F
  
  sumF_r=matrix(0,kmax,kmax)
  for (t in 1:T) {
    sumF_r=sumF_r+Ftilde[t,,]%*%t(Ftilde[t,,])
  }
  eigR=eigen(sumF_r%*%t(sumF_r))$values[1:(1+kmax)]
  k1=which.max(eigR[1:kmax]/(eigR[2:(1+kmax)]))
  
  sumF_c=matrix(0,kmax,kmax)
  for (t in 1:T) {
    sumF_c=sumF_c+t(Ftilde[t,,])%*%Ftilde[t,,]
  }
  eigC=eigen(sumF_c%*%t(sumF_c))$values[1:(1+kmax)]
  k2=which.max(eigC[1:kmax]/(eigC[2:(1+kmax)]))
  
  return(list(k1=k1,k2=k1))
}
