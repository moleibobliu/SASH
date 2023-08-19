seed <- 1
iter <- 1
n <- 500
N <- 1000
d <- 10
M <- 2
setting <- "strong"
magnitude <- 1


set.seed(seed)


if(setting == "weak"){
  para1 <- 3
  para2 <- 0.75
} else if(setting == "strong"){
  para1 <- 5
  para2 <- 0.85
}

cat(para1, para2, "\n")


X_lst <- vector('list', M)
Y_lst <- vector('list', M)
S_lst <- vector('list', M)
N_lst <- rep(N, M)


######### generate X #########

### Pois  
X <- NULL
for(m in 1:M){
  mu <- rep(5, d)
  r <- 0.25
  Sigma <- matrix(r, d, d) + diag(rep(1 - r, d))
  
  X_lst[[m]] <- rmvpois(N, mu, Sigma)
  
  X_lst[[m]] <- log(X_lst[[m]] + 1)
  
  X <- rbind(X, X_lst[[m]])
}


######### define beta #########

alpha_star <- 0
beta_star <- c(magnitude*1, -magnitude,
               rep(magnitude/2, 1), rep(-magnitude/2, 1),
               rep(magnitude/4, 1), rep(-magnitude/4, 1),
               rep(magnitude/8, 1), rep(-magnitude/8, 1))

beta_star <- matrix(c(beta_star, rep(0, d - length(beta_star))))
beta_star

prob <- NULL
for(m in 1:M){
  prob <- c(prob, logit(alpha_star + X_lst[[m]] %*% beta_star))
}

Y <- rbinom(sum(N_lst), size = 1, prob)
# size:	number of trials (zero or more).


######### generate S #########

S <- NULL
for(m in 1:(M/2)){
  Y_m <- Y_lst[[m]] <- Y[((m-1)*N + 1):(m*N)]
  mu_1 <- para1
  mu_0 <- exp(0)
  S <- c(S, rpois(N, lambda = (mu_1)^Y_m * (mu_0)^(1-Y_m)))
}
S <- log(S + 1)

for(m in (M/2 + 1):M){
  Y_m <- Y_lst[[m]] <- Y[((m-1)*N + 1):(m*N)]
  mu_1 <- para2
  mu_0 <- 1 - para2
  print(mu_1)
  print(mu_0)
  S <- c(S, rbinom(N, 1, (mu_1)^Y_m * (mu_0)^(1-Y_m)))
}

for(m in 1:M){
  S_lst[[m]] <- S[((m-1)*N + 1):(m*N)]
}

X_label <- X[1:n, ]
Y_label <- Y[1:n]






######### SASH #########


SASH_res = SASH(X_label,Y_label,X_lst,S_lst)


######### IPD #########

IPD_res = IPD(X_label,Y_label,X_lst,S_lst)




else if(method == 2){
  #### 2 METHOD ######
  #################### comparing - no Y ###################
  
  final=NULL
  Grad=Hess=NULL
  beta_debias=NULL
  variance=NULL
  
  # loss_step3=function(b,lambda_step3){
  #   # b: length d-1
  #   b=matrix(b)
  #   sum=0
  #   for(m in 1:M){
  #     sum=sum+t(rbind(1,b))%*%Hess[[m]]%*%rbind(1,b)-2*t(rbind(1,b))%*%Grad[[m]]
  #   }  
  #   # print(sum/M)
  #   # print(lambda_step3*sum(abs(b)))
  #   return(sum/M+lambda_step3*sum(abs(b)))
  #   
  # }
  ### step 1 #### 
  lambda_final=sigma.square_final=NULL
  for(m in 1:M){
    cat("m =",m,"\n")
    X_m=X_lst[[m]]
    S_m=S_lst[[m]]
    
    cat("iter",iter,"\n")
    
    
    beta=c(1,rep(0,d-1))
    
    beta_pool=NULL
    
    for (t in 1:iter){
      cat("t =",t,"\n")
      # if(t<4){
      #   h=(log(d)/n)^(1/4)
      # }else{
      #   h=(log(d)/N)^(1/4)
      # }
      BIC_step2=NULL
      beta_step2=NULL
      # h_lst=c((5)*0.2*((log(d)/N_lst[m])^(1/4)))
      h_lst=c((maxh)*sqrt(log(N_lst[m])/N_lst[m]))
      
      if (t==1){
        h=h_lst[1]
        print(Sys.time())
        A=apply(X_lst[[m]],1,function(o){f_hat(X_lst[[m]],o,beta,S_lst[[m]],h)})
        
        print(rbind(A[1:20],S_lst[[m]][1:20],Y_lst[[m]][1:20]))
        print( mean(A[Y_lst[[m]]==1],na.rm = T))
        print(mean(A[Y_lst[[m]]==0],na.rm = T))
        print(mean(A,na.rm = T))
        print(mean(S_lst[[m]],na.rm = T))
        
        
        B=apply(X_lst[[m]],1,function(o){f_hat_grad(X_lst[[m]],o,beta,S_lst[[m]],h)})
        # S_m-A-c(t(rbind(1,b)-beta)%*%B)
        # w=1/A
        if(m<=(M/2)){
          w=rep(1,length(A))
        }else{
          w=1/(A*(1-A))
        }
        w[abs(w)>100]=100
        sigma.square=mean(w*(S_lst[[m]]-A)^2)
        print("sigma.square")
        print(sigma.square)
        
        
        sigma.square_final=c(sigma.square_final,sigma.square)
      }
      
      for(j in 1:length(h_lst)){
        
        h=h_lst[j]
        print(Sys.time())
        A=apply(X_m,1,function(o){f_hat(X_m,o,beta,S_m,h)})
        
        print(rbind(A[1:20],S_lst[[m]][1:20],Y_lst[[m]][1:20]))
        print( mean(A[Y_lst[[m]]==1],na.rm = T))
        print(mean(A[Y_lst[[m]]==0],na.rm = T))
        
        B=apply(X_m,1,function(o){f_hat_grad(X_m,o,beta,S_m,h)})
        # S_m-A-c(t(rbind(1,b)-beta)%*%B)
        # w=1/A
        if(m<=(M/2)){
          w=rep(1,length(A))
        }else{
          w=1/(A*(1-A))
        }
        w[abs(w)>100]=100
        cv.result=cv.glmnet(t(B[-1,]),c(S_m-A+t(beta)%*%B-B[1,]), 
                            weights=w,family="gaussian",
                            lambda=0.005*(1:120)*sqrt(log(d)/N),
                            alpha=1,intercept = FALSE)
        # plot(cv.result)
        cat("lambda",cv.result$lambda.min/sqrt(log(d)/N),"\n")
        beta_new=coef(cv.result, s="lambda.min")[-1]
        
        
        # beta=c(1,beta_new)
        # print(Sys.time())
        print(round(c(1,beta_new),3)[1:15])
        
        
        print(Sys.time())
        
        BIC_step2=c(BIC_step2,
                    loss_step2(beta_new,
                               c(S_lst[[m]]-A+t(beta)%*%B-B[1,]),t(B[-1,]),
                               lamb=cv.result$lambda.min/sqrt(log(d)/N)
                               ,N_lst[[m]],w,sigma.square))
        beta_step2=cbind(beta_step2,beta_new)
        
      } # end for h
      print(round(cbind(beta_step2,
                        beta_star[-1]/beta_star[1]),3))
      beta_h=beta_step2[,which.min(BIC_step2)]
      print(which.min(BIC_step2))
      
      beta=c(1,beta_h)
      
      # loss(gamma_initial[-1])
      # loss(beta[-1])
      # loss(beta_star[-1]/beta_star[2])
      
    } # end iter
    # print(round(beta_pool,4))
    
    Grad[[m]]=grad(X_m,S_m,beta,h,m,N)
    Hess[[m]]=hess(X_m,S_m,beta,h,m,N)
    
    final=cbind(final,beta)
    
    cat("\n")
    
    
    
  } # end m
  
  rownames(final)=1:d
  ind=apply(abs(final),1,sum)>0.00001
  print(cbind(round(final[ind,],3), beta_star[ind]/beta_star[1]))
  
  
  dim(Grad[[1]])
  dim(Hess[[1]])
  #print(Grad)
  #print(Hess)
  
  Grad_pool=Hess_pool=0
  for(m in 1:M){
    Grad_pool=Grad_pool+Grad[[m]]
    Hess_pool=Hess_pool+Hess[[m]]
  }  
  
  
  
  
  ##### noY1 #####
  
  vec=eigen(Hess_pool[-1,-1])$vectors
  
  value=eigen(Hess_pool[-1,-1])$values
  value[abs(value)<1e-10]=0
  
  XX=diag(sqrt(value))%*%t(vec)
  YY=solve(t(XX))%*%(Grad_pool[-1]-Hess_pool[-1,1])
  
  cv.result=cv.glmnet(XX,YY,family="gaussian",
                      alpha=1,intercept = FALSE)
  # plot(cv.result)
  cat("\n","lambda",cv.result$lambda.min/sqrt(log(d)/N),"\n")
  
  model=glmnet(XX,YY,family="gaussian", 
               lambda =cv.result$lambda.min,intercept = FALSE)
  print(as.vector(round(model$beta,4))[1:20])
  ## our1
  beta_dagger=NULL
  beta_dagger=rbind(beta_dagger,c(1,as.vector(model$beta)))
  
  cat(c(seed,iter,n,N,d,M,"noY1",c(0,beta_dagger[1,])),"\n",file=file,append=TRUE)
  
  
  ### BIC ####
  BIC=NULL
  beta_step3=NULL
  for(lambda in 0.0025*(1:600)*sqrt(log(d)/N) ){
    model=glmnet(XX,YY,family="gaussian", 
                 lambda =lambda,intercept = FALSE)
    # print(as.vector(round(model$beta,4)))
    bb=model$beta
    beta_step3=cbind(beta_step3,bb)
    
    
    
    if(sum(abs(bb)>1e-8)>0){
      cat(M*N*loss_step3(model$beta,
                         Grad_pool,Hess_pool, lambda,
                         0, N_lst),
          " ")
      cat((length(bb)-sum(abs(bb)<1e-8))*log(M*N),"\n")
    }else{
      break
    }
    
    BIC=c(BIC,(length(bb)-sum(abs(bb)<1e-8))*log(M*N)+
            M*N*loss_step3(model$beta,
                           Grad_pool,Hess_pool, lambda,
                           0, N_lst))
    
  }
  print(which.min(BIC))
  print(c(beta_step3[,which.min(BIC)])[1:25])
  ## ## our2
  
  beta_dagger=rbind(beta_dagger,c(1,beta_step3[,which.min(BIC)]))
  
  print(beta_dagger)
  # print(BIC_dagger)
  
  
  cat(c(seed,iter,n,N,d,M,"noY2",c(0,beta_dagger[2,])),"\n",file=file,append=TRUE)
  
  # for (m in 1:M){
  #   fit_last=glm(Y_label~X_label%*%final[,m],family="binomial")
  #   print(fit_last[[1]])
  #   print(c(fit_last[[1]][1],fit_last[[1]][2]*final[,m])[1:20])
  #   
  #   cat(c(seed,iter,n,N,d,M,paste("noy_Site",m,sep=""),c(fit_last[[1]][1],fit_last[[1]][2]*final[,m])),"\n",file=file,append=TRUE)
  # }
  
}else if(method == 3){
  #### 3 METHOD #####
  
  
  cvfit=cv.glmnet(X_label,Y_label,
                  # lambda=seq(0.05,1.5,0.01)*sqrt(log(d)/n),
                  family="binomial",alpha=1)
  # plot(cvfit)
  # sqrt(log(d)/n)
  cvfit$lambda.min/sqrt(log(d)/n)
  # cvfit$lambda.1se
  print(coef(cvfit, s="lambda.min")[1:20])
  print(coef(cvfit, s="lambda.min")[2:20]/coef(cvfit, s="lambda.min")[2])
  
  if(coef(cvfit, s="lambda.min")[2] != 0) {
    index0 = 1
  } else {
    index0=order(abs(coef(cvfit, s="lambda.min"))[2:3],
                 decreasing = T)[1]
  }
  index0
  index=1:d
  index[index0]=1
  index[1]=index0
  index[1:10]
  
  beta_initial=coef(cvfit, s="lambda.min")[-1][index]
  
  gamma_initial=matrix(beta_initial/ifelse(beta_initial[1]==0,1,beta_initial[1]))
  gamma_initial[1]=1
  c(gamma_initial)[1:10]
  dim(gamma_initial)
  
  beta_initial = c(coef(cvfit, s="lambda.min")[1], beta_initial)
  c(beta_initial)[1:10]
  cat(c(seed,iter,n,N,d,M,"index",c(0,index)),"\n",file=file,append=TRUE)
  
  
  X_label=X_label[,index]
  for(m in 1:M){
    X_lst[[m]]=X_lst[[m]][,index]
    X_lst[[m]]=as.matrix(X_lst[[m]])
  }
  
  X= X[,index]
  ######## no sharing constraint ###################
  
  beta=c(gamma_initial)
  beta_pool=NULL
  cat("iter",iter,"\n")
  for (t in 1:(iter)){
    cat("t =",t,"\n")
    
    
    
    print(Sys.time())
    w=AA=BB=NULL
    lambda_final=sigma.square_final=NULL
    for (m in 1:M){
      cat("m=",m, "\n") 
      
      h=c((maxh)*sqrt(log(N_lst[m])/N_lst[m]))
      
      X_m=X_lst[[m]]
      S_m=S_lst[[m]]
      A=apply(X_m,1,function(o){f_hat(X_m,o,beta,S_m,h)})
      B=apply(X_m,1,function(o){f_hat_grad(X_m,o,beta,S_m,h)})
      AA=c(AA,A)
      BB=cbind(BB,B)
      if(m<=(M/2)){
        ww=rep(1,length(A))
      }else{
        ww=(1/(A*(1-A)))
      }
      w=c(w,ww)
      w[abs(w)>100]=100
      print(Sys.time())
      
      sigma.square=mean(ww*(S_m-A)^2)
      sigma.square_final=c(sigma.square_final,
                           sigma.square)
      
      print(length(ww))
      print(length(w))
      print(length(AA))
      print(dim(BB))
    }
    
    print(length(w))
    print(length(AA))
    print(dim(BB))
    
    cat("***","\n")
    
    XX=diag(sqrt(w))%*%t(BB[-1,])
    dim(XX)
    cat("****","\n")
    YY = diag(sqrt(w))%*%c(S-AA+t(beta)%*%BB-BB[1,])
    dim(YY)
    cat("*****","\n")
    print(Sys.time())
    
    loss_IPD=function(b,YYY,XXX,lamb,N_lst,weight,sigma.square_final){
      # b: length d-1
      b=matrix(b)
      pro=weight*(YYY-XXX%*%b)^2
      result_IPD=NULL
      for(m in 1:length(N_lst)){
        avg=mean(pro[((m-1)*N+1):(m*N)])
        print(avg)
        result_IPD=c(result_IPD,
                     N_lst[m]*(avg)/sigma.square_final[m])
      }
      
      # print(sum/M)
      # print(lambda_step3*sum(abs(b)))
      # if((length(b)-sum(abs(b)<1e-8))>0){
      #   cat(Nm*(avg)/sigma.square," ",(length(b)-sum(abs(b)<1e-8))*log(Nm),"\n")
      # }
      
      return(sum(result_IPD)+
               (length(b)-sum(abs(b)<1e-8))*log(sum(N_lst)))
      
    }
    
    
    BIC_step2=NULL
    beta_step2=NULL
    lambda_lst=0.005*(1:400)*sqrt(log(d)/N)
    B1=B2=NULL
    for (lambda in lambda_lst){
      
      cv.result=glmnet(t(BB[-1,]),c(S-AA+t(beta)%*%BB-BB[1,]),
                       weights=w,family="gaussian",
                       lambda=lambda,
                       alpha=1,intercept = FALSE)
      # plot(cv.result)
      
      beta_new=coef(cv.result, s="lambda.min")[-1]
      # 
      # avg=mean(w*(c(S_lst[[m]]-A+t(beta)%*%B-B[1,])-t(B[-1,])%*%beta_new)^2)
      # # B1=c(B1,N_lst[[m]]*(avg)/sigma.square)
      # B2=c(B2,(length(beta_new)-sum(abs(beta_new)<1e-8))*log(N_lst[[m]]))
      # 
      BIC_step2=c(BIC_step2,
                  loss_IPD(beta_new,
                           c(S-AA+t(beta)%*%BB-BB[1,]),
                           t(BB[-1,]),
                           lamb=lambda
                           ,N_lst,w,sigma.square_final))
      beta_step2=cbind(beta_step2,beta_new)
      
    } # end for lambda
    # print(round(cbind(beta_step2,
    #                   gamma_initial[-1],
    #                   beta_star[-1]/beta_star[1]),3))
    print(head(beta_step2[,which.min(BIC_step2)],15))
    print(which.min(BIC_step2))
    print("lambda")
    IPD_lambda=lambda_lst[which.min(BIC_step2)]
    print(lambda_lst[which.min(BIC_step2)]/sqrt(log(d)/N))
    beta=c(1,beta_step2[,which.min(BIC_step2)])
    
    # print(Sys.time())
    print(round(c(beta),3)[1:15])
    
    beta_pool=cbind(beta_pool,beta)
    
    
    cat("t=",t,"\n")
    print(c(beta_pool))
    # beta
    fit_last=glm(Y_label~X_label%*%beta,family="binomial")
    print(fit_last[[1]])
    print(round(c(fit_last[[1]][1],fit_last[[1]][2]*beta),4))
    
    cat(c(seed,iter,n,N,d,M,paste0("POOLED",t),
          c(fit_last[[1]][1],fit_last[[1]][2]*beta[index])),"\n",file=file,append=TRUE)
    
    POOLED=c(fit_last[[1]][1],fit_last[[1]][2]*beta)
    
  } # end iter
  dim(XX)
  dim(YY)
  
  
  Grad=Hess=NULL
  
  for (m in 1:M){
    Grad[[m]]=grad(X_lst[[m]],S_lst[[m]],beta,h,m,N)
    Hess[[m]]=hess(X_lst[[m]],S_lst[[m]],beta,h,m,N)
    
  } # end m
  
  
  dim(Grad[[1]])
  dim(Hess[[1]])
  
  
  Grad_pool=Hess_pool=0
  for(m in 1:M){
    print(median(diag(Hess[[m]]))*N)
    Grad_pool=Grad_pool+Grad[[m]]*N
    Hess_pool=Hess_pool+Hess[[m]]*N
  }
  
  Hess_pool[1:5,1:5]
  eigen(Hess_pool)$values
  
  
  #####
  N_lst = rep(N,M)
  # loss_step3=function(b,Grad_final,Hess_final,lambda){
  #   # b: length d-1
  #   b=matrix(b)
  #   sum=t(rbind(1,b))%*%Hess_final%*%rbind(1,b)-
  #     2*t(rbind(1,b))%*%Grad_final
  #   # print(sum/M)
  #   # print(lambda_step3*sum(abs(b)))
  #   return(sum/(sum(N_lst)+n)+lambda/(sum(N_lst)+n)*sum(abs(b)))
  #   
  # }
  
  
  fit=glm(Y_label~X_label%*%beta,
          #family="binomial",
          family="binomial")
  print(fit[[1]])
  
  beta_tilde=beta
  for(t in 1:(iter)){
    a_tilde=fit[[1]][1]
    b_tilde=fit[[1]][2]
    cat(a_tilde,b_tilde,"\n")
    H_tilde=0
    for(j in 1:n){
      o=X_label[j,]
      H_tilde=H_tilde+
        c(b_tilde^2*g_grad(a_tilde+b_tilde*(o%*%beta_tilde)))*(o%*%t(o))
    }
    H_tilde[1:6,1:6]
    eigen(H_tilde)$values
    H_tilde=0.5*H_tilde
    
    xi_tilde=2*H_tilde%*%beta_tilde
    for(j in 1:n){
      o=X_label[j,]
      xi_tilde=xi_tilde+
        c(Y_label[j]-g(a_tilde+b_tilde*(o%*%beta_tilde)))*b_tilde*o
    }
    
    Hess_final=Hess_pool+H_tilde
    Grad_final=Grad_pool+0.5*xi_tilde
    
    vec=eigen(Hess_final[-1,-1])$vectors
    
    value=eigen(Hess_final[-1,-1])$values
    value[abs(value)<1e-10]=0
    
    XX=diag(sqrt(value))%*%t(vec)
    YY=solve(t(XX))%*%(Grad_final[-1]-Hess_final[-1,1])
    
    BIC=NULL
    beta_step3=NULL
    for(lambda in 0.003*(1:2000)*sqrt(log(d)/N) ){
      model=glmnet(XX,YY,family="gaussian",
                   lambda =lambda,intercept = FALSE)
      # print(as.vector(round(model$beta,4)))
      bb=model$beta
      bb[1:15]
      beta_step3=cbind(beta_step3,bb)
      # print(length(bb)-sum(abs(bb)<1e-8))
      
      if(sum(abs(bb)>1e-8)>0){
        cat((n+sum(N_lst))*loss_step3(model$beta,
                                      Grad_final,Hess_final,lambda,n, N_lst),
            " ")
        cat((length(bb)-sum(abs(bb)<1e-8))*log(n+sum(N_lst)),"\n")
        
      }
      
      
      BIC=c(BIC,
            (length(bb)-sum(abs(bb)<1e-8))*log(n+sum(N_lst))+
              (n+sum(N_lst))*(loss_step3(model$beta,
                                         Grad_final,Hess_final,lambda, n, N_lst)))
      
    }
    print(which.min(BIC))
    print(round(c(beta_step3[,which.min(BIC)]),4)[1:30])
    
    
    beta_tilde=c(1,c(beta_step3[,which.min(BIC)]))
    
    fit=glm(Y_label~X_label%*%beta_tilde,
            family="binomial")
    print(fit[[1]])
    
  }
  
  
  
  a_final=fit[[1]][1]
  b_final=fit[[1]][2]
  beta_dir_final=beta_tilde
  matrix=NULL
  matrix=rbind(matrix,c(a_final,b_final*beta_dir_final))
  
  beta_dagger=NULL
  beta_dagger=rbind(beta_dagger,beta_dir_final)
  
  cat(c(seed,iter,n,N,d,M,"IPD",
        c(matrix[1,1],matrix[1,-1][index])),"\n",file=file,append=TRUE)
  
  ##### debias #####
  deb_res = debias2(beta_fit =  t(matrix), 
                    beta_star,
                    X_label = X_label,
                    n_lst = n,  tau_lst = 0.001 * 1:1000 * sqrt((log(d)) / n))  
  
  as.numeric(deb_res[[3]])
  
  cat(c(seed,iter,n,N,d,M,"IPD-test", 
        as.numeric(c(deb_res[[3]][1], deb_res[[3]][-1][index]))),
      "\n",file=file,append=TRUE)
  
  cat(c(seed,iter,n,N,d,M,"IPD-cover", 
        as.numeric(c(deb_res[[4]][1], deb_res[[4]][-1][index]))),
      "\n",file=file,append=TRUE)
  
  cat(c(seed,iter,n,N,d,M,"IPD-deb", 
        as.numeric(c(deb_res[[1]][1], deb_res[[1]][-1][index]))),
      "\n",file=file,append=TRUE)
  
  cat(c(seed,iter,n,N,d,M,"IPD-var", 
        as.numeric(c(deb_res[[2]][1], deb_res[[2]][-1][index]))),
      "\n",file=file,append=TRUE)
  
  print("done debias!")
  
  
  # end method 3  
}else if(method == 4){
  
  #### 4 METHOD #####
  
  cvfit=cv.glmnet(X_label,Y_label,
                  # lambda=seq(0.05,1.5,0.01)*sqrt(log(d)/n),
                  family="binomial",alpha=1)
  # plot(cvfit)
  # sqrt(log(d)/n)
  cvfit$lambda.min/sqrt(log(d)/n)
  # cvfit$lambda.1se
  print(coef(cvfit, s="lambda.min")[1:20])
  
  if(coef(cvfit, s="lambda.min")[2] != 0) {
    index0 = 1
  } else {
    index0=order(abs(coef(cvfit, s="lambda.min"))[2:3],
                 decreasing = T)[1]
  }
  
  index0
  index=1:d
  index[index0]=1
  index[1]=index0
  print(index[1:11])
  
  beta_initial=coef(cvfit, s="lambda.min")[-1][index]
  
  gamma_initial=matrix(beta_initial/ifelse(beta_initial[1]==0,1,beta_initial[1]))
  gamma_initial[1]=1
  c(gamma_initial)[1:10]
  dim(gamma_initial)
  
  beta_initial = c(coef(cvfit, s="lambda.min")[1], beta_initial)
  c(beta_initial)[1:10]
  
  cat(c(seed,iter,n,N,d,M,"index",c(0,index)),"\n",file=file,append=TRUE)
  
  
  X_label=X_label[,index]
  for(m in 1:M){
    X_lst[[m]]=X_lst[[m]][,index]
    X_lst[[m]]=as.matrix(X_lst[[m]])
  }
  
  X= X[,index]
  
  beta=c(gamma_initial)
  beta_pool=NULL
  cat("iter",iter,"\n")
  for (t in 1:(iter)){
    cat("t =",t,"\n")
    
    
    
    print(Sys.time())
    w=AA=BB=NULL
    lambda_final=sigma.square_final=NULL
    for (m in 1:M){
      cat("m=",m, "\n") 
      
      h=c((maxh)*sqrt(log(N_lst[m])/N_lst[m]))
      
      X_m=X_lst[[m]]
      S_m=S_lst[[m]]
      A=apply(X_m,1,function(o){f_hat(X_m,o,beta,S_m,h)})
      B=apply(X_m,1,function(o){f_hat_grad(X_m,o,beta,S_m,h)})
      AA=c(AA,A)
      BB=cbind(BB,B)
      if(m<=(M/2)){
        ww=rep(1,length(A))
      }else{
        ww=(1/(A*(1-A)))
      }
      w=c(w,ww)
      w[abs(w)>100]=100
      print(Sys.time())
      
      sigma.square=mean(ww*(S_m-A)^2)
      sigma.square_final=c(sigma.square_final,
                           sigma.square)
      
      print(length(ww))
      print(length(w))
      print(length(AA))
      print(dim(BB))
    }
    
    print(length(w))
    print(length(AA))
    print(dim(BB))
    
    cat("***","\n")
    
    XX=diag(sqrt(w))%*%t(BB[-1,])
    dim(XX)
    cat("****","\n")
    YY = diag(sqrt(w))%*%c(S-AA+t(beta)%*%BB-BB[1,])
    dim(YY)
    cat("*****","\n")
    print(Sys.time())
    
    loss_IPD=function(b,YYY,XXX,lamb,N_lst,weight,sigma.square_final){
      # b: length d-1
      b=matrix(b)
      pro=weight*(YYY-XXX%*%b)^2
      result_IPD=NULL
      for(m in 1:length(N_lst)){
        avg=mean(pro[((m-1)*N+1):(m*N)])
        print(avg)
        result_IPD=c(result_IPD,
                     N_lst[m]*(avg)/sigma.square_final[m])
      }
      
      # print(sum/M)
      # print(lambda_step3*sum(abs(b)))
      # if((length(b)-sum(abs(b)<1e-8))>0){
      #   cat(Nm*(avg)/sigma.square," ",(length(b)-sum(abs(b)<1e-8))*log(Nm),"\n")
      # }
      
      return(sum(result_IPD)+
               (length(b)-sum(abs(b)<1e-8))*log(sum(N_lst)))
      
    }
    
    
    BIC_step2=NULL
    beta_step2=NULL
    lambda_lst=0.005*(1:400)*sqrt(log(d)/N)
    B1=B2=NULL
    for (lambda in lambda_lst){
      
      cv.result=glmnet(t(BB[-1,]),c(S-AA+t(beta)%*%BB-BB[1,]),
                       weights=w,family="gaussian",
                       lambda=lambda,
                       alpha=1,intercept = FALSE)
      # plot(cv.result)
      
      beta_new=coef(cv.result, s="lambda.min")[-1]
      # 
      # avg=mean(w*(c(S_lst[[m]]-A+t(beta)%*%B-B[1,])-t(B[-1,])%*%beta_new)^2)
      # # B1=c(B1,N_lst[[m]]*(avg)/sigma.square)
      # B2=c(B2,(length(beta_new)-sum(abs(beta_new)<1e-8))*log(N_lst[[m]]))
      # 
      BIC_step2=c(BIC_step2,
                  loss_IPD(beta_new,
                           c(S-AA+t(beta)%*%BB-BB[1,]),
                           t(BB[-1,]),
                           lamb=lambda
                           ,N_lst,w,sigma.square_final))
      beta_step2=cbind(beta_step2,beta_new)
      
    } # end for lambda
    # print(round(cbind(beta_step2,
    #                   gamma_initial[-1],
    #                   beta_star[-1]/beta_star[1]),3))
    print(head(beta_step2[,which.min(BIC_step2)],15))
    print(which.min(BIC_step2))
    print("lambda")
    IPD_lambda=lambda_lst[which.min(BIC_step2)]
    print(lambda_lst[which.min(BIC_step2)]/sqrt(log(d)/N))
    beta=c(1,beta_step2[,which.min(BIC_step2)])
    
    # print(Sys.time())
    print(round(c(beta),3)[1:15])
    
    beta_pool=cbind(beta_pool,beta)
    
    
    cat("t=",t,"\n")
    print(c(beta_pool))
    # beta
    fit_last=glm(Y_label~X_label%*%beta,family="binomial")
    print(fit_last[[1]])
    print(round(c(fit_last[[1]][1],fit_last[[1]][2]*beta),4))
    
    # cat(c(seed,iter,n,N,d,M,paste0("POOLED",t),
    #       c(fit_last[[1]][1],fit_last[[1]][2]*beta[index])),"\n",file=file,append=TRUE)
    # 
    POOLED=c(fit_last[[1]][1],fit_last[[1]][2]*beta)
    
  } # end iter
  dim(XX)
  dim(YY)
  
  
  Grad=Hess=NULL
  
  for (m in 1:M){
    Grad[[m]]=grad(X_lst[[m]],S_lst[[m]],beta,h,m,N)
    Hess[[m]]=hess(X_lst[[m]],S_lst[[m]],beta,h,m,N)
    
  } # end m
  
  
  dim(Grad[[1]])
  dim(Hess[[1]])
  
  
  Grad_pool=Hess_pool=0
  for(m in 1:M){
    print(median(diag(Hess[[m]]))*N)
    Grad_pool=Grad_pool+Grad[[m]]*N
    Hess_pool=Hess_pool+Hess[[m]]*N
  }
  
  Hess_pool[1:5,1:5]
  eigen(Hess_pool)$values
  
  
  #####
  N_lst = rep(N,M)
  loss_step3=function(b,Grad_final,Hess_final,lambda){
    # b: length d-1
    b=matrix(b)
    sum=t(rbind(1,b))%*%Hess_final%*%rbind(1,b)-
      2*t(rbind(1,b))%*%Grad_final
    # print(sum/M)
    # print(lambda_step3*sum(abs(b)))
    return(sum/(sum(N_lst)+n)+lambda/(sum(N_lst)+n)*sum(abs(b)))
    
  }
  
  
  fit=glm(Y_label~X_label%*%beta,
          #family="binomial",
          family="binomial")
  print(fit[[1]])
  
  beta_tilde=beta
  for(t in 1:(2*iter)){
    a_tilde=fit[[1]][1]
    b_tilde=fit[[1]][2]
    cat(a_tilde,b_tilde,"\n")
    H_tilde=0
    for(j in 1:n){
      o=X_label[j,]
      H_tilde=H_tilde+
        c(b_tilde^2*g_grad(a_tilde+b_tilde*(o%*%beta_tilde)))*(o%*%t(o))
    }
    H_tilde[1:6,1:6]
    eigen(H_tilde)$values
    H_tilde=0.5*H_tilde
    
    xi_tilde=2*H_tilde%*%beta_tilde
    for(j in 1:n){
      o=X_label[j,]
      xi_tilde=xi_tilde+
        c(Y_label[j]-g(a_tilde+b_tilde*(o%*%beta_tilde)))*b_tilde*o
    }
    
    Hess_final=Hess_pool+H_tilde
    Grad_final=Grad_pool+0.5*xi_tilde
    
    vec=eigen(Hess_final[-1,-1])$vectors
    
    value=eigen(Hess_final[-1,-1])$values
    value[abs(value)<1e-10]=0
    
    XX=diag(sqrt(value))%*%t(vec)
    YY=solve(t(XX))%*%(Grad_final[-1]-Hess_final[-1,1])
    
    BIC=NULL
    beta_step3=NULL
    for(lambda in 0.003*(1:2000)*sqrt(log(d)/N) ){
      model=glmnet(XX,YY,family="gaussian",
                   lambda =lambda,intercept = FALSE)
      # print(as.vector(round(model$beta,4)))
      bb=model$beta
      bb[1:15]
      beta_step3=cbind(beta_step3,bb)
      # print(length(bb)-sum(abs(bb)<1e-8))
      
      if(sum(abs(bb)>1e-8)>0){
        cat((n+sum(N_lst))*loss_step3(model$beta,
                                      Grad_final,Hess_final,lambda),
            " ")
        cat((length(bb)-sum(abs(bb)<1e-8))*log(n+sum(N_lst)),"\n")
        
      }
      
      
      BIC=c(BIC,
            (length(bb)-sum(abs(bb)<1e-8))*log(n+sum(N_lst))+
              (n+sum(N_lst))*(loss_step3(model$beta,
                                         Grad_final,Hess_final,lambda)))
      
    }
    print(which.min(BIC))
    print(round(c(beta_step3[,which.min(BIC)]),4)[1:30])
    
    
    beta_tilde=c(1,c(beta_step3[,which.min(BIC)]))
    
    fit=glm(Y_label~X_label%*%beta_tilde,
            family="binomial")
    print(fit[[1]])
    
  }
  
  
  
  a_final=fit[[1]][1]
  b_final=fit[[1]][2]
  beta_dir_final=beta_tilde
  matrix=NULL
  matrix=rbind(matrix,c(a_final,b_final*beta_dir_final))
  
  beta_dagger=NULL
  beta_dagger=rbind(beta_dagger,beta_dir_final)
  
  cat(c(seed,iter,n,N,d,M,paste0("POOLED2-y",t),
        c(matrix[1,1],matrix[1,-1][index])),"\n",file=file,append=TRUE)
  
  
  
}else if (method == 5) {
  
  ###### 5 METHOD ######  
  # extreme S
  for (Quantile in c(1,2,4,5)) {
    temp = quantile(S, c(Quantile/100, 1 - Quantile/100))
    temp
    ind = (S <= temp[1])|(S >= temp[2])
    sum(ind)
    Quantile/100*2*N*M
    S_extremte = S[ind]
    Y_extreme = c(as.numeric(S[ind] >= temp[2]),Y_label)
    X_extreme = rbind(X[ind, ], X_label)
    str(X_extreme)
    str(Y_extreme)
    
    cvfit5 = cv.glmnet(y=Y_extreme,x=X_extreme,
                       family = "binomial",alpha = 1)
    # plot(cvfit)
    # cvfit$lambda.1se
    print(coef(cvfit5, s="lambda.min")[1:10])
    
    gamma = (coef(cvfit5, s="lambda.min")[2:(d+1)]/coef(cvfit5, s="lambda.min")[2])
    
    
    fit=glm(Y_label~X_label%*%gamma,
            #family="binomial",
            family="binomial")
    print(fit[[1]])
    res=c(fit[[1]][1],fit[[1]][2]*gamma)
    print(res)
    
    cat(c(seed,iter,n,N,d,M,paste0("extreme-",Quantile),
          c(fit[[1]][1],fit[[1]][2]*gamma)
    ),"\n",file=file,append=TRUE)
    
    ##### debias #####
    deb_res = debias2(beta_fit =  t(t(res)), 
                      beta_star,
                      X_label = X_label,
                      n_lst = n,  tau_lst = 0.001 * 1:1000 * sqrt((log(d)) / n))  
    
    as.numeric(deb_res[[3]])
    as.numeric(deb_res[[4]])
    sum(deb_res[[4]])
    
    cat(c(seed,iter,n,N,d,M,paste0("extreme-",Quantile,"-test"), 
          as.numeric(deb_res[[3]])),
        "\n",file=file,append=TRUE)
    
    cat(c(seed,iter,n,N,d,M,paste0("extreme-",Quantile,"-cover"), 
          as.numeric(deb_res[[4]])),
        "\n",file=file,append=TRUE)
    
    cat(c(seed,iter,n,N,d,M,paste0("extreme-",Quantile,"-deb"), 
          as.numeric(deb_res[[1]])),
        "\n",file=file,append=TRUE)
    
    cat(c(seed,iter,n,N,d,M,paste0("extreme-",Quantile,"-var"), 
          as.numeric(deb_res[[2]])),
        "\n",file=file,append=TRUE)
    
    print("done debias!")
    
  }  
}else if (method == 6) {
  
  ###### 6 METHOD ###### 
  # PASS
  # least square
  print(str(X))
  print(str(S))
  cvfit6 = cv.glmnet(x=X, y=S, 
                     family = "gaussian",alpha = 1)
  # plot(cvfit)
  # cvfit$lambda.1se
  print(coef(cvfit6, s="lambda.min"))
  print(coef(cvfit6, s="lambda.min")[2:20]/coef(cvfit6, s="lambda.min")[2])
  
  cat(c(seed,iter,n,N,d,M,paste0("PASS"),
        as.vector(coef(cvfit6, s="lambda.min"))
  ),"\n",file=file,append=TRUE)
  
  
  
  
}else if(method == 8){
  ##### method 8 #######
  print("RMRCE")
  print(Sys.time())
  res = NULL
  for(m in 1:M) {
    cat("m =",m,"\n")
    gamma <- RMRCE(X_lst[[m]][,1:10],S_lst[[m]],0.01,5)
    
    round(gamma,3)
    gamma = c(gamma,rep(0,d-length(gamma)))
    
    fit=glm(Y_label~X_label%*%gamma,
            #family="binomial",
            family="binomial")
    print(fit[[1]])
    res=cbind(res,c(fit[[1]][1],fit[[1]][2]*gamma))
    
    print(Sys.time())
    
  }
  
  print(round(res[1:5,],3))
  
  index = rowSums(abs(res)>1e-6) >= M/2
  final_res = rowMeans(res)
  final_res[!index] = 0
  
  print("save result")
  cat(c(seed,iter,n,N,d,M,"RMRCE_0.01_5_10",
        final_res
  ),"\n",file=file,append=TRUE)
  
  
  
  print("RMRCE_CV")
  
  # cv <- RMRCE_cv(X[,1:10],S,c(0.001,0.01,0.1),c(1,5))
  # optlambda <- cv[which.min(cv[,3]),1]
  # optalpha <- cv[which.min(cv[,3]),2]
  # optlambda
  # optalpha
  # 
  # gamma <- RMRCE(XX[,1:10],S,optlambda,optalpha)
  # 
  # round(gamma,3)
  # 
  # fit=glm(Y_label~X_label%*%gamma,
  #         #family="binomial",
  #         family="binomial")
  # print(fit[[1]])
  # res=c(fit[[1]][1],fit[[1]][2]*gamma)
  # print(res)
  # 
  # cat(c(seed,iter,n,N,d,M,"RMRCE_CV",
  #       c(fit[[1]][1],fit[[1]][2]*gamma)
  # ),"\n",file=file,append=TRUE)
  
  
  
}  else if(method == 9){
  #####  9 #######
  print("RMRCE")
  print(Sys.time())
  
  index = order(abs(cor(S,X)),decreasing = TRUE)[1:10]
  print(index)
  
  res = NULL
  for(m in 1:M) {
    cat("m =",m,"\n")
    gamma_fit <- RMRCE(X_lst[[m]][,index],S_lst[[m]],0.01,5)
    
    round(gamma_fit,3)
    
    gamma = rep(0,d)
    gamma[index] = gamma_fit
    round(gamma,3)
    
    fit=glm(Y_label~X_label%*%gamma,
            #family="binomial",
            family="binomial")
    print(fit[[1]])
    res=cbind(res,c(fit[[1]][1],fit[[1]][2]*gamma))
    
    print(Sys.time())
    
  }
  
  print(round(res[1:10,],3))
  
  # index = rowSums(abs(res)>1e-6) >= M/2
  final_res = rowMeans(res)
  # final_res[!index] = 0
  round(final_res,3)
  
  print("save result")
  cat(c(seed,iter,n,N,d,M,"RMRCE_PreScreen",
        final_res
  ),"\n",file=file,append=TRUE)
  
  
  
  print("done")
  
  
  
}  









