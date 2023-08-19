
# Loading necessary libraries
library(glmnet)
library(MASS)
library(RMRCE)

# Defining the logit function
logit <- function(a) {
  (exp(a)) / (1 + exp(a))
}

# Defining the g function
g <- function(a) {
  (exp(a)) / (1 + exp(a))
}

# Defining the gradient of g function
g_grad <- function(a) {
  (exp(a)) / (1 + exp(a))^2
}

# Defining the K_h function
K_h <- function(X, h) {
  dnorm(X / h) / h
}

# Defining the gradient of K function
K_grad <- function(X) {
  -X * dnorm(X)
}

# Defining the K_h gradient function
K_h_grad <- function(X, h) {
  K_grad(X / h) / h^2
}

# Defining the f_hat function
f_hat <- function(X, x, beta, S, h) {
  xx <- matrix(rep(x, nrow(X)), nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  sum(K_h((X - xx) %*% beta, h) * S) / sum(K_h((X - xx) %*% beta, h))
}

# Defining the gradient of f_hat function
f_hat_grad <- function(X, x, beta, S, h) {
  xx <- matrix(rep(x, nrow(X)), nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  A <- t(K_h_grad((X - xx) %*% beta, h) * S) %*% (X - xx)
  B <- sum(K_h((X - xx) %*% beta, h))
  C1 <- sum(K_h((X - xx) %*% beta, h) * S)
  C2 <- t(K_h_grad((X - xx) %*% beta, h)) %*% (X - xx)
  A / B - (C2 * C1) / B^2
}




loss <- function(b, lambda) {
  # b: length d-1
  b <- matrix(b)
  A <- apply(X_m, 1, function(o) {f_hat(X_m, o, beta, S_m, h)})
  B <- apply(X_m, 1, function(o) {f_hat_grad(X_m, o, beta, S_m, h)})
  
  if(m <= (M / 2)) {
    w <- 1
  } else {
    w <- 1 / (A * (1 - A))
  }
  
  w[abs(w) > 100] <- 100
  
  mean(w * (S_m - A - c(t(rbind(1, b) - beta) %*% B))^2) + lambda * sum(abs(b))
}

hess <- function(X_m, S_m, beta, h, m, N) {
  temp <- apply(X_m, 1, function(o) {f_hat_grad(X_m, o, beta, S_m, h)})
  A <- apply(X_m, 1, function(o) {f_hat(X_m, o, beta, S_m, h)})
  
  if(m <= (M / 2)) {
    w <- rep(1, length(A))
  } else {
    w <- 1 / (A * (1 - A))
  }
  
  w[abs(w) > 100] <- 100
  
  (temp %*% diag(w) %*% t(temp)) / N
}

grad <- function(X_m, S_m, beta, h, m, N) {
  A <- apply(X_m, 1, function(o) {f_hat(X_m, o, beta, S_m, h)})
  B <- apply(X_m, 1, function(o) {f_hat_grad(X_m, o, beta, S_m, h)})
  
  if(m <= (M / 2)) {
    w <- rep(1, length(A))
  } else {
    w <- 1 / (A * (1 - A))
  }
  
  w[abs(w) > 100] <- 100
  
  B %*% diag(w) %*% matrix(S_m - A + c(t(beta) %*% B)) / N
}

AR_cov <- function(p, acorr) {
  cov_mat <- diag(rep(1, p))
  cor_vec <- acorr^c(1:p)
  for (i in 1:(p - 1)) {
    cov_mat[i, (i + 1):p] <- cor_vec[1: (p - i)]
    cov_mat[(i + 1):p, i] <- cor_vec[1: (p - i)]
  }
  
  return(cov_mat)
}

rmvpois <- function(n, mu, Sigma) {
  Cor <- Sigma
  SDs <- sqrt(diag(Sigma))
  d <- length(SDs)
  normd <- mvrnorm(n, rep(0, d), Cor)
  unif <- pnorm(normd)
  data <- t(qpois(t(unif), mu))
  
  return(data)
}

loss_step2 <- function(b, YYY, XXX, lamb, Nm, weight, sigma.square) {
  # b: length d-1
  b <- matrix(b)
  avg <- mean(weight * (YYY - XXX %*% b)^2)
  
  return(Nm * (avg) / sigma.square + 
           (length(b) - sum(abs(b) < 1e-8)) * log(Nm))
}

loss_step3 <- function(gamma, Grad_final, Hess_final, lambda, n, N_lst) {
  # gamma: length d-1
  gamma <- matrix(gamma)
  sum_val <- t(rbind(1, gamma)) %*% Hess_final %*% rbind(1, gamma) -
    2 * t(rbind(1, gamma)) %*% Grad_final
  
  return(sum_val / (sum(N_lst) + n) + lambda / (sum(N_lst) + n) * sum(abs(gamma)))
}







debias2 <- function(beta_fit, beta_star, X_label, n_lst, tau_lst = NULL){
  
  ss=0
  X_obs = cbind(1,X_label)
  for(i in 1:(dim(X_obs)[1])){
    ss=ss+t(t(X_obs[i,]*c(Y_label[i]-g(X_obs[i,]%*%beta_fit))))%*%
      t(X_obs[i,]*c(Y_label[i]-g(X_obs[i,]%*%beta_fit)))
  }
  ss[1:5,1:5]
  
  HH = t(cbind(1,X_label))%*%diag(Vectorize(g_grad)(cbind(1,X_label)%*%beta_fit))%*%cbind(1,X_label)  
  
  GG = t(cbind(1,X_label))%*%matrix((Y_label-Vectorize(g)(cbind(1,X_label)%*%beta_fit))) +
    HH %*% beta_fit
  
  H_lst = list(HH)
  g_lst = list(GG)
  
  M <- 1
  n <- n_lst
  p <- length(H_lst[[1]][1, ]) - 1
  
  
  
  
  X_lst <- vector('list', M)
  Y_lst <- vector('list', M)
  
  for (i in 1:M){
    H <- H_lst[[i]]
    d <- g_lst[[i]]
    mat_all <- cbind(H, d)
    mat_all <- rbind(mat_all, t(c(d, max(mat_all) + n)))
    
    svd_result <- svd(mat_all)
    s_value <- svd_result$d
    s_mat <- diag(sqrt(s_value))[ ,1:(min(p + 1, n_lst[[i]]) + 1)]
    data_all <- svd_result$u %*% s_mat
    
    X <- t(data_all[-length(data_all[ ,1]),])
    Y <- data_all[length(data_all[ ,1]),]
    X_lst[[i]] <- X
    Y_lst[[i]] <- Y
  }
  
  (t(X)%*%X)[1:5,1:5]
  HH[1:5,1:5]
  Y[1:5]
  
  # (t(X)%*%X/n)[1:5,1:5]
  
  beta_deb_mat = var_mat = matrix(0, p + 1, M)
  for (j in 1:(p + 1)) {
    ej <- rep(0, p + 1)
    ej[j] <- 1
    m=1
    X_mj <- X_lst[[m]][, -j]
    X_j <- X_lst[[m]][,j]
    
    
    set.seed(1)
    if (is.null(tau_lst)){
      
      model_fit <- cv.glmnet(X_mj, X_j,  
                             # lambda = tau_lst,
                             intercept = F,
                             standardize = F)
      tau <- model_fit$lambda.min
      model_fit <- glmnet(X_mj, X_j, lambda = tau, intercept = F,
                          standardize = F)
      
    } else {
      model_fit <- cv.glmnet(X_mj , X_j ,  
                             lambda = tau_lst,
                             intercept = F,
                             standardize = F)
      tau <- model_fit$lambda.min
      model_fit <- glmnet(X_mj , X_j , lambda = tau, intercept = F,
                          standardize = F)
    }
    
    round(as.vector(model_fit$beta),5)
    
    u_j <- rep(1, p + 1)
    u_j[-j] <- - model_fit$beta
    sigma2 <- mean((X_j - X_mj %*% model_fit$beta)^2) + tau * sqrt(n / n_lst[m]) * sum(abs(model_fit$beta))
    sigma2
    u_m <- as.vector(u_j / sigma2)
    
    beta_deb <- beta_fit[j, m] + 
      t(u_m) %*% (g_lst[[m]] - H_lst[[m]] %*% beta_fit[,m]) / n_lst[m]
    beta_deb
    beta_deb_mat[j, m] <- beta_deb
    
    var_mat[j, m] = t(u_m)%*%(ss/(dim(X_obs)[1])) %*% u_m
    
 }
  
  testresult=abs(beta_deb_mat) > (qnorm(0.975)*sqrt(var_mat)/sqrt(n))
  cover = (c(beta_deb_mat - qnorm(0.975)*sqrt(var_mat)/sqrt(n)) <= c(alpha_star,beta_star)) & 
    (c(beta_deb_mat + qnorm(0.975)*sqrt(var_mat)/sqrt(n)) >= c(alpha_star,beta_star))
  cover = as.numeric(cover)
  
  LB = c(beta_deb_mat - qnorm(0.975)*sqrt(var_mat)/sqrt(n))
  UB = c(beta_deb_mat + qnorm(0.975)*sqrt(var_mat)/sqrt(n))
  
  
  return(list(beta_deb_mat, var_mat, testresult, cover, LB, UB))
  

}





SASH <- function(X_label, Y_label, X_lst, S_lst) {
  
  #### SL ####
  cvfit <- cv.glmnet(
    X_label,
    Y_label,
    family = "binomial",
    alpha = 1
  )
  
  
  
  if (coef(cvfit, s = "lambda.min")[2] != 0) {
    index0 <- 1
  } else {
    index0 <- order(abs(coef(cvfit, s = "lambda.min"))[2:3], decreasing = T)[1]
  }
  
  index <- 1:d
  index[index0] <- 1
  index[1] <- index0
  
  beta_initial <- coef(cvfit, s = "lambda.min")[-1][index]
  
  gamma_initial <- matrix(
    beta_initial / ifelse(beta_initial[1] == 0, 1, beta_initial[1])
  )
  gamma_initial[1] <- 1
  
  beta_initial <- c(coef(cvfit, s = "lambda.min")[1], beta_initial)
  
  X_label <- X_label[, index]
  for (m in 1:M) {
    X_lst[[m]] <- X_lst[[m]][, index]
    X_lst[[m]] <- as.matrix(X_lst[[m]])
  }
  
  X <- X[, index]
  
  SL <- c(beta_initial[1], beta_initial[-1][index])
  
  #### debias ####
  deb_res <- debias2(
    beta_fit = t(t(beta_initial)),
    beta_star,
    X_label = X_label,
    n_lst = n,
    tau_lst = 0.001 * 1:1000 * sqrt((log(d)) / n)
  )
  
  SL_test <- c(deb_res[[3]][1], deb_res[[3]][-1][index])
  SL_cover <- as.numeric(c(deb_res[[4]][1], deb_res[[4]][-1][index]))
  SL_deb <- as.numeric(c(deb_res[[1]][1], deb_res[[1]][-1][index]))
  SL_var <- as.numeric(c(deb_res[[2]][1], deb_res[[2]][-1][index]))
  SL_LB <- as.numeric(c(deb_res[[5]][1], deb_res[[5]][-1][index]))
  SL_UB <- as.numeric(c(deb_res[[6]][1], deb_res[[6]][-1][index]))
  
  beta_label <- c(coef(cvfit, s = "lambda.min")[1], coef(cvfit, s = "lambda.min")[-1][index])
  
  final <- NULL
  Grad <- Hess <- NULL
  beta_debias <- variance <- NULL
  lambda_final <- sigma.square_final <- NULL
  
  
  
    #### step 1 ####
    for(m in 1:M){
      cat("m =",m,"\n")

      cat("iter",iter,"\n")
      
      
      beta=c(gamma_initial)
      
      beta_pool=NULL
      
      for (t in 1:iter){
        cat("t =",t,"\n")
     
        
        BIC_step2=NULL
        beta_step2=NULL
    
        
        
        curr_h = (log(N_lst[m])/N_lst[m])^(1/5)
        h_lst=c(curr_h)
        
        
        curr_min_BIC = 1000000
        curr_beta = NULL
        curr_h = 0
        
        
        for(j in 1:length(h_lst)){
          
          h=h_lst[j]
         
          A=apply(X_lst[[m]],1,function(o){f_hat(X_lst[[m]],o,beta,S_lst[[m]],h)})
          
      
          
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
        
          if (t==iter) {
            sigma.square_final=c(sigma.square_final,sigma.square)
          }
          BIC_step2=NULL
          beta_step2=NULL
          lambda_lst=0.005*(1:600)*sqrt(log(d)/N)
          B1=B2=NULL
          for (lambda in lambda_lst){
            
            cv.result=glmnet(t(B[-1,]),c(S_lst[[m]]-A+t(beta)%*%B-B[1,]), 
                             weights=w,family="gaussian",
                             lambda=lambda,
                             alpha=1,intercept = FALSE)
            
            
            beta_new=coef(cv.result, s="lambda.min")[-1]
            
            avg=mean(w*(c(S_lst[[m]]-A+t(beta)%*%B-B[1,])-t(B[-1,])%*%beta_new)^2)

            BIC_step2=c(BIC_step2,
                        loss_step2(beta_new,
                                   c(S_lst[[m]]-A+t(beta)%*%B-B[1,]),t(B[-1,]),
                                   lamb=cv.result$lambda.min/sqrt(log(d)/N)
                                   ,N_lst[[m]],w,sigma.square))
            beta_step2=cbind(beta_step2,beta_new)
            
          } 
          
          beta=c(1,beta_step2[,which.min(BIC_step2)])

          
          if(min(BIC_step2) < curr_min_BIC){
            curr_min_BIC = min(BIC_step2)
            curr_beta = beta
            curr_h = h
          }
          
        } 
      } # end iter
     
      lambda_final=c(lambda_final,lambda_lst[which.min(BIC_step2)])
      
      beta = curr_beta
      h = curr_h
      
   
      
      Grad[[m]]=grad(X_lst[[m]],S_lst[[m]],beta,h,m,N)
      Hess[[m]]=hess(X_lst[[m]],S_lst[[m]],beta,h,m,N)
      
      final=cbind(final,beta)
      
      cat("\n")
      
    } # end m
    
    rownames(final)=1:d

    

    
    G=Grad
    H=Hess
    
    Grad_pool=Hess_pool=0
    for(m in 1:M){
      
      Grad_pool=Grad_pool+Grad[[m]]*N
      Hess_pool=Hess_pool+Hess[[m]]*N
    }
    
    Hess_pool[1:5,1:5]
    eigen(Hess_pool)$values
    beta_dagger=NULL
    
    
    
    
    #### step 3 ####
    N_lst=rep(N,M)

    
    beta_tilde=apply(final,1,mean)
    fit=glm(Y_label~X_label%*%beta_tilde,
            #family="binomial",
            family="binomial")
  
    
    for(t in 1:iter){
      a_tilde=fit[[1]][1]
      b_tilde=fit[[1]][2]
      
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
        
        bb=model$beta
        bb[1:15]
        beta_step3=cbind(beta_step3,bb)

        
        
        BIC=c(BIC,
              (length(bb)-sum(abs(bb)<1e-8))*log(n+sum(N_lst))+
                (n+sum(N_lst))*(loss_step3(model$beta,
                                           Grad_final,Hess_final,lambda,
                                           n, N_lst)))
        
      }
   
      
      beta_tilde=c(1,c(beta_step3[,which.min(BIC)]))
      
      fit=glm(Y_label~X_label%*%beta_tilde,
              family="binomial")
      
      
    }
    
    
    
    a_final=fit[[1]][1]
    b_final=fit[[1]][2]
    beta_dir_final=beta_tilde
    matrix=NULL
    matrix=rbind(matrix,c(a_final,b_final*beta_dir_final))
    
    beta_dagger=NULL
    beta_dagger=rbind(beta_dagger,beta_dir_final)
    
    SASH =  c(matrix[1,1],matrix[1,-1][index])
    
    ##### debias #####
    deb_res = debias2(beta_fit =  t(matrix), 
                      beta_star,
                      X_label = X_label,
                      n_lst = n,  tau_lst = 0.001 * 1:1000 * sqrt((log(d)) / n))  
    
    
    SASH_test=as.numeric(c(deb_res[[3]][1], deb_res[[3]][-1][index]))
    SASH_cover=as.numeric(c(deb_res[[4]][1], deb_res[[4]][-1][index]))
    SASH_deb=as.numeric(c(deb_res[[1]][1], deb_res[[1]][-1][index]))
    SASH_var=as.numeric(c(deb_res[[2]][1], deb_res[[2]][-1][index]))
    SASH_LB=as.numeric(c(deb_res[[5]][1], deb_res[[5]][-1][index]))
    SASH_UB=as.numeric(c(deb_res[[6]][1], deb_res[[6]][-1][index]))
    
    
    gamma_initial = beta_tilde
  
    
    return(list(
      SL, SL_test, SL_cover, SL_deb, SL_var, SL_LB, SL_UB,
      SASH, SASH_test, SASH_cover, SASH_deb, SASH_var, SASH_LB, SASH_UB
    ))
  
}





IPD <- function(X_label, Y_label, X_lst, S_lst) {
  
  cvfit=cv.glmnet(X_label,Y_label,
                  # lambda=seq(0.05,1.5,0.01)*sqrt(log(d)/n),
                  family="binomial",alpha=1)

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
      
      
      sigma.square=mean(ww*(S_m-A)^2)
      sigma.square_final=c(sigma.square_final,
                           sigma.square)
      
  
    }
    
 
    
    XX=diag(sqrt(w))%*%t(BB[-1,])
  
    YY = diag(sqrt(w))%*%c(S-AA+t(beta)%*%BB-BB[1,])
 
    
    loss_IPD=function(b,YYY,XXX,lamb,N_lst,weight,sigma.square_final){
      # b: length d-1
      b=matrix(b)
      pro=weight*(YYY-XXX%*%b)^2
      result_IPD=NULL
      for(m in 1:length(N_lst)){
        avg=mean(pro[((m-1)*N+1):(m*N)])
       
        result_IPD=c(result_IPD,
                     N_lst[m]*(avg)/sigma.square_final[m])
      }
      
      
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

      BIC_step2=c(BIC_step2,
                  loss_IPD(beta_new,
                           c(S-AA+t(beta)%*%BB-BB[1,]),
                           t(BB[-1,]),
                           lamb=lambda
                           ,N_lst,w,sigma.square_final))
      beta_step2=cbind(beta_step2,beta_new)
      
    } 
    

    IPD_lambda=lambda_lst[which.min(BIC_step2)]
    # print(lambda_lst[which.min(BIC_step2)]/sqrt(log(d)/N))
    beta=c(1,beta_step2[,which.min(BIC_step2)])
    
    
    beta_pool=cbind(beta_pool,beta)
    
  
  
    # beta
    fit_last=glm(Y_label~X_label%*%beta,family="binomial")


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
   
    Grad_pool=Grad_pool+Grad[[m]]*N
    Hess_pool=Hess_pool+Hess[[m]]*N
  }
  
  Hess_pool[1:5,1:5]
  eigen(Hess_pool)$values
  

  N_lst = rep(N,M)

  
  fit=glm(Y_label~X_label%*%beta,
          #family="binomial",
          family="binomial")
  print(fit[[1]])
  
  beta_tilde=beta
  for(t in 1:(iter)){
    a_tilde=fit[[1]][1]
    b_tilde=fit[[1]][2]
    
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
      
      
      
      BIC=c(BIC,
            (length(bb)-sum(abs(bb)<1e-8))*log(n+sum(N_lst))+
              (n+sum(N_lst))*(loss_step3(model$beta,
                                         Grad_final,Hess_final,lambda, n, N_lst)))
      
    }
 
    
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
  
  IPD=c(matrix[1,1],matrix[1,-1][index])
  
  

  deb_res = debias2(beta_fit =  t(matrix), 
                    beta_star,
                    X_label = X_label,
                    n_lst = n,  tau_lst = 0.001 * 1:1000 * sqrt((log(d)) / n))  
  
  as.numeric(deb_res[[3]])
  
  IPD_test=as.numeric(c(deb_res[[3]][1], deb_res[[3]][-1][index]))
  
  IPD_cover=as.numeric(c(deb_res[[4]][1], deb_res[[4]][-1][index]))
  
  IPD_deb=as.numeric(c(deb_res[[1]][1], deb_res[[1]][-1][index]))
  
  IPD_var=as.numeric(c(deb_res[[2]][1], deb_res[[2]][-1][index]))

  return(list(
    IPD, IPD_test, IPD_cover, IPD_deb, IPD_var
  ))
  
}



