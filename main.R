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


