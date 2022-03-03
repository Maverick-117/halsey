# Problem 1

library(bench)

lgtc <- function(x){
  y <- 1/(1+exp(-x));
  return(y)
}

logres_GD <- function(x,y,theta,alpha,T_end) {
  LL_history = c();
  N = nrow(x);
  visited_thetas = matrix(rep(0,T_end*length(theta)),nrow=length(theta))
  for (t in 1:T_end) {
    p = 1/(1+exp(-(x %*% theta)))

    loss = -(((t(y) %*% log(p)) + (t(1-y) %*% log(1-p)) ));
    LL_history = c(LL_history, loss);
    new_theta = theta - ( alpha/N * ( t(x) %*% (p - y) ))
    theta = new_theta
    visited_thetas[,t] = theta;
  }
  conf_int = rep(0,length(theta));
  for (i in 1:length(theta)){
    conf_int[i] = 1.96 * sd(visited_thetas[i,(T_end-20):T_end])/sqrt(21);
  }
  return(list(theta, conf_int, LL_history))
}

logres_FS <- function(x,y,theta,alpha,T_end) {
  LL_history = c();
  conf_int = matrix(rep(0,T_end*length(theta)),nrow=length(theta));
  for (t in 1:T_end) {
    p = 1/(1+exp(-(x %*% theta)))

    loss = -(((t(y) %*% log(p)) + (t(1-y) %*% log(1-p)) ));
    LL_history = c(LL_history, loss);
    # the gradient descent part
    # new_theta = theta - ( alpha/N * ( t(x) %*% (p - y) ))

    # fisher scoring part
    var_mat = diag(as.vector(p)*as.vector((1-p)));
    f_i = t(x) %*% var_mat %*% x; f_i_inv = chol2inv(chol(f_i));
    new_theta = theta - alpha * f_i_inv %*% ( t(x) %*% (p - y) );
    theta = new_theta
    conf_int[,t] = 1.96 * diag(sqrt(f_i_inv))
  }

  return(list(theta, conf_int, LL_history))
}
# compare versus
# dat = read.table('SAheart.data.txt',header = TRUE, sep=",")
# glm(formula = chd ~ sbp+ tobacco + ldl+ adiposity+ famhist+ typea+ obesity+ alcohol+ age,data=dat,family=binomial)

dat = read.table('SAheart.data.txt',header = TRUE, sep=",")
pre_x = as.matrix(dat[,-c(1,ncol(dat))])

N = nrow(pre_x);
d = ncol(pre_x);
for (i in 1:N) {
  if (pre_x[i,"famhist"] == "Present") {
    pre_x[i,"famhist"] = 1
  } else {
    pre_x[i,"famhist"] = 0
  }
}

x = matrix(as.numeric(pre_x),ncol=d)
y = dat[,ncol(dat)]
T_end=100;
alpha = 0.1
theta = rep(0,each=ncol(pre_x))
est_GD <- logres_GD(x,y,rep(0,each=ncol(pre_x)),0.0001,T_end*100)
est_FS <- logres_FS(x,y,rep(0,each=ncol(pre_x)),0.1,T_end)
plot(est_GD[[3]])
plot(est_FS[[3]])

