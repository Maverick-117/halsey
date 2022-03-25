# Problem 3
# update HMMforward and HMMbackward to avoid roundoff errors? (see Note 4.5 in hmm_handout.pdf)
HMMforward <- function(y,P,E,v0) {
  n = length(y);
  s = length(v0);
  result = rep(1,s*n,nrow=s);
  a1 = c(v0[1]*E[1,y[1]],v0[2]*E[2,y[1]]);
  a = matrix(rep(0,n*2),nrow=2); a[,1] = a1;
  for (t in 2:n) {
    for (i in 1:s) {
      a[i,t] = E[i,y[t]] * sum(a[,t-1]*P[i,])
    }
  }
  return(a)
}

HMMbackward <- function(y,P,E,v0) {
  n = length(y);
  s = length(v0);
  result = rep(1,s*n,nrow=s);
  b = matrix(rep(0,n*2),nrow=2); b[,n] = rep(1,s);
  for (t in seq(from=n-1,to=1,by=-1)) {
    for (i in 1:s) {
      b[i,t] = sum(P[,i]*E[,y[t+1]]*b[,t+1])
    }
  }
  return(b)
}

hiddenStateInference <- function(a,b,T_end){
  N_s = nrow(a);
  mP = rep(0,T_end);
  Pb = rep(0,T_end);
  marg = matrix(rep(0,2*T_end),nrow=2);
  for (t in 1:T_end) {
    denom = sum(a[,t]*b[,t]);
    for (s in 1:N_s) {
      marg[s,t] = (a[s,t]*b[s,t])/denom
    }
    # mP[t] = which.max(pre_mP)
    # Pb[t] = max(pre_mP)
  }
  return(marg)
}

baum <- function(em_row, y, T_end, K_end) {
  # em_row = Et[1,]; y = att; K_end = 100;
  em_row_rand = rep(1/6,6);
  tp = matrix(rep(0.5,4),nrow=2);
  em = matrix(c(em_row, em_row_rand));
  v0 = rep(1/2,2);
  v = matrix(rep(0,K_end*2),nrow=2);
  P = array(rep(0,K_end*4),c(2,2,K_end));
  E = array(rep(0,K_end*12),c(2,6,K_end));

  for (k in 1:K_end) {
    E[1,,k] = em_row;
  }

  v[,1] = v0; P[,,1] = tp; E[2,,1] = em_row_rand;

  g = array(rep(0,T_end*4),c(2,2,T_end));

  for (k in 1:(K_end-1)){
    a = HMMforward(y,P[,,k], E[,,k], v[,k]);
    b = HMMbackward(y,P[,,k], E[,,k], v[,k]);
    gamma = hiddenStateInference(a,b,T_end);
    # calculating g
    for (t in 2:T_end) {
      for (i in 1:2){
        for (j in 1:2) {
          g[i,j,t] = a[i,t-1] * P[i,j,k] * E[j,y[t],k] * b[j,t];
        }
      }
      g[,,t] = g[,,t]/sum(g[,,t])
    }
    # calculating v_k and P_k
    for (i in 1:2) {
      v[i,k+1] = gamma[i,1];
      for (j in 1:2) {
        P[i,j,k+1] = sum(g[i,j,]);
      }
      P[i,,k+1] = P[i,,k+1]/sum(P[i,,k+1])
    }
    # calculating E_k for the unfair state
    for (l in 1:6) {
      E[2,l,k+1] = sum(gamma[2,]*(y == l));
    }
    E[2,,k+1] = E[2,,k+1]/sum(gamma[2,]);
  }
  return(list(v,P,E))
}

# Problem 3a

P = matrix(c(0.98,0.02,0.05,0.95),nrow=2);
E = t(matrix(c(rep(1,6)/6,c(1,1,5,1,1,1)/10),ncol=2));
v0 = rep(1,2)/2;
N_s = 2;
# simulating HMM
T_end = 100;
#v = matrix(rep(1,N_s * T_end),nrow = N_s);

q = rep(-1,T_end)
a = rep(-1,T_end)
v = v0;
for (t in 1:T_end){
  q[t] = sample(x = c(1,2), size = 1,
                prob = v, replace=T)
  a[t] = sample(x = c(1,2,3,4,5,6), size = 1,
    prob = E[q[t],], replace=T)
  v = P %*% v;
}

par(mfcol=c(2,1), mar = numeric(4), oma = c(4, 4, .5, .5),
    mgp = c(2, .6, 0))
plot(1:T_end,q,xlab='',ylab='',xaxt='n',yaxt='n',ylim=c(0,3),type="b")
legend("top",legend=c("Fairness State"),col=1,pch=1)
axis(2,at = c(1,2),las=1,labels=c("Fair","Unfair"))
plot(1:T_end,a,col=2,ylab='States',las=1,ylim=c(0,8),yaxt='n',type="b")
axis(2,at=c(1,2,3,4,5,6),las=1)
mtext("Iterations", col = "black", adj=0.5, padj=21)
legend("top",legend=c("Dice Roll State"),col=2,pch=1)

bR = HMMbackward(a,P,E,v0)
aR = HMMforward(a,P,E,v0)

aRn = sum(aR[,T_end])
bR0 = sum(v0*E[,a[1]]*bR[,1])

infrd = hiddenStateInference(aR,bR,T_end)

# Problem 3b

best_est = rep(0,T_end); for(t in 1:T_end) {best_est[t] = which.max(infrd[,t]);}
pltState = matrix(c(q,best_est),ncol=2);
matplot(pltState,pch=1,type="b",yaxt = "n",xaxt='n',xlab="Iterations",ylab="States",las=1,ylim=c(0.5,3))
legend("top",legend=c("True State","Inferred State"),col=1:2,pch=1)
axis(2,at = c(1,2),las=1,labels=c("Fair","Unfair"),ylab="States")
axis(1,labels="")
plot(infrd[1,],xlab='Iterations',yaxt='n',ylim=c(-0.01,1.3),ylab="Prob",col=3)
axis(2,las=1,at=c(0,0.2,0.4,0.6,0.8,1.0),las=1)
legend("top",legend=c("Marginal Probability of Fair State"),col=3,pch=1)
mtext("Iterations", col = "black", adj=0.5, padj=21)

nmz <- function(x) {
  return(x/sum(x));
}
K_end = 100;

# Problem 3c

result <- baum(E[1,],a,100,100)
print("final v:")
print(result[[1]][,K_end])
print("final P:")
print(result[[2]][,,K_end])
print("final E:")
print(result[[3]][,,K_end])
