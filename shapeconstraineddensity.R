require("pracma")
require("ggplot2")
rm(list=ls())
t=seq(0,1,1/99); # grid

### generating simulated data ########

ft=0.75*dnorm(t,0.3,0.05) +0.25*dnorm(t,0.75,1/8);
#normcons= 0.75*(pnorm(1,0.3,0.2) - pnorm(0,0.3,0.2))+0.25*(pnorm(1,0.75,1/8) - 
#                                                               pnorm(0,0.75,1/8))

#ft=(1/3)*dbeta(t,1,3) + (1/3)*dbeta(t,1,4) + (1/3)*dbeta(t,3,15);
#ft=(4/5)*dnorm(t,0.5,4) + (1/5)*dnorm(t,0.5,0.1);
#ft=dnorm(t,0.5,0.1);
#truedens=ft/normcons
ft=ft/sum(ft*mean(diff(t)));
truedens=ft;
Ft=cumsum(ft)*mean(diff(t));
n=100;
X=matrix(0,nrow=1,ncol=n); 


u=runif(n,0,1);
for (i in 1:n){
  X[i]=t[min(which(Ft>u[i]))]# SIMULATED DATA
}
#######################################

M=2; # NUMBER OF MODES
zz=6; # NO. OF BASIS ELEMENTS chosen ad-hoc, but can be optimized by running a loop.
T=2*pi*t;

### BASIS SET CREATION #######
phi=matrix(0,nrow=zz,ncol=length(t));
for (i in 1:(zz/2)){
  phi[i,]=sqrt(2)*sin(i*T);
  phi[(zz/2)+i,]=sqrt(2)*cos(i*T);
}
##################################
A=matrix(0,nrow=2*M -3, ncol=zz+2*M-2); # Ax<B constraint setup
if (M>1){
  for (i in 1:((2*M)-3)){
    A[i,]=c(rep(0,zz +i-1), (-1)^(i-1), (-1)^i, rep(0,2*M-i-3));
  }
}
B=matrix(0,nrow=2*M-3,ncol=1); 

####################################

### FORM GAMMA FROM COEFFICIENTS: FUNCTION ######
gam=function(c,phi){
  N=dim(phi)[2];
  v=c%*%phi;
  nc=norm(v,'f');
  
  if(nc==0){
    q=matrix(1,nrow=1,ncol=N);
  }else{
    q=cos(nc)*matrix(1,nrow=1,ncol=N) + sin(nc)*v/nc;#Exponential map
  }
  
  gam=cumsum(q*q);
  gam=(gam - gam[1])/(gam[N] - gam[1]);
  return(gam)
}

#### LIKELIHOOD COMPUTATION: FUNCTION #######
#L=function(d,X,phi,t,M){
  L=function(d){
    if(dim(phi)[1]>=length(d)){
      l=vector();#Usually happens when unimodal density, i.e. M=1
    } else{
  l=c(d[(dim(phi)[1]+1):(length(d))]);#height ratio vector
    }
    
  g_ind=floor(c(seq(1,100,100/(2*M)),100));#mode location in template
  
  
  gval=matrix(0,nrow=1,ncol=length(g_ind));
  
  gval[1]=0;
  gval[2]=1;
  gval[length(gval)]=0;
  if(length(gval)>3){
  gval[3:(length(gval)-1)]=l;
  }
  c1=d[1:dim(phi)[1]];
  gam0=gam(c1,phi);
  
  fp=interp1(t[g_ind],as.vector(gval),t,method="cubic");#cubic interpolation. One can also use spline. Results are almost identical.
  
  fn=interp1(t,fp,gam0,method="linear");#Can also use spline. Results are almost indistinguishable.
  fn=fn/(sum(fn)*mean(diff(t)));
  yy=interp1(t,fn,X,method="linear");
  L1=-sum(log(yy));# Use a penalized version , AIC if optimizing over the number of basis elements.
  L=L1;
  return(L);
}

################ DENSITY FROM PARAMETER VECTOR : FUNCTION ######
densest=function(d,phi,t,M){
  if(dim(phi)[1]>=length(d)){
    l=vector();
  } else{
    l=c(d[(dim(phi)[1]+1):(length(d))]);#estimated height ratio vector
  }
  g_ind=floor(c(seq(1,100,100/(2*M)),100));
  
  
  gval=matrix(0,nrow=1,ncol=length(g_ind));
  
  gval[1]=0;gval[2]=1;
  gval[length(gval)]=0;
  if(length(gval)>3){
    gval[3:(length(gval)-1)]=l;
  }
  c1=d[1:dim(phi)[1]];#estimated basis coefficients
  gam0=gam(c1,phi);#estimated warping function
  
  fp=interp1(t[g_ind],as.vector(gval),t,method="cubic");#template with estimated height ratio vector
  
  fn=interp1(t,fp,gam0,method="linear");#estimated density
  densest=fn/(sum(fn)*mean(diff(t)));
  return(densest);
}

start=c(rep(1,zz),0.5,1);
lower=c(rep(-10,zz),rep(.001,2*M -2));
if (M>1){
 upper=c(rep(10,zz), 0.99,10*sqrt(n)*rep(1,2*M-3));
}else{
upper=rep(10,zz);
}
coefest=fmincon(start,L,gr = NULL, method = "SQP",
                A = A, b = B, Aeq = NULL, beq = NULL,
                lb = lower, ub = upper, hin = NULL, heq = NULL,
                tol = 1e-06, maxfeval = 10000, maxiter = 5000)
d1=coefest$par
fest=densest(d1,phi,t,M)
plot(t,fest,"l")
lines(t,ft)