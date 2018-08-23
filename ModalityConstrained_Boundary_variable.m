close all;
%%BOUNDARY NOT CONSTRAINED TO BE ZERO HERE...USES GLOBALSEARCH TOOLBOX, SO
%%SLOWER ALGORITHM.
n=1000;%sample size
M=2;%the modality
%N = 100;
% tt0 = (0:N)/N;
% dt=1/N;
% ft =0.75*normpdf(tt0,0.3,0.2).*(tt0>=0 & tt0<=1)+ 0.25*normpdf(tt0,0.75,1/8).*(tt0>=0 & tt0<=1);


[X,tt0,ft] = GenerateData(n);% Must have the density mentioned above selected in the file GenerateData, or a valid density with support on [0,1]

 A=0;B=1; % ASSUMED [0,1] support
X=(X-A)/(B-A);
X = sort(X);
%%
%%Find an unconstrained density estimator to obtain the modal height
T=100;%grid length
t1=0:1/(T-1):1;%dense grid on [0,1]
del = mean(diff(t1));

%%
%%% Nonparametric part starts here
%%fourier basis%%
K = 5;
T1 = 2*pi*t1;
Phi(1:K,:) = sqrt(2)*sin((1:K)'*T1);
Phi(K+1:2*K,:) = sqrt(2)*cos((1:K)'*T1);

options = optimset('MaxFunEvals',60000,'MaxIter',90000,'Display','off');
KK=0:1:K;
Z=zeros(1,1);
start(1,:)=[zeros(1,2) repmat([0.5 1],1,(M-1)) 0.5 0.5];
r=2;
ll(1)=100;ll(2)=0;
Phiopt=Phi;

B1=-10^(-7)*ones(2*M -2,1);
opts = optimoptions(@fmincon,'Algorithm','interior-point','MaxFunEvals',60000,'MaxIter',90000);
while (iter <=length(KK) ) %& (ll(iter)<ll(iter-1))
    p1=Phi(1:KK(iter),:);
    p2=Phi(K+1:K+KK(iter),:);
    p=[p1;p2];
    m_temp=size(p,1);
  
   lb=[-Inf*ones(1,m_temp) zeros(1,2*(M))];
if (M>1)
    ub=[Inf*ones(1,m_temp) 1 Inf*ones(1,2*M-3) 0.9 0.9];
else
    ub=[Inf*ones(1,m_temp) 0.9 0.9];
end
if (M>1)
    A1=zeros(2*M -2,m_temp+ 2*(M));
      for i=1:(2*M -3)
         A1(i,:)=[zeros(1,m_temp+i-1) (-1)^(i-1) (-1)^i zeros(1,2*M - i -1)];
      end
    A1(2*M -2,:)=[zeros(1,m_temp + 2*M -3) -1 1 0];
else
    A1=[];
    B1=[];
end
%% USES GLOBALSEARCH TOOLBOX
  problem=createOptimProblem('fmincon','objective',@(c)FormMmodalLikelihoodboundaryFromC(c,X,T,p,t1,M),'x0',start(1,:),'Aineq',A1,'bineq',B1,'lb',lb,'ub',ub,'options',opts);
  gs = GlobalSearch('Display','off');
        [temp,tempval]=run(gs,problem);
    dcell{iter}=temp;
    ll(iter)=tempval;
    dd=dcell{iter};
    sd=length(dd(1:m_temp));
    d1=dd(1:sd/2);
    d2=dd(((sd/2)+1):m_temp);
    ddd=dcell{iter}(1:m_temp);
    clear start;clear temp;clear tempval;
    start(1,:)=[d1,Z,d2,Z repmat([0.5 1],1,(M-1)) 0.5 0.5];
    iter=iter+1;
end
[~,i0]=min(ll(2:end));
i0=i0+1;
p1=Phi(1:KK(i0),:);
p2=Phi(K+1:K+KK(i0),:);
Phiopt=[p1;p2];
m=size(Phiopt,1);
d=dcell{i0};
gamEst = FormGammaFromC(d(1:m),Phiopt);
l=[1 d((m+1):end)];
g_ind=floor([(1:(100/(2*M)):100) 100]);
gval=ones(1,length(g_ind));
gval(1)=l(end);
gval(end)=l(end-1);
gval(2:(end-1))=l(1:(end-2));

g=interp1(t1(g_ind),gval,t1,'pchip');
g=round(g.*(10^5))/(10^5);
fn = interp1(t1, g, (t1(end)-t1(1)).*gamEst + t1(1));
fn=fn/(sum(fn)/T);
fn=fn/(B-A);
fn=fn/(sum(fn)*(B-A)/T);
t=t1*(B-A)+A;



figure(1);
clf;
plot(t,fn,'b:','Linewidth',3);hold on;
%plot(datast(a,:),zeros(1,length(datast(a,:))),'r*','Linewidth',3);
set(gca,'fontsize',18);
