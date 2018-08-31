%% Assumption that Boundary Values are zero
n=100;%sample size
M=2;%the modality

[X,tt0,ft] = GenerateData(n);% Simulated examples in GenerateData file
   A=min(X) - std(X)/sqrt(n); % Estimated boundaries

  B=max(X) + std(X)/sqrt(n);  % Estimated boundaries

X=(X-A)/(B-A);
X = sort(X);
%%
%%Find an unconstrained density estimator to obtain the modal height
T=100;%grid length
t1=0:1/(T-1):1;%dense grid between the bounds
del = mean(diff(t1));

%%
%%% Nonparametric part starts here
%%fourier basis%%
K = 8; % No. of basis elements/2; 
T1 = 2*pi*t1;
Phi(1:K,:) = sqrt(2)*sin((1:K)'*T1);%sine components
Phi(K+1:2*K,:) = sqrt(2)*cos((1:K)'*T1);%cosine components

%%%%%%%%%%%%%%%%%%%%
options = optimset('MaxFunEvals',60000,'MaxIter',90000,'Display','off');
KK=0:1:K;
Z=zeros(1,1);
start(1,:)=[zeros(1,2) repmat([0.5 1],1,(M-1))];

iter=2;
ll(1)=100;ll(2)=0;
Phiopt=Phi;
B1=-0.5*ones(2*M -3,1);

while (iter <=length(KK) ) 
    p1=Phi(1:KK(iter),:);
    p2=Phi(K+1:K+KK(iter),:);
    p=[p1;p2];
    m_temp=size(p,1);
   
    lb=[-Inf*ones(1,m_temp) 10^(-3)*ones(1,2*(M-1))];
if (M>1)
 
    ub=[Inf*ones(1,m_temp) 0.99 10*sqrt(n)*ones(1,2*M-3)];
else
    ub=Inf*ones(1,m_temp);
end
A1=zeros(2*M -3,m_temp+ 2*(M-1));
for i=1:(2*M -3)
    A1(i,:)=[zeros(1,m_temp+i-1) (-1)^(i-1) (-1)^i zeros(1,2*M - i -3)];
end
    for startiter=1:1 
         temp(startiter,:)=fmincon(@(c)FormMmodalLikelihoodFromC(c,X,T,p,t1,M),start(startiter,:),A1,B1,[],[],lb,ub,[],options);
         tempval(startiter)=FormMmodalLikelihoodFromC(temp(startiter,:),X,T,p,t1,M);
    end
    [~,strtpt]=min(tempval);
    dcell{iter}=temp(strtpt,:);
    ll(iter)=tempval(strtpt);
    dd=dcell{iter};
    sd=length(dd(1:m_temp));
    d1=dd(1:sd/2);
    d2=dd(((sd/2)+1):m_temp);
    ddd=dcell{iter}(1:m_temp);
    clear start;clear temp;clear tempval;
    start(1,:)=[d1,Z,d2,Z repmat([0.5 1],1,(M-1))];
    iter=iter+1;
end
[~,i0]=min(ll(2:end));
i0=i0+1;
p1=Phi(1:KK(i0),:);
p2=Phi(K+1:K+KK(i0),:);
Phiopt=[p1;p2];
m=size(Phiopt,1);
d=dcell{i0};
gamEst = FormGammaFromC(d(1:m),Phiopt); % Construct \gamma from coefficient vector
l=[1 d((m+1):end)]; %height ratio vector
g_ind=floor([(1:(100/(2*M)):100) 100]); %location of critical point in template
gval=ones(1,length(g_ind));
gval(1)=eps; %setting to very small value (close to 0)
gval(end)=eps;
gval(2:(end-1))=l;

g=interp1(t1(g_ind),gval,t1,'pchip'); %can also use linear or spline interpolation
g=round(g.*(10^5))/(10^5);% for numerical stability
fn = interp1(t1, g, (t1(end)-t1(1)).*gamEst + t1(1));
fn=fn/(sum(fn)/T); %normalize to a density on [0,1]
fn=fn/(B-A);%scaled to the density on [A,B]
fn=fn/(sum(fn)*(B-A)/T);%re-normalize to miniize numerical error 
t=t1*(B-A)+A;

figure(1);
clf;
plot(t,fn,'b:','Linewidth',3);hold on;
plot(datast(a,:),zeros(1,length(datast(a,:))),'r*','Linewidth',3);
set(gca,'fontsize',18);

