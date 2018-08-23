function fn=Modalconstrainedcde(X,Y,y0,A,B,M,h,choice)
%choice is for the basis set to use
n=length(X);
X=(X-A)/(B-A);
X = sort(X);
%%
%%Find an unconstrained density estimator to obtain the modal height
T=100;%grid length
t1=0:1/(T-1):1;%dense grid between the bounds
% [test1,~,hh]=ksdensity(Y,y0);
% h=hh/sqrt(test1(1));
%%
%%% Nonparametric part starts here
%%
%%fourier basis%%
if (choice==1)
K = 2;
T1 = 2*pi*t1;
Phi(1:K,:) = sqrt(2)*sin((1:K)'*T1);
Phi(K+1:2*K,:) = sqrt(2)*cos((1:K)'*T1);

%options = optimset('MaxFunEvals',60000,'MaxIter',90000,'Display','off');
opts = optimoptions(@fmincon,'Algorithm','interior-point','MaxFunEvals',60000,'MaxIter',90000);

KK=0:1:K;
Z=zeros(1,1);
start=[zeros(1,2) repmat([0.5 1],1,(M-1))];

iter=2;
ll(1)=10^6;ll(2)=0;
Phiopt=Phi;

B1=-eps*ones(2*M -3,1);
while (iter <=length(KK) ) 
    p1=Phi(1:KK(iter),:);
    p2=Phi(K+1:K+KK(iter),:);
    p=[p1;p2];
    m_temp=size(p,1);
    lb=[-Inf*ones(1,m_temp) 10^(-6)*n^(-.25)*ones(1,2*(M-1))];
    %lb=[-Inf*ones(1,m_temp) repmat([10^(-6) 0.1],1,M-1)];
if (M>1)
    ub=[Inf*ones(1,m_temp) 1 10*sqrt(n)*ones(1,2*M-3)];
else
    ub=Inf*ones(1,m_temp);
end
A1=zeros(2*M -3,m_temp+ 2*(M-1));
for i=1:(2*M -3)
    A1(i,:)=[zeros(1,m_temp+i-1) (-1)^(i-1) (-1)^i zeros(1,2*M - i -3)];
end
    problem=createOptimProblem('fmincon','objective',@(c)FormwtedMmodalLikelihoodFromC(c,X,Y,y0,T,p,t1,M,h),'x0',start,'Aineq',A1,'bineq',B1,'lb',lb,'ub',ub,'options',opts);
          gs = GlobalSearch('Display','off');
         [temp,tempval]=run(gs,problem);
  %   temp=fmincon(@(c)FormwtedMmodalLikelihoodFromC(c,X,Y,y0,T,p,t1,M,h),start,A1,B1,[],[],lb,ub,[],options);
   %  tempval=FormwtedMmodalLikelihoodFromC(temp,X,Y,y0,T,p,t1,M,h);
      dcell{iter}=temp;
    ll(iter)=tempval;
    dd=dcell{iter};
    sd=length(dd(1:m_temp));
    d1=dd(1:sd/2);
    d2=dd((sd/2)+1:m_temp);
    clear start;clear temp;clear tempval;
    start=[d1,Z,d2,Z repmat([0.5 1],1,(M-1))];
    
    iter=iter+1;
end
 [~,i0]=min(ll(2:end));
i0=i0+1;
p1=Phi(1:KK(i0),:);
p2=Phi(K+1:K+KK(i0),:);
Phiopt=[p1;p2];
else %choice=2;Meyer basis set
    options = optimset('MaxFunEvals',60000,'MaxIter',90000,'Display','off');
    J=[1,2,2,3,3,4,4];
    K=[4,4,5,5,6,6,7];
    start(1,:)=[ones(1,K(1)+1)*eps repmat([0.5 1],1,(M-1))];
    start(2,:)=[ones(1,K(1)+1)*eps repmat([10^(-5) 10^(-3)],1,(M-1))];
    start(3,:)=[ones(1,K(1)+1)*eps repmat([5 10],1,(M-1))];
    %start(4,:)=[zeros(1,2) repmat([5 10],1,(M-2)) 10^(-5) 10^(-3)];
    %
    iter=1;
    ll(1)=0;
    
    B1=-1*ones(2*M -3,1);
    while (iter <=length(K) )
        
        p=meyerbasisgenerator(J(iter),K(iter),t1);
        m_temp=size(p,1);
        lb=[-Inf*ones(1,m_temp) 10^(-6)*n^(-.25)*ones(1,2*(M-1))];
        if (M>1)
            ub=[Inf*ones(1,m_temp) 1 10*sqrt(n)*ones(1,2*M-3)];
        else
            ub=Inf*ones(1,m_temp);
        end
        A1=zeros(2*M -3,m_temp+ 2*(M-1));
        for i=1:(2*M -3)
            A1(i,:)=[zeros(1,m_temp+i-1) (-1)^(i-1) (-1)^i zeros(1,2*M - i -3)];
        end
        for startiter=1:3
            temp(startiter,:)=fmincon(@(c)FormwtedMmodalLikelihoodFromC(c,X,Y,y0,T,p,t1,M,h),start(startiter,:),A1,B1,[],[],lb,ub,[],options);
            tempval(startiter)=FormwtedMmodalLikelihoodFromC(temp(startiter,:),X,Y,y0,T,p,t1,M,h);
        end
        [~,strtpt]=min(tempval);
        dcell{iter}=temp(strtpt,:);
        ll(iter)=FormwtedMmodalLikelihoodFromC(dcell{iter},X,Y,y0,T,p,t1,M,h)
        dd=dcell{iter};
        dd=dd(1:(m_temp-1));% Leave out the scaling function
        clear start;clear temp;clear tempval;
        
        if (iter<length(K))
            if (K(iter+1)>K(iter))
                index=[];
                for indloop=0:((J(iter)-1))
                    index=[index (1+ indloop*K(iter+1)):(indloop*K(iter+1)+K(iter))];
                end
                Z=zeros(1,1+J(iter+1)*K(iter+1));
                Z(index)=dd;
                start(1,:)=[Z repmat([0.5 1],1,(M-1))];% 3 different starting points considered.. can be skipped if using GlobalSearch anyway
                start(2,:)=[Z repmat([10^(-5) 10^(-3)],1,(M-1))];
                start(3,:)=[Z repmat([5 10],1,(M-1))];
                iter=iter+1;
            else
                Z=zeros(1,1+J(iter+1)*K(iter+1));
                Z(1:length(dd))=dd;
                start(1,:)=[Z repmat([0.5 1],1,(M-1))];
                start(2,:)=[Z repmat([10^(-5) 10^(-3)],1,(M-1))];
                start(3,:)=[Z repmat([5 10],1,(M-1))];
                iter=iter+1;
            end
        else
            iter=iter+1;
        end
    end
    [~,i0]=min(ll);
    Phiopt=meyerbasisgenerator(J(i0),K(i0),t1);
end
%%
m=size(Phiopt,1);
d=dcell{i0};
gamEst = FormGammaFromC(d(1:m),Phiopt);
l=[1 d((m+1):end)];
g_ind=floor([(1:(100/(2*M)):100) 100]);
gval=ones(1,length(g_ind));
gval(1)=eps;
gval(end)=eps;
gval(2:(end-1))=l;

g=interp1(t1(g_ind),gval,t1,'pchip');
g=round(g.*(10^5))/(10^5);
fn = interp1(t1, g, (t1(end)-t1(1)).*gamEst + t1(1));
fn=fn/(sum(fn)/T); % Returns density on [0,1]
% fn=fn/(B-A);
% fn=fn/(sum(fn)*(B-A)/T);

