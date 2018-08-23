%%We assume boundary values to be zero in this code.
n=100;

%Y=rand(1,n);%predictor

 Y=normrnd(0,1,1,n);% Y is the predictor, in this code.
% u=rand(1,n);
% X=normrnd(Y-1.5,0.5).*(u<0.5) + normrnd(Y+1.5,0.5).*(u>=0.5);%response

%X=normrnd(Y,1);
X=(2*Y -1).^1+laprnd(1,n,0,1); % X is the response 
M=1;

 
y0=norminv(0.5,0,1);
[test1,test2,hh]=ksdensity(Y,y0);
h=hh/sqrt(test1(1));
diff=Y-ones(1,length(X))*y0;
normdiff=sqrt(sum(diff.^2,1));
wt=normpdf(normdiff/h,0,1);% weights for the weighted likelihood fn.
ind=find(wt>quantile(wt,0.5)); %nearest 50\% of the observations

wt(ind)=wt(ind)/sum(wt(ind));

A=min(X(ind))- (std(X(ind))/sqrt(length(ind)));
B=max(X(ind))+ (std(X(ind))/sqrt(length(ind)));
   
T=100;
t0=0:(1/(T-1)):1;
t1=t0*(B-A) +A;

fn=Modalconstrainedcde(X(ind),Y(ind),y0,A,B,M,h,1);
fn=fn/(B-A);
%ft=0.5*normpdf(t1,y0-1.5,0.5) + 0.5*normpdf(t1,y0+1.5,0.5);
 ft=exp(-abs(t1-((2*y0 -1)^1))/(1/sqrt(2)))/(sqrt(2));
%ft=normpdf(t1,y0,1);

figure(1);
clf;
plot(t1,fn,'Linewidth',3);
set(gca,'fontsize',18);
xlim([min(grid(i0,:))-1  max(grid(i0,:))+1]);
xlabel('X=0'),ylabel('conditional density');
