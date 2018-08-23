function L = FormMmodalLikelihoodFromC(c,X,N,p,t1,M)

%%%%%%%%%%%%%%%Unpenalized likelihood%%%%%%%%%%%%%%%%%%%%%%%%%%
l0=length(c);
l=[1 c(((size(p,1))+1):end)];
%%
%l=[1 1];% flat spot
%%
c=c(1:size(p,1));
gam0 = FormGammaFromC(c,p);
gam = (gam0-gam0(1))/(gam0(end)-gam0(1));

g_ind=floor([(1:(100/(2*M)):100) 100]);
%g_ind=[1 40 60 100];% for flat spots
gval=ones(1,length(g_ind));
gval(1)=eps;
gval(end)=eps;
gval(2:(end-1))=l;
fp=interp1(t1(g_ind),gval,t1,'pchip');
fp=round(fp.*(10^5))/(10^5);
fn = interp1(t1, fp, (t1(end)-t1(1)).*gam + t1(1),'linear','extrap') ;
fn=fn/(sum(fn)/N);
yy=interp1(t1,fn,X,'linear','extrap');
 L = -2*sum(log(yy))+2*l0; %AIC


