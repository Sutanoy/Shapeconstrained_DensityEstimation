function L = FormMmodalLikelihoodboundaryFromC(c,X,N,p,t1,M)

% c is the parameter vector
% X is the SCALED observations.
% N : number of grid points(usually 100)
% p : Basis matrix
% t1 : the grid on [0,1], the parameter N can be removed as length(t1)
% M : number of modes.
% This code assumes antimodes on the boundaries but not restricted to be
% zero at the boundaries



%%%%%%%%%%%%%%%Penalized likelihood%%%%%%%%%%%%%%%%%%%%%%%%%%
l0=length(c); 
l=[1 c(((size(p,1))+1):(end-2))]; % interior critical point heights
g_ind=floor([(1:(100/(2*M)):100) 100]);% location of critical points 
                                       %(chosen to be equally spaced)
gval=ones(1,length(g_ind)); % initialization
gval(1)=c(end); % last value is the left boundary value
                %( helped me format the constraints in fmincon, satisfies
                %value<1 if antimode on boundary.)
gval(end)=c(end-1); % 2nd last parameter is right boundary.
                    %(< c(end -2) if antimode on right boundary)
gval(2:(end-1))=l; % interior values
c1=c(1:size(p,1));% basis coefficients of \gamma
gam0 = FormGammaFromC(c1,p); %forms gamma given basis matrix & coefficients
gam = (gam0-gam0(1))/(gam0(end)-gam0(1));%removes numerical errors.

fp=interp1(t1(g_ind),gval,t1,'pchip'); %the initial template function
fp=round(fp.*(10^5))/(10^5); % rounded upto 5 decimal places
fn = interp1(t1, fp, (t1(end)-t1(1)).*gam + t1(1),'linear','extrap') ;
fn=fn/(sum(fn)/N); % the density estimate
yy=interp1(t1,fn,X,'linear','extrap');
 L = -2*sum(log(yy))+2*l0; % penalized likelihood function


