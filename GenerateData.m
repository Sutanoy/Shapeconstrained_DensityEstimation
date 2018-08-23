
%function [X,a0,b0,t,g0] = GenerateData(n)
function [X,t,f0] = GenerateData(n)
%%%%%% generate from mixture beta 1 %%%%%%%
% u=rand(1,n);
% X=betarnd(1,3,1,n).*(u<1/3)+ betarnd(1,4,1,n).*(u>1/3 & u<2/3) +betarnd(3,47,1,n).*(u>2/3);
% 
% 
% N = 100;
% t = (1:N)/N;
% 
% f0 = (betapdf(t,1,3) + betapdf(t,1,4) + betapdf(t,3,47))/3;
%%%%%% generate from mixture beta 2 %%%%%%%
% u=rand(1,n);
% X=betarnd(1,3,1,n).*(u<1/3)+ betarnd(1,4,1,n).*(u>1/3 & u<2/3) +betarnd(47,3,1,n).*(u>2/3);
% 
% 
% % N = 100;
% % t = (1:N)/N;
%  t = -1:0.001:2;
%  dt=0.001;
% %f0 = (betapdf(t,1,3) + betapdf(t,1,4) + betapdf(t,47,3))/3;
% f0=betapdf(t,10,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Generate from mixture normal %%%%%%%%%%%
%   t = -15:0.001:15;
%   dt=0.001;
% % % x=sqrt(1/10);
% % % % %% CLAW DENSITY %%INCLUDE last 3 lines
% % %f0 = 0.5*normpdf(t,0,1) + x*normpdf(t,-1,x^2)+ x*normpdf(t,0.5-1,x^2)+ x*normpdf(t,1-1,x^2)+ x*normpdf(t,1.5-1,x^2)+ x*normpdf(t,2-1,x^2);
% % %% UNIMODAL DENSITIES
% % %f0=normpdf(t,0,1);
% % % % % %% kurtotic unimodal density %%
% % %f0=(2/3)*normpdf(t,0,0.5) + (1/3)*normpdf(t,0,10);
% % %f0=(4/5)*normpdf(t,0,4) + (1/5)*normpdf(t,0,.5);
% % %%
% % %Gamma densities
% % % t=0:.01:100;
% % % dt=0.01;
% % % f0=gampdf(t,30,0.5);
% % % % % %%
% % %  f0=normpdf(t,-1,0.6)*0.5 + normpdf(t,1,0.6)*0.5;
% %  % f0=normpdf(t,-5,1)*.2 + normpdf(t,-3,0.3)*.3 + normpdf(t,0,0.5)/4 + normpdf(t,2,0.5)/4;
% % % f0=normpdf(t,-1,0.25)/3 + normpdf(t,2,0.3)/3 + normpdf(t,0,0.25)/3;
% % %  f0=normpdf(t,-1,1)/3 + 2*normpdf(t,1,0.3)/3;
% % %f0=normpdf(t,0,0.5)*0.95 + normpdf(t,3,1)*0.05;
% % %% flat density%%
%  f0=0.*(t<=0) + 0.*(t>=1) + t.*(t>0 & t<=1/3) + (1/3).*(t>1/3 & t<2/3) +(1-t).*(t<1 & t>=2/3);
% % %f0=normpdf(t,1/2,1);
% % % % % 
%   f0 = f0/(sum(f0)*dt);
%   f = f0;
% % % 
%   F = cumsum(f)*dt; 
% % % 
% for i = 1:n
%     u = rand;
%     X(i) = t(min(find(F > u)));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Generate from density in Tokdar paper%%%%
 N = 100;
 t = (0:N)/N;
 dt=1/N;
% t = -2:0.01:2;
%   dt=0.001;
% x=0.2;
% f0 =0.75*3*exp((-3.*t)) + 0.25*sqrt(32/pi)*exp(-32.*((t-.75).^2)); 
% f0 =0.75*normpdf(t,0.3,0.2).*(t>=0 & t<=1)+ 0.25*sqrt(32/pi)*exp(-32.*((t-.75).^2)).*(t>=0 & t<=1); 
 f0 =0.75*normpdf(t,0.3,0.2).*(t>=0 & t<=1)+ 0.25*normpdf(t,0.75,1/8).*(t>=0 & t<=1);
 % f0=normpdf(t,0.1,0.2).*(t<0.3)/2 + normpdf(t,0.7,0.2).*(t>.6)/2;
% f0=normpdf(t,0,0.4).*(t>=0 &t<=1);
f0 = f0/(sum(f0)*dt);

f = f0;

F = cumsum(f)*dt; 
for i = 1:n
    u = rand;
    X(i) = t(min(find(F > u)));
 end
