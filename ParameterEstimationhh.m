function [ k,v,sigma,Sigma_H ,k_sigma,k_cov] = ParameterEstimationhh( x,L,d )

[row,col]=size(x);
N=numel(x);

% log moments and log cumulants
m1=ones(1,row)*log(x)*ones(col,1)./N;
m2=ones(1,row)*(log(x).^2)*ones(col,1)./N;
m3=ones(1,row)*(log(x).^3)*ones(col,1)./N;

k1_HW =m1;
k2_HW =m2 - m1^2;
k3_HW =m3 - 3*m1*m2 + 2*m1^3;

% The log cumulants of Generalized Gamma distribution
k2_GG = k2_HW-psi(1,L)  ;
k3_GG = k3_HW-psi(2,L)  ;
%-----------------------------------------------------
% parameters evaluation of Generalized Gamma distribution

%  evaluation of kappa
m=k2_GG.^3./k3_GG.^2;
f=@(y) psi(1,abs(y)).^3./(psi(2,abs(y)).^2)-m;
k=abs(fzero(f,10));

%  evaluation of nu
v_norm= (psi(1,k)/k2_GG).^(1/2) ;
if k3_GG>0
    v=v_norm;
else
    v=-v_norm;
end
%计算kwishart的参数 k——sigama
f_k=@(k_sigma)norm(psi(1,k_sigma)-k2_GG).^2;
[k_sigma,err]=fminsearch(f_k,0);
 

%  evaluation of sigma
%  under the condition of unity of texture
sigma=gamma(k)/gamma(k+(1/v));

%  evaluation of Sigma_H
Sigma_H=1/sigma*L*exp(m1-psi(0,k)/v-psi(0,L));
%计算k-Wishart的方差
k_cov=exp(m1-psi(0,k_sigma)-psi(0,L)+log(k_sigma*L))





end



