function y=HWpdf(x,setting,res)

%--------------------------------------------
% Preprocessing
%--------------------------------------------
k=setting.k;
v=setting.v;
L=setting.L;
d=setting.d;
sigma=setting.sigma;
Sigma_H=setting.Sigma_H;

N=numel(x);

if nargin<3
    res=1e-7;
end
%--------------------------------------------
% Computing the p.d.f. of x
%--------------------------------------------

for i=1:N
    f=@(t) t.^(k.*v-L.*d-1).*exp(-L./t./Sigma_H.*x(i)-(t./sigma).^v)  ;
    I(i)=quadgk(f,0,inf,'RelTol',res);
end

y=abs(v).*L.^(L.*d).*x.^(L-d)./sigma.^(k.*v)./gamma(k)./gamma(L)./Sigma_H.^L.*I;



end
