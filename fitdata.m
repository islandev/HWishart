function [ setting ] = fitdata( data ,L)
%对数据进行拟合 参数估计等
%   输入为 数据
%输出为 一个结构体里面包含有 参数和方差
k=1;
data=k*data;
N=numel(data);

    
d=1;

%---------------------------------------------------
% data truncation
%---------------------------------------------------
data_des=sort(data',2);
maxsample=max(data_des);

N_trun=floor(1*N);
data_trun=data_des(1:N_trun);

maxsample=max(data_trun,[],2);


%---------------------------------------------------------------------
%parameters estimation
%---------------------------------------------------------------------
[  k,v,sigma,Sigma_H,k_sigma ,k_cov] = ParameterEstimationhh( data_trun,L,d );

%---------------------------------------------------------------------
%fittiing
%---------------------------------------------------------------------
setting.k=k;
setting.v=v;
setting.L=L;
setting.d=d;
setting.sigma=sigma;
setting.Sigma_H=Sigma_H;
setting.k_sigma=k_sigma;
x=linspace(1e-6,1*maxsample,1e3);
%kwishart的pdf
y_k=kwishartpdf( k_sigma,k_cov,L,d ,maxsample);
y=HWpdf(x,setting);

%---------------------------------------------------------------------
%Plot
%---------------------------------------------------------------------

Nbin=1000;
[nelements,xcenters]=hist(data_trun,Nbin);
out=hist(data_trun,Nbin);
out=out/sum(out);

figure;
bar(xcenters,nelements/N_trun/(xcenters(2)-xcenters(1)),'y');
hold on;
%H-Wishart
plot(x,y,'b-x');
hold on;
%k-wishart
plot(x,y_k,'r')
xlabel('Texture'); %标记横坐标
ylabel('Pd'); %标记纵坐标
legend( 'data hist-diagram ','H-Wishart','K-distribution')
y=y.*maxsample/Nbin;

setting.ksd=ksdistance(out,y);
%计算MSE
y_k=y_k.*maxsample/Nbin;
k_mse=calMSE(out,y_k);
h_mse=calMSE(out,y);
setting.kmse=k_mse;
setting.hmse=h_mse;

end

