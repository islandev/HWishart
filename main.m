clear all;clc

% fid=fopen('Ocean_C11.bin','r');
% [data,N]=fread(fid,'float');
% fclose(fid);

% save  Ocean_SanFrancisco data;

load Ocean_SanFrancisco;
k=1;
data=k*data;
N=numel(data);

L=4;
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
[  k,v,sigma,Sigma_H ] = ParameterEstimationhh( data_trun,L,d );

%---------------------------------------------------------------------
%fittiing
%---------------------------------------------------------------------
setting.k=k;
setting.v=v;
setting.L=L;
setting.d=d;
setting.sigma=sigma;
setting.Sigma_H=Sigma_H;

x=linspace(1e-6,1*maxsample,1e3);
y=HWpdf(x,setting);

%---------------------------------------------------------------------
%Plot
%---------------------------------------------------------------------

Nbin=128;
[nelements,xcenters]=hist(data_trun,Nbin);

figure;
bar(xcenters,nelements/N_trun/(xcenters(2)-xcenters(1)),'y');
hold on;

plot(x,y,'b-x');










