clear all;clc;

filePath = 'C:\Users\Administrator\Desktop\data_polsar\AIRSAR_SanFrancisco\sea\C3\';

row=127;
col=127;
[k,sigma,v,covhh,covhv,covvv]=ParameterEstimation(filePath,row,col,4,1);
%[k,sigma,v,covhh]=ParameterEstimationhh(filePath,row,col); 




%[k,sigma,v]= mlcparameter( filePath,row,col ) %“Ï÷ Œ∆¿Ì