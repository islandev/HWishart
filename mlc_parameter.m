clear all;clc;

filePath = 'C:\Users\Administrator\Desktop\data_polsar\ESAR_Oberpfaffenhofen\forest\C3\';

row=126;
col=121 ;
%���ú��� ����ھ���Ķ����ۻ����Ĳ���
[k,sigma,v,covhh,covhv,covvv]=ParameterEstimation(filePath,row,col); %��ͨ��
%[k,sigma,v,covhh]=ParameterEstimationhh(filePath,row,col); 



%[k,sigma,v]= mlcparameter( filePath,row,col ) %��������

    




