function [ pdfk ] = kwishartpdf( sigma,c,L,d ,maxsample)
%计算KWishart的PDF
%   Detailed explanation goes here
%sigma  is the shape parameter 
%c  是方差
%L  the look of the sar img
% d  the dimision of the sar img




L=4;
d=1;
y=[];                      
j=1;
x=linspace(1e-6,1*maxsample,1e3);
for i=1:1000
    spa=2*x(i)^(L-d)*(L*sigma)^((sigma+L*d)/2);
    spb=c^L*gamma(sigma)*gamma(L);
    spc=(x(i)/c)^((sigma-L*d)/2);
    spd=besselk(sigma-L*d,2*((sigma*L*x(i)/c)^0.5));
    y(i)=spa*spc*spd/spb;
   
end



pdfk=y;




end

