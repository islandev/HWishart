function [ ksd ] = ksdistance( c1,c2 )
%计算两组数据的ks距离
%输入 两个数组 长度一致 sum为1
%输出为两个数组的ks距离


ksd=0;
%存放现在的cdf
cur_val1=0;
cur_val2=0;
%读取数组求 其cdf ，和数组数据的ks距离
for i=1:size(c1,2)
   cur_val1=cur_val1+c1(i);
   cdf1(i)=cur_val1;
  
   cur_val2=cur_val2+c2(i);
    cdf2(i)=cur_val2;
   cur_ksd=abs(cur_val1-cur_val2);
   if ksd<=cur_ksd
       ksd=cur_ksd;
   else 
       ksd=ksd;
   end
end

    
end

