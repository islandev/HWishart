function [ ksd ] = ksdistance( c1,c2 )
%�����������ݵ�ks����
%���� �������� ����һ�� sumΪ1
%���Ϊ���������ks����


ksd=0;
%������ڵ�cdf
cur_val1=0;
cur_val2=0;
%��ȡ������ ��cdf �����������ݵ�ks����
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

