function [ mse_val ] = calMSE( c1,c2 )
%�������������MSE
%   Detailed explanation goes here
%����Ϊ�����ȳ�����
%���Ϊmse
    m_size=size(c1,2);
%ʹ��ѭ������MSE    
    mse_val=0;
    for i=1:m_size
        mse_val=mse_val+(c1(i)-c2(i))^2;
    end
    mse_val=(mse_val/m_size)^1/2;
end

