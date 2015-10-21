function [ mse_val ] = calMSE( c1,c2 )
%计算两个数组的MSE
%   Detailed explanation goes here
%输入为两个等长数组
%输出为mse
    m_size=size(c1,2);
%使用循环计算MSE    
    mse_val=0;
    for i=1:m_size
        mse_val=mse_val+(c1(i)-c2(i))^2;
    end
    mse_val=(mse_val/m_size)^1/2;
end

