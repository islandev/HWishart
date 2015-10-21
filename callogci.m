filePath='C:\Users\Administrator\Desktop\data_polsar\AIRSAR_DeathValley\vge\C3\';
row=76;
col=74 ;
data = zeros(row,col,9);
%从文件夹读取数据
d=3;
L=4;
fIn = fopen([filePath 'C11.bin'],'r');
data(:,:,1) = fread(fIn,[col,row],'float').';     fclose(fIn);
fIn = fopen([filePath 'C22.bin'],'r');
data(:,:,2) = fread(fIn,[col,row],'float').';     fclose(fIn);
fIn = fopen([filePath 'C33.bin'],'r');
data(:,:,3) = fread(fIn,[col,row],'float').';     fclose(fIn);
fIn = fopen([filePath 'C12_real.bin'],'r');
data(:,:,4) = fread(fIn,[col,row],'float').';   fclose(fIn);
fIn = fopen([filePath 'C13_real.bin'],'r');
data(:,:,5) = fread(fIn,[col,row],'float').';   fclose(fIn);
fIn = fopen([filePath 'C23_real.bin'],'r');
data(:,:,6) = fread(fIn,[col,row],'float').';   fclose(fIn);
fIn = fopen([filePath 'C12_imag.bin'],'r');
data(:,:,7) = fread(fIn,[col,row],'float').';   fclose(fIn);
fIn = fopen([filePath 'C13_imag.bin'],'r');
data(:,:,8) = fread(fIn,[col,row],'float').';   fclose(fIn);
fIn = fopen([filePath 'C23_imag.bin'],'r');
data(:,:,9) = fread(fIn,[col,row],'float').';   fclose(fIn); 
%对矩阵进行构造
c11=data(:,:,1); %c11
c22=data(:,:,2); %c22
c33=data(:,:,3); %c33

c12=data(:,:,4)+i*data(:,:,7); %c12
c13=data(:,:,5)+1i*data(:,:,8); %c13
c23=data(:,:,6)+i*data(:,:,9); %c23
c21=data(:,:,4)-i*data(:,:,7); %c12
c31=data(:,:,5)-i*data(:,:,8); %c13
c32=data(:,:,6)-i*data(:,:,9); %c23
num=row*col;
%求基于矩阵的对数累积量
m1=0;
m2=0;
m3=0;
m4=0;
m5=0;
for s=1:row
    for j=1:col
        z_matrix=[c11(s,j),c12(s,j),c13(s,j);c21(s,j),c22(s,j),c23(s,j);c31(s,j),c32(s,j),c33(s,j)];
       
        m1=log(det(z_matrix))+m1;
        m2=log(det(z_matrix))^2+m2;
        m3=log(det(z_matrix))^3+m3;
        m4=log(det(z_matrix))^4+m4;
        m5=log(det(z_matrix))^5+m5;
        
    end
end

m1=m1/num
m2=m2/num
m3=m3/num
m4=m4/num
m5=m5/num

%计算H-Wishart的对数累积量


c2 =real( m2 - m1^2)-sumpsi(d,L,1)
c3 = real(m3 - 3*m1*m2 + 2*m1^3)-sumpsi(d,L,2)
c4=real(m4-4*m1*m3-3*m2^2+12*m1^2-6*m1^4)-sumpsi(d,L,3)
c5=real(m5-5*m1*m4-10*m2*m3+20*m1^2*m3+30*m1*m2^2-60*m1^3*m2+24*m1^5)-sumpsi(d,L,4)
