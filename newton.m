function z=newton (f,x0,e)
syms x y a b,
Jacobi=jacobian(f,[x,y,a,b]);
x=x0(1);
y=x0(2);
a=x0(3);
b=x0(4)
n=1;
x1=x0-inv(eval(Jacobi))*eval(f);
while norm(x1-x0)>e & n<1000
x=x0(1);
y=x0(2);
a=x0(3);
b=x0(4); 
n=n+1;
x1=x0-inv(eval(Jacobi))*eval(f);
end
n
z=x1
