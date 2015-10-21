function f=myfun
syms x; 
syms y;
syms a;
syms b;
f=[psi(1,x)/y^2+psi(1,a)/b^2;psi(2,x)/y^3+psi(2,a)/b^3;psi(3,x)/y^4+psi(3,a)/b^4;psi(4,x)/y^5+psi(4,a)/b^5];
x0=[  5.0485; 0.8036;-2.3059e+004; -106.1066]
end