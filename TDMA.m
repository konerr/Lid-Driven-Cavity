function[]=TDMA(n, a, b, c, f)
%Input arguments
%n=no.of grid points
%a,b,c=vectors of principal,lower and upper diagonal 
%elements respectively

global x
%% LU Decomposition
    alp(1)=a(1);
    gam(1)=c(1)/alp(1);
    for i=2:n
        alp(i)=a(i)-b(i-1)*gam(i-1);
        gam(i)=c(i-1)/alp(i);
    end
%% U*x=z
    z(1)=f(1)/alp(1);
    for j=2:n
        z(j)=(f(j)-b(j-1)*z(j-1))/alp(j);
    end
%% Back Substitution 
    x(n)=z(n);
    for k=n-1:-1:1
        x(k)=z(k)-gam(k)*x(k+1);
    end
%     tx=x;
end
