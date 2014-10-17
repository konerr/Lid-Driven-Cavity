function[]=ppe(n, h, dt, u_star, v_star)
global p
% The following code solves a 2-D Poisson Equation
% Author: Rahul Babu Koneru

%% Input Parametres
iter=8000;
tol=0.01;

% For SOR Only
% psi=(1-2*(sin(pi/2*n)));
% omg=(2/(1+(1-psi^2)^0.5));
% omg = 2/(1+(pi/h));
% omg=1.6;

%% Initial condition
p=zeros(n,n);
f=zeros(n,n);
    
%% Source Terms
j = 2:n-1;
k = 2:n-1;

f(j,k)=(h/(dt))*(u_star(j+1,k) - u_star(j,k) + v_star(j,k+1) - v_star(j,k));
  
%% Main Loop %%

i = 2:n-1;
j = 2:n-1;
for k=1:iter
 p_old=p;  
 
 p(i,j)=0.25 * (p_old(i+1,j) + p(i-1,j) + p_old(i,j+1) + p(i,j-1) - f(i,j)); %Gauss-Seidel%
%         p(i,j)=0.25*(p_old(i+1,j)+p_old(i-1,j)+p_old(i,j+1)+p_old(i,j-1)-f(i,j)*h^2); %Jacobi%
%         p(i,j)=(1-omg)*p_old(i,j)+omg*(0.25*(p_old(i+1,j)+p(i-1,j)+p_old(i,j+1)+p(i,j-1)-f(i,j)*h^2)); %SOR%

%Pressure Boundary Condition
p(1,(2:n-1)) = p(2,(2:n-1));
p(n,(2:n-1)) = p((n-1),(2:n-1));
p((2:n-1),1) = p((2:n-1),2);
p((2:n-1),n) = p((2:n-1),(n-1));
p(1,1) = p(1,2);
p(n,n) = p(n,n-1);
p(1,n) = p(1,n-1);
p(n,1) = p(n-1,1);

 %% Error
 err=norm(p(i,j)-p_old(i,j),inf);
%  fprintf('\n%d   %f\n',k,err);
 
 if err<=tol
     break
 end
end

fprintf('\nError after %d iterations: %f',k,err);
end
