%The following code solve a 2D inompressible Navier Stokes equation for   %
%a driven cavity problem. Non-linear convection terms are discretized     %
%using a 2nd order AB scheme while the diffusion terms are treated        %
%implicitly using Crank Nicolson scheme. Alternate Direction Implicit(ADI)%
%is used to handle the implicit part.%
%TDMA(n,a,b,c,f) is the tri-diagonal solver(Thomas Algorithm)             %
%ppe(n,h,dt,u_star,v_star) is the Poisson equation solver which uses Gauss%
%Seidel iterative technique.                                              %
%u_frac, v_frac are the intermediate velocities in the ADI scheme         %
%u_star, v_star are the fractional step velocities from operator splitting%

clc
clear all
close all
tic

%% Initialization
ax = 0;                   %Initial point on the grid
bx = 1;                   %End point on the grid
n = 81;                   %Grid size
dt = 0.0005;              %Time step
h = (bx-ax)/(n-1);        %Grid spacing
X = ax:h:bx;
Y = X;
Re = 1000;                %Reynold's number
[lx,ly] = meshgrid(X,Y);  %2D Mesh

a = zeros(1,n-2);
b = zeros(1,n-3);
c = zeros(1,n-3);
f = zeros(1,n-2);
ae_I = zeros(1,n-4);
aw_I = zeros(1,n-4);
f1_I = zeros(1,n-4);
an_I = zeros(1,n-4);
as_I = zeros(1,n-4);
f2_I = zeros(1,n-4);
u_frac1 = zeros(n-2,n-2);
v_frac1 = zeros(n-2,n-2);
u_star1 = zeros(n-2,n-2);
v_star1 = zeros(n-2,n-2);
u_new1 = zeros(n-2,n-2);
v_new1 = zeros(n-2,n-2);
err = u_new1;

global x                %Global variable output from TDMA solver
global p                %Global declaration of output from ppe solver

%Initial condition
u = zeros(n,n);         %Velocity in X-directiom
v = zeros(n,n);         %Velocity in Y-directiom

%Boundary Conditions
u(1,:) = 0;             %Left wall
v(1,:) = 0;
u(n,:) = 0;             %Right wall
v(n,:) = 0;
u(:,n) = 1;             %Top wall
v(:,n) = 0;
u(:,1) = 0;             %Bottom wall
v(:,1) = 0;

u_frac = u;
u_old = u;
v_frac = v;
v_old = v;
u_star = u;
v_star = v;
u_new = u;
v_new = v;

time = 5;

%% 
for i = dt:dt:time
%% ADI
    
%-----------------------------------------------------------------------%    
%---------------------------Velocity in X-Direction---------------------%
%-----------------------------------------------------------------------%

%----------------------------------X-Sweep------------------------------%
    for k = 2:n-1

        for j = 2:n-1
            if (j==2)
                ap_n = 1+(dt/(Re*h^2));
                as_n = dt/(2*Re*h^2);         
                an = 0;
                f2_n = u(j,k) +(dt/(2*Re*h^2)) * u(j-1,k)... 
                      + (1/2) * (dt/(2*h)) * (((v_old(j+1,k)...
                      + v_old(j,k))*(u_old(j+1,k)+u_old(j,k))...
                      - (v_old(j-1,k)+v_old(j,k))*(u_old(j-1,k)+u_old(j,k)))/2)...
                      - (3/2) * (dt/(2*h)) * (((v(j+1,k)+v(j,k))...
                      *(u(j+1,k)+u(j,k)) - (v(j-1,k)+v(j,k))*(u(j-1,k)+u(j,k)))/2)...
                      + (dt/(2*Re*h^2)) * (u(j,k+1) - 2*u(j,k) + u(j,k-1));            

            elseif (j==n-1)
                ap_s = 1+(dt/(Re*h^2));
                as = 0;
                an_s = dt/(2*Re*h^2);      
                f2_s = u(j,k) +(dt/(2*Re*h^2)) * u(j+1,k)... 
                          + (1/2) * (dt/(2*h)) * (((v_old(j+1,k)...
                          + v_old(j,k))*(u_old(j+1,k)+u_old(j,k))...
                          - (v_old(j-1,k)+v_old(j,k))*(u_old(j-1,k)...
                          + u_old(j,k)))/2)...
                          - (3/2) * (dt/(2*h)) * (((v(j+1,k)+v(j,k))...
                          *(u(j+1,k)+u(j,k)) - (v(j-1,k)+v(j,k))*(u(j-1,k)+u(j,k)))/2)...
                          + (dt/(2*Re*h^2)) * (u(j,k+1) - 2*u(j,k) + u(j,k-1));

            elseif (j~=2 && j~=n-1)
                ap = 1+(dt/(Re*h^2));
                an_I(1,j-2) = dt/(2*Re*h^2);         
                as_I(1,j-2) = dt/(2*Re*h^2);
                f2_I(1,j-2) = u(j,k)... 
                          + (1/2) * (dt/(2*h)) * (((v_old(j+1,k)+v_old(j,k))...
                          * (u_old(j+1,k)+u_old(j,k))...
                          - (v_old(j-1,k)+v_old(j,k))*(u_old(j-1,k)...
                          + u_old(j,k)))/2)...
                          - (3/2) * (dt/(2*h)) * (((v(j+1,k)+v(j,k))...
                          * (u(j+1,k)+u(j,k)) - (v(j-1,k)+v(j,k))*(u(j-1,k)+u(j,k)))/2)...
                          + (dt/(2*Re*h^2)) * (u(j,k+1) - 2*u(j,k) + u(j,k-1));                        
            end
        end    


        a(1,:) = ap;
        c = -[an_I(1,:) an_s];
        b = -[as_n as_I(1,:)];
        f = [f2_n f2_I(1,:) f2_s];
        TDMA(n-2,a,b,c,f)
        u_frac1(:,k-1)=x;

    end
    
    for j = 2:n-1
            for k = 2:n-1
                u_frac(j,k) = u_frac1(j-1,k-1);
            end
    end

%---------------------------------Y-Sweep--------------------------------%
    for j = 2:n-1

        for k = 2:n-1
            if (k==2)
                ap_w = 1+(dt/(Re*h^2));
                ae_w = dt/(2*Re*h^2);         
                aw = 0;
                f1_w = u_frac(j,k) +(dt/(2*Re*h^2)) * u(j,k-1)...
                       + (1/2) * (dt/(2*h)) * (((u_old(j,k+1) + u_old(j,k))^2 ...
                       - (u_old(j,k-1)+u_old(j,k))^2)/2)...
                       - (3/2) * (dt/(2*h)) * (((u(j,k+1) + u(j,k))^2 ...
                       - (u(j,k-1) + u(j,k))^2)/2)...
                       + (dt/(2*Re*h^2)) * (u_frac(j+1,k) - ...
                       2*u_frac(j,k) + u_frac(j-1,k));            

            elseif (k==n-1)
                ap_e = 1+(dt/(Re*h^2));
                ae = 0;
                aw_e = dt/(2*Re*h^2);      
                f1_e = u_frac(j,k) +(dt/(2*Re*h^2)) * u(j,k+1)...
                          + (1/2) * (dt/(2*h)) * (((u_old(j,k+1) + u_old(j,k))^2 ...
                          - (u_old(j,k-1)+u_old(j,k))^2)/2)...
                          - (3/2) * (dt/(2*h)) * (((u(j,k+1) + u(j,k))^2 ...
                          - (u(j,k-1) + u(j,k))^2)/2)...
                          + (dt/(2*Re*h^2)) * (u_frac(j+1,k) - ...
                          2*u_frac(j,k) + u_frac(j-1,k));

            elseif (k~=2 && k~=n-1)
                ap = 1+(dt/(Re*h^2));
                ae_I(1,k-2) = dt/(2*Re*h^2);         
                aw_I(1,k-2) = dt/(2*Re*h^2);
                f1_I(1,k-2) = u_frac(j,k)...
                          + (1/2) * (dt/(2*h)) * (((u_old(j,k+1) + u_old(j,k))^2 ...
                          - (u_old(j,k-1)+u_old(j,k))^2)/2)...
                          - (3/2) * (dt/(2*h)) * (((u(j,k+1) + u(j,k))^2 - ...
                          (u(j,k-1) + u(j,k))^2)/2)...
                           + (dt/(2*Re*h^2)) * (u_frac(j+1,k) - ...
                           2*u_frac(j,k) + u_frac(j-1,k));                        
            end
        end

        a(1,:) = ap;
        b = -[aw_I(1,:) aw_e];
        c = -[ae_w ae_I(1,:)];
        f = [f1_w f1_I(1,:) f1_e];
        TDMA(n-2,a,b,c,f)
        u_star1(j-1,:) = x;

    end
    
    for j = 2:n-1
            for k = 2:n-1
                u_star(j,k) = u_star1(j-1,k-1);
            end
    end

%------------------------------------------------------------------------%    
%----------------------Velocity in Y-Direction---------------------------%
%------------------------------------------------------------------------%

%------------------------------X-Sweep-----------------------------------%
    for k = 2:n-1

        for j = 2:n-1

            if (j==2)
                ap_n = 1+(dt/(Re*h^2));
                as_n = dt/(2*Re*h^2);         
                an = 0;
                f2_n = v(j,k) +(dt/(2*Re*h^2)) * v(j-1,k)...
                       + (1/2) * (dt/(2*h)) * (((v_old(j+1,k)...
                       + v_old(j,k))^2 - (v_old(j-1,k) + v_old(j,k))^2)/2)...
                       - (3/2) * (dt/(2*h)) * (((v(j+1,k) + v(j,k))^2 ...
                       - (v(j-1,k) + v(j,k) )^2)/2)...
                       + (dt/(2*Re*h^2)) * (v(j,k+1) - ...
                       2*v(j,k) + v(j,k-1));    

            elseif (j==n-1)
                ap_s = 1+(dt/(Re*h^2));
                as = 0;
                an_s = dt/(2*Re*h^2);      
                f2_s = v(j,k) +(dt/(2*Re*h^2)) * v(j+1,k)...
                       + (1/2) * (dt/(2*h)) * (((v_old(j+1,k) + v_old(j,k))^2 ...
                       - (v_old(j-1,k) + v_old(j,k))^2)/2)...
                       - (3/2) * (dt/(2*h)) * (((v(j+1,k) + v(j,k))^2 ...
                       - (v(j-1,k) + v(j,k) )^2)/2)...
                       + (dt/(2*Re*h^2)) * (v(j,k+1) - ...
                       2*v(j,k) + v(j,k-1));

            elseif (j~=2 && j~=n-1)
                ap = 1+(dt/(Re*h^2));
                an_I(1,j-2) = dt/(2*Re*h^2);         
                as_I(1,j-2) = dt/(2*Re*h^2);
                f2_I(1,j-2) = v(j,k)...
                          + (1/2) * (dt/(2*h)) * (((v_old(j+1,k) ...
                          + v_old(j,k))^2 - (v_old(j-1,k) + v_old(j,k))^2)/2)...
                          - (3/2) * (dt/(2*h)) * (((v(j+1,k) + v(j,k))^2 ...
                          - (v(j-1,k) + v(j,k) )^2)/2)...
                          + (dt/(2*Re*h^2)) * (v(j,k+1) - 2*v(j,k) + v(j,k-1));                        
            end
        end    

        a(1,:) = ap;
        c = -[an_I(1,:) an_s];
        b = -[as_n as_I(1,:)];
        f = [f2_n f2_I(1,:) f2_s];
        TDMA(n-2,a,b,c,f)
        v_frac1(:,k-1)=x;

    end
    
    for j = 2:n-1
            for k = 2:n-1
                v_frac(j,k) = v_frac1(j-1,k-1);
            end
    end

%-------------------------------Y-Sweep----------------------------------%
    for j = 2:n-1

        for k = 2:n-1
            if (k==2)
                ap_w = 1+(dt/(Re*h^2));
                ae_w = dt/(2*Re*h^2);         
                aw = 0;
                f1_w = v_frac(j,k) +(dt/(2*Re*h^2)) * v(j,k-1)...
                          + (1/2) * (dt/(2*h)) * ((u_old(j,k+1)...
                          + u_old(j,k))*(v_old(j,k+1)+v_old(j,k))...
                          - (u_old(j,k-1)+u_old(j,k-1))*(v_old(j,k-1)...
                          + v_old(j,k))/2)...
                          - (3/2) * (dt/(2*h)) * ((u(j,k+1)+u(j,k))*(v(j,k+1)+v(j,k))...
                          - (u(j,k-1)+u(j,k-1))*(v(j,k-1)+v(j,k))/2)...
                          + (dt/(2*Re*h^2)) * (v_frac(j+1,k)...
                          - 2*v_frac(j,k) + v_frac(j-1,k));            

            elseif (k==n-1)
                ap_e = 1+(dt/(Re*h^2));
                ae = 0;
                aw_e = dt/(2*Re*h^2);      
                f1_e = v_frac(j,k) +(dt/(2*Re*h^2)) * v(j,k+1)...
                       + (1/2) * (dt/(2*h)) * ((u_old(j,k+1)...
                       + u_old(j,k))*(v_old(j,k+1)+v_old(j,k))...
                       - (u_old(j,k-1)+u_old(j,k-1))*(v_old(j,k-1)+v_old(j,k))/2)...
                       - (3/2) * (dt/(2*h)) * ((u(j,k+1)+u(j,k))*(v(j,k+1)+v(j,k))...
                       - (u(j,k-1)+u(j,k-1))*(v(j,k-1)+v(j,k))/2)...
                       + (dt/(2*Re*h^2)) * (v_frac(j+1,k) ...
                       - 2*v_frac(j,k) + v_frac(j-1,k));

            elseif (k~=2 && k~=n-1)
                ap = 1+(dt/(Re*h^2));
                ae_I(1,k-2) = dt/(2*Re*h^2);         
                aw_I(1,k-2) = dt/(2*Re*h^2);
                f1_I(1,k-2) = v_frac(j,k)...
                              + (1/2) * (dt/(2*h)) * ((u_old(j,k+1)...
                              + u_old(j,k))*(v_old(j,k+1)+v_old(j,k))...
                              - (u_old(j,k-1)+u_old(j,k-1))*(v_old(j,k-1)+v_old(j,k))/2)...
                              - (3/2) * (dt/(2*h)) * ((u(j,k+1)+u(j,k))*(v(j,k+1)+v(j,k))...
                              - (u(j,k-1)+u(j,k-1))*(v(j,k-1)+v(j,k))/2)...
                              + (dt/(2*Re*h^2)) * (v_frac(j+1,k)...
                              - 2*v_frac(j,k) + v_frac(j-1,k));                        
            end
        end

        a(1,:) = ap;
        b = -[aw_I(1,:) aw_e];
        c = -[ae_w ae_I(1,:)];
        f = [f1_w f1_I(1,:) f1_e];
        TDMA(n-2,a,b,c,f)
        v_star1(j-1,:) = x;

    end
    
    for j = 2:n-1
            for k = 2:n-1
                v_star(j,k) = v_star1(j-1,k-1);
            end
    end

    %% Pressure Poisson Equation Solver
    ppe(n,h,dt,u_star,v_star)

    %% Velocity Correction
    for j = 2:n-1
       for k = 2:n-1
           u_new1(j-1,k-1) = u_star(j,k) - (dt/(2*h))*(p(j+1,k) - p(j-1,k));
           v_new1(j-1,k-1) = v_star(j,k) - (dt/(2*h))*(p(j,k+1) - p(j,k-1));
       end
    end   

    u_old = u;
    v_old = v;
    
    for j = 2:n-1
        for k = 2:n-1
            u_new(j,k) = u_new1(j-1,k-1);
            v_new(j,k) = v_new1(j-1,k-1);
        end
    end

    %% Check Residue From Continuity
    for j = 2:n-1
        for k = 2:n-1
            err(j-1,k-1) = (1/(2*h))*((u_new(j+1,k) - u_new(j-1,k))...
                                        + (v_new(j,k+1) - v_new(j,k-1)));
        end
    end

    u = u_new;
    v = v_new;

end
% Centerline Velocities

U = u_new(round((n-1)/2),:);
V = v_new(round((n-1)/2),:);
L = lx(1,:);

%% Plots
% box on
% axis square
% plot(lx,T)
% figure(1)
% z1 = quiver(lx,ly,u_new',v_new');
% streamslice(lx,ly,u_new',v_new');
% figure(2)
% z = contourf(lx,ly,p');
toc
