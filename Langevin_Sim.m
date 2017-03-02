
close all;
clear all;

q=1; % 0 for Langevin dynamics, 1 for Brownian dynamics

T=1; % room temperature
N=5e4; % number of time steps
t=1; % total time
dt=t/N; % time step
kb=1; % Boltzmann constant
m=1; % particle mass
D=2; % diffusion coefficient
l=10*sqrt(2*D*t); % box length

X=zeros(N,2); % column vector for all particle positions
V=zeros(2,2); % row vector for initial (1,:) and final (2,:) velocity
V(1,:)=sqrt(kb*T/m).*randn(1,2); % initialisation of thermal velocities
X(1,:)=l/2; % initial position at centre of box




    for i=2:N % loop over all time steps           

        % integration of equations of motion 

        noise=sqrt(2*(kb.*T).^2./D/dt).*randn(1,2); % thermal noise force

        if q==1 % Brownian dynamics
            
            V(2,:)=D./(kb*T).*noise; % new velocity
            
        else % Langevin dynamics
            
            V(2,:)=V(1,:)+(-kb.*T./D.*V(1,:)+noise)/m*dt; % new velocity
            
        end

        X(i,:)=X(i-1,:)+V(1,:)*dt; % new position


     %%   

        % boundary conditions for wall reflections

        % for motion in x
        if (X(i,1)<0) % specular reflection

            X(i,1)=-X(i,1); % position after reflection
            V(1,1)=-V(2,1); % reverse velocity

        elseif (X(i,1)>l) % specular reflection

                X(i,1)=2*l-X(i,1); % position after reflection
                V(1,1)=-V(2,1); % reverse velocity

        else

            V(1,1)=V(2,1); % keep final velocity for next time step

        end
        
        % for motion in y
        if (X(i,2)<0) % specular reflection

            X(i,2)=-X(i,2); % postion after reflection
            V(1,2)=-V(2,2); % reverse velocity

        elseif (X(i,2)>l) % specular reflection

                X(i,2)=2*l-X(i,2); % postion after reflection
                V(1,2)=-V(2,2); % reverse velocity

        else

            V(1,2)=V(2,2); % keep final velocity for next time step

        end

%         plot(X(i,1),X(i,2),'o','MarkerSize',6,'LineWidth',2,'Color','r');
%         xlim([0,l])
%         ylim([0,l])
%         drawnow

    end
% tic
% [MSD,Dxy] = Langevin_MSD(X,dt,D);
% toc
tic
msd = msd_routine_alt(dt,X,'taumax',dt*1e4);
toc