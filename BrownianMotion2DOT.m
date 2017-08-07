close all; clc;

%Parameter declaration

dt = 1e-4; %Time step
R = 0.75e-6; %Radius of bead
eta = 0.001; %Viscosity
T = 300; %Temperature
k = [17.38e-6 15.7e-6]; %Trap elastic constant
r_eq = [0 0]; %Trap Equilibrium position
N = 1e+6; %Number of steps
run = 10; %Number of trajectories


r0 = [0 0]; %Initial position

% Pre-calculation coefficients
kB=1.38e-23; %Boltzmann constant
kBT=kB*T;
D = kBT/(6*pi*eta*R);
gamma = kBT/D;
C1 = 1 - k/gamma*dt;
C2 = k.*r_eq/gamma*dt;
h = randn(N,2*run); % noise
hs = sqrt(2*D*dt)*h; % scaled noise


% Inizialization
r = zeros(N,2*run);
r(1,:) = repmat(r0,1,run);
t = (0:1:N-1)'*dt;


% Main loop
% It simulates a Brownian Particle in a solvent according to the overdamped
% Langevin equation. 
for n = 2:1:N
    r(n,:) = repmat(C1,1,run).*r(n-1,:) + repmat(C2,1,run) + hs(n,:);
    if mod(n,10000)==0
        disp(n)
    end
end

%Shuffles x & y coordinates to separate them
r = [r(:,1:2:end),r(:,2:2:end)];

%%
hsup = 0.25*sqrt(2*D*dt)*randn(N,2*run); %added noise
reff = r + hsup;
[msd_x,~]=msd_routine(dt,reff(:,1:10));
[msd_y,tau]=msd_routine(dt,reff(:,11:20));
loglog(tau,msd_x,'bo')
hold on
loglog(tau,msd_y,'ro')
hold off
legend('MSD,x','MSD,y')
