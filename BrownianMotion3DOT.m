close all; clc;

%Parameter declaration

dt = 1e-3; %Time step
R = 1e-6; %Radius of bead
eta = 0.001; %Viscosity
T = 300; %Temperature
k = [1e-6 1e-6 1e-7]; %Trap elastic constant
r_eq = [0 0 0]; %Trap Equilibrium position
N = 1e+5; %Number of steps

r0 = [0 0 0]; %Initial position

% Pre-calculation coefficients
kB=1.38e-23; %Boltzmann constant
kBT=kB*T;
D = kBT/(6*pi*eta*R);
gamma = kBT/D;
C1 = 1 - k/gamma*dt;
C2 = k.*r_eq/gamma*dt;
h = randn(N,3); % noise
hs = sqrt(2*D*dt)*h; % scaled noise

% Inizialization
r = zeros(N,3);
r(1,:) = r0;
t = (0:1:N-1)'*dt;


% Simulation
% It simulates a 
for n = 2:1:N
    r(n,:) = C1.*r(n-1,:) + C2 + hs(n,:);
    if mod(n,10000)==0
        disp(n)
    end
end

