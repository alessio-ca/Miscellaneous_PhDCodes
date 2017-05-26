%The first order integrator of the Langevin Equation is based on
%Volpe (Am. J. Phys., 81, 3, 2013)

close all; clc;

%Parameter declaration

dt = 1e-8; %Time step
R = 1e-6; %Radius of bead
rho = 2650; %Density of the bead
eta = 0.001; %Viscosity
T = 300; %Temperature
k = [1e-6 1e-6]; %Trap elastic constant
r_eq = [0 0]; %Trap Equilibrium position
N = 1e+8; %Number of steps

r0 = [0 0]; %Initial position

% Pre-calculation coefficients
m = (4./3.)*pi*R.^3*rho; 
kB=1.38e-23; %Boltzmann constant
kBT=kB*T;
D = kBT/(6*pi*eta*R);
gamma = kBT/D;
sigma = sqrt(2*kBT*gamma)/m;
C1 = 1  + dt*(gamma/m) + dt.^2*(k/m);
C2 = (2  + dt*(gamma/m))./C1;
h = randn(N,2);
hs = sqrt(2*D*dt)*h; % scaled noise

% Inizialization
r = zeros(N,2);
v = zeros(N,2);
rOV = zeros(N,2);

r(1,:) = r0;
rOV(1,:) = r0;
v(1,:) = [0 0];
t = (0:1:N-1)'*dt;


% Main loop
% It simulates a Brownian Particle in a solvent according to the
% non-overdamped Langevin equation. 
for n = 3:1:N   
    r(n,:) = C2.*r(n-1,:) - (1./C1).*r(n-2,:) + sqrt(2*kBT*gamma)./(m.*C1).*dt.^(1.5).*h(n,:);  
    if mod(n,10000)==0
        disp(n)
    end
end


C1 = 1 - k/gamma*dt;
C2 = k.*r_eq/gamma*dt;
% Main loop
% It simulates a Brownian Particle in a solvent according to the overdamped
% Langevin equation. 
for n = 2:1:N
    rOV(n,:) = C1.*rOV(n-1,:) + C2 + hs(n,:);
    if mod(n,10000)==0
        disp(n)
    end
end

