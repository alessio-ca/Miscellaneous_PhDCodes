%The first order integrator of the Langevin Equation is based on
%Volpe (Am. J. Phys., 81, 3, 2013)

close all; clc;

%Parameter declaration

dt = 1e-6; %Time step
R = 0.115e-6; %Radius of bead
rho = 1.04*1e3; %Density of the bead (in kg/m^3)
eta = 0.0008872; %Fluid viscosity (in Pa.s)
T = 298; %Temperature
k = [0e-6 0e-6 0e-6]; %Trap elastic constant
r_eq = [0 0 0]; %Trap Equilibrium position
N = 1e+7; %Number of steps
run = 10; %Numer of different runs

r0 = [0 0 0]; %Initial position

% Pre-calculation coefficients
m = (4./3.)*pi*R.^3*rho; 
kB=1.38e-23; %Boltzmann constant
kBT=kB*T;
D = kBT/(6*pi*eta*R);
gamma = kBT/D;
sigma = sqrt(2*kBT*gamma)/m;
C1 = 1  + dt*(gamma/m) + dt.^2*(k/m);
C2 = (2  + dt*(gamma/m))./C1;
h = randn(N,3*run);
hs = sqrt(2*D*dt)*h; % scaled noise

% Inizialization
r = zeros(N,3*run);
v = zeros(N,3*run);
rOV = zeros(N,3*run);

r(1,:) = repmat(r0,1,run);
v(1,:) = zeros(1,3*run);
rOV(1,:) = repmat(r0,1,run);
t = (0:1:N-1)'*dt;


% Main loop
% It simulates a Brownian Particle in a solvent according to the
% non-overdamped Langevin equation.
for n = 3:1:N
    r(n,:) = repmat(C2,1,run).*r(n-1,:) - (1./repmat(C1,1,run)).*r(n-2,:) + sqrt(2*kBT*gamma)./(m.*repmat(C1,1,run)).*dt.^(1.5).*h(n,:);
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
    rOV(n,:) = repmat(C1,1,run).*rOV(n-1,:) + repmat(C2,1,run) + hs(n,:);
    if mod(n,10000)==0
        disp(n)
    end
end

