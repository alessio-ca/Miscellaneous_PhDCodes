tic
%The second order integrator of the Langevin Equation is based on
%Vanden-Eijden and Ciccotti (Chem. Phys. Lett. 429, 310316 (2006))

close all; clc;

%Parameter declaration
N = 2e6; %Number of steps
dt = 5e-2; %Time step
dt0_5 = sqrt(dt);
dt1_5 = dt.^(1.5);
dt2 = dt.^2;
R = 0.5; %Radius of bead
kBT = 1; %Temperature
m = 1; %Particle mass
D = 1; %Diffusion coefficient

k = [0.01 0.01]; %Trap elastic constant

r_eq = [0 0]; %Trap Equilibrium position
r0 = [0 0]; %Initial position


gamma = kBT/D;
gamma = gamma/m; %scaled friction coefficient
k = k/m; %scaled elastic constant
sigma = sqrt(2*kBT*gamma/m);
C1 = k.*r_eq*(0.5*dt.^2);
eps = randn(N,2); % noise, epsilon
th = randn(N,2); % noise, theta

% Inizialization
r = zeros(N,2);
v = zeros(N,2);
r(1,:) = r0;
v(1,:) = [0 0];
t = (0:1:N-1)'*dt;

% Main loop
% It simulates a Brownian Particle in a solvent according to the
% non-overdamped Langevin equation. 
for n = 2:1:N
    C = (0.5*dt2)*(-k.*r(n-1,:)- gamma*v(n-1,:)) + C1 + sigma*dt1_5*(0.5*eps(n,:) + (0.5/sqrt(3.))*th(n,:));
    r(n,:) = r(n-1,:) + dt*v(n-1,:) + C;
    v(n,:) = v(n-1,:) + 0.5*dt*(-k.*r(n,:) - k.*r(n-1,:)) + 2*C1/dt - dt*gamma*v(n-1,:) + sigma*dt0_5*eps(n,:) - gamma*C;
    if mod(n,100000)==0
        disp(n)
    end
end
    
%If trap is symmetric, both coordinates obey the same Langevin Equation
%MSD and VACF can be calculated by averaging over X and Y
[msd,t_msd]=msd_routine(dt,r,'TauMax',dt*length(t)/10);
[vacf,t_vacf]=acf_routine(dt,v,'TauMax',dt*length(t)/100);

%If trap is not symmetric, coodrinates do not obey the same Langevin Equation
%MSD and VACF have to be calculated separately for X and Y
%%
figure(1)
loglog(t_msd,msd)
hold on
loglog(t_msd,2*D*t_msd)
hold off
title('Mean Square Displacement')
xlabel('Time')
ylabel('MSD')

figure(2)
semilogx(t_vacf,vacf)
hold on
semilogx(t_vacf,exp(-gamma*t_vacf))
hold off
title('Velocity Autocorrelation Function')
xlabel('Time')
ylabel('VACF')

toc