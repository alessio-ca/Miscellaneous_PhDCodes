%The second order integrator of the Generalized Langevin Equation (GLE) is
%based on Baczweski and Bond (J. Chem. Phys. 129, 044107 (2013))


close all; clc;

tic
%Parameter declaration (LJ style)

sigma = 1e-6; %Distance scale
m = 1e-15; %Mass scale
kB=1.38e-23; %Boltzmann constant
epsilon = kB*300; %Energy scale
tau = sqrt(sigma.^2*m/epsilon); %Time scale

R = 0.5;
rho = 2650*sigma.^3/m; %Density of silica in RU
eta = 0.001*sigma.^3/(epsilon*tau); %Viscosity of water in RU
T = 300*kB/epsilon; %Temperature
D = T/(3*pi*eta); %Diffusion coefficient of the bead in RU
k=[0 0 0]; %Trap elastic constant
r_eq = [0 0 0]; %Trap Equilibrium position
N = 5e+5; %Number of steps
run = 100; %Number of trajectories
r0 = [0 0 0]; %Initial position

%%
% Maxwell's model coefficients
tau_rel = 1e-3;
gamma = T/D;
Omega = sqrt(gamma/tau_rel - 1/(4*tau_rel.^2));

disp(['Omega (RU - SIU) = ', num2str(Omega), ' - ', num2str(Omega/tau),' rad/s'])
disp(['Tau_p (RU - SIU) = ', num2str(1/gamma), ' - ', num2str(1e6*tau/gamma),' us'])
disp(['Tau_M (RU - SIU) = ', num2str(1/(tau_rel.^(-1) + gamma)),' - ',num2str(1e6*tau/(tau_rel.^(-1) + gamma)),' us'])

% Inizialization
dt = 1e-5; %Time step
theta = exp(-dt/tau_rel);
alpha = (1-theta)/sqrt(dt);
r = zeros(N,3*run);
v = zeros(N,3*run);
S = zeros(N,3*run);
W = randn(N,3*run); % noise

r(1,:) = repmat(r0,1,run);
v(1,:) = zeros(1,3*run);
S(1,:) = zeros(1,3*run);


%%Equilibrate
for n = 2:1:1000
    vhalf = v(1,:) + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + S(1,:));
    r(1,:) = r(1,:) + dt*vhalf;
    S(1,:) = theta*S(1,:) - (1 - theta)*gamma*vhalf + alpha*sqrt(2*T*gamma)*W(n-1,:);
    Favg = mean(S(1,:)-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(1,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + S(1,:) - Favg);
end
W = randn(N,3*run); % noise
%%Sample
for n = 2:1:N
    vhalf = v(n-1,:) + dt/(2)*(-repmat(k,1,run).*(r(n-1,:)-repmat(r_eq,1,run)) + S(n-1,:));
    r(n,:) = r(n-1,:) + dt*vhalf;
    S(n,:) = theta*S(n-1,:) - (1 - theta)*gamma*vhalf + alpha*sqrt(2*T*gamma)*W(n-1,:);
    Favg = mean(S(n,:)-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(n,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run)) + S(n,:) - Favg);
    if mod(n,1e4)==0
        disp(n)
    end
end
toc
tic
[vacf,t_vacf,vacf_err]=acf_routine(dt,v,'TauMax',1e3*dt);
[msd,t_msd,msd_err]=msd_routine(dt,r,'TauMax',1e4*dt);
toc
%% 
%Plotting (in real units)
%Conversions
gammaREAL = (3*pi*eta*(epsilon*tau)/sigma.^2);
massREAL = rho*(4/3)*pi*R.^3*m;
msdREAL = sigma.^2 * msd;
msd_errREAL = sigma.^2 * msd_err;
t_vacfREAL = tau*t_vacf;
t_msdREAL = tau*t_msd;
kBT = T*epsilon;
D_REAL = (sigma^2/tau)*D; 

figure(1)
errorbar(t_vacfREAL(1:25:end),vacf(1:25:end),vacf_err(1:25:end),'o')
hold on
plot(t_vacfREAL,exp(-t_vacf/(2*tau_rel)).*(cos(Omega*t_vacf) + 1/(2*tau_rel*Omega)*sin(Omega*t_vacf)))
plot(tau*t_vacf,exp(-gammaREAL*t_vacfREAL/massREAL))
hold off
legend('Simulation','Theory','Plain Langevin')
xlabel('Time (s)')
ylabel('VACF (a.u.)')

figure(2)
index=unique(floor(logspace(0,log10(1e4))));
errorbar(t_msdREAL(index),1e+12*msdREAL(index),1e+12*msd_errREAL(index),'o')
hold on
plot(t_msdREAL,(1e+12)*(2*kBT/gammaREAL)*(t_msdREAL - (massREAL/gammaREAL)*(1 - exp(-gammaREAL*t_msdREAL/massREAL))));
plot(t_msdREAL,(1e+12)*2*D_REAL*(t_msdREAL));
hold off
set(gca,'XScale','log','YScale','log');
legend('Simulation','Plain Langevin','Plain Brownian')
xlabel('Time (s)')
ylabel('MSD (um^2)')
