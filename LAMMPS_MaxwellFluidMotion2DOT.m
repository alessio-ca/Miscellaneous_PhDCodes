%The second order integrator of the Generalized Langevin Equation (GLE) is
%based on Baczweski and Bond (J. Chem. Phys. 129, 044107 (2013))
tic

close all; clc;

%Parameter declaration

dt = 1e-8; %Time step
R = 0.5e-6; %Radius of bead
rho = 2650; %Density of the bead
eta = 0.001; %Viscosity (fictituous)
tau_rel = 5e-7; %Relaxation time
T = 300; %Temperature
k = [1e-6 1e-6 1e-6]; %Trap elastic constant
k=[0 0 0];
r_eq = [0 0 0]; %Trap Equilibrium position
N = 1e+4; %Number of steps
run = 100; %Number of trajectories
r0 = [0 0 0]; %Initial position

% Pre-calculation coefficients
m = (4./3.)*pi*R.^3*rho; 
kB=1.38e-23; %Boltzmann constant
kBT=kB*T;
D = kBT/(6*pi*eta*R);
gamma = kBT/D;
theta = exp(-dt/tau_rel);
alpha = (1-theta)/sqrt(dt);

% Inizialization

r = zeros(N,3*run);
v = zeros(N,3*run);
S = zeros(N,3*run);
W = randn(N,3*run); % noise

r(1,:) = repmat(r0,1,run);
v(1,:) = zeros(1,3*run);
S(1,:) = zeros(1,3*run);
Omega = sqrt(gamma/(m*tau_rel) - 1/(4*tau_rel.^2));
disp(['Omega = ', num2str(Omega)])
disp(['Tau_p = ', num2str(m/gamma)])
disp(['Tau_M = ', num2str(1/(tau_rel.^(-1) + gamma/m))])

%%Equilibrate%%
for n = 2:1:100
    vhalf = v(1,:) + dt/(2*m)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + S(1,:));
    r(1,:) = r(1,:) + dt*vhalf;
    S(1,:) = theta*S(1,:) - (1 - theta)*gamma*vhalf + alpha*sqrt(2*kBT*gamma)*W(n-1,:);
    Favg = mean(S(1,:)-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(1,:) = vhalf + dt/(2*m)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + S(1,:) - Favg);
end
W = randn(N,3*run); % re-instate noise

%%Sample%%
tic
%OnTheFly_CoarseGrainACF(v(1,:),dt,20,0);
for n = 2:1:N
    vhalf = v(n-1,:) + dt/(2*m)*(-repmat(k,1,run).*(r(n-1,:)-repmat(r_eq,1,run)) + S(n-1,:));
    r(n,:) = r(n-1,:) + dt*vhalf;
    S(n,:) = theta*S(n-1,:) - (1 - theta)*gamma*vhalf + alpha*sqrt(2*kBT*gamma)*W(n-1,:);
    Favg = mean(S(n,:)-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(n,:) = vhalf + dt/(2*m)*(-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run)) + S(n,:) - Favg);
    %OnTheFly_CoarseGrainACF(v(n,:),dt,20,1);
end
toc
%[vacfCG,t_vacfCG,vacfCG_err]=OnTheFly_CoarseGrainACF(v(1,:),dt,20,2,'TauMax',1e3*dt);
tic
[vacf,t_vacf,vacf_err]=acf_routine(dt,v,'TauMax',1e3*dt);
[msd,t_msd,msd_err]=msd_routine(dt,r);

toc
figure(1)
%errorbar(t_vacfCG,vacfCG,vacfCG_err/sqrt(run),'o')
errorbar(t_vacf,vacf,vacf_err/sqrt(run),'o')
hold on
plot(t_vacf,exp(-t_vacf/(2*tau_rel)).*(cos(Omega*t_vacf) + 1/(2*tau_rel*Omega)*sin(Omega*t_vacf)))
plot(t_vacf,exp(-gamma*t_vacf/m))
hold off
xlabel('Time')
ylabel('VACF')

figure(2)
errorbar(t_msd,1e+12*msd,1e+12*msd_err/sqrt(run),'o')
hold on
plot(t_msd,1e+12*2*D*t_msd)
hold off
set(gca,'XScale','log','YScale','log')
xlabel('Time')
ylabel('MSD')