%The second order integrator of the Generalized Langevin Equation (GLE) is
%based on Baczweski and Bond (J. Chem. Phys. 129, 044107 (2013))


close all; clc; clear all;

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
omega = 1.4;
k = (sigma.^2/epsilon)*[1e-6 1e-6]; %Trap elastic constant
k = [omega.^2 omega.^2];
r_eq = [0 0]; %Trap Equilibrium position
N = 1e+5; %Number of steps
run = 100; %Number of trajectories
r0 = [0 0]; %Initial position
dt = 1e-2; %Time step


%Power law kernel
lambda = 0.5;
gamma_l = 1;
K = @(x) gamma_l/gamma(1-lambda)*x.^(-0.5);

%% Power law kernel, Mode Study
tau_min = dt/10;
tau_max = 100;
t = linspace(tau_min,tau_max,ceil(1e4/dt));
figure(1)
plot(log10(t),log10(K(t)),'o','MarkerSize',1)
hold on
y = K(t)';
t = t';
modes = [2 4 6 8];
for i=1:length(modes)
    tau_k = logspace(log10(tau_min),log10(tau_max),modes(i));
    X = exp(-t./tau_k);
    a = X\y;
    plot(log10(t),log10(X*a))
end
hold off
legend('Exact','2 Modes','4 Modes','6 Modes','8 Modes')
xlabel('Time (Log)')
ylabel('K(t) (Log)')
drawnow

%% Power Law, 2 Modes

% Pre-calculation coefficients
nmodes=2;
dt = 1e-2; %Time step
tau_rel = logspace(log10(tau_min),log10(tau_max),nmodes);
X = exp(-t./tau_rel);
gamma = (X\y).*tau_rel';
theta = exp(-dt./tau_rel)';
alpha = (1-theta)./sqrt(dt)';

% Inizialization

r = zeros(N,2*run);
v = zeros(N,2*run);
S = zeros(nmodes*N,2*run);
W = randn(nmodes*N,2*run); % noise


r(1,:) = repmat(r0,1,run);
v(1,:) = zeros(1,2*run);
S(1:nmodes,:) = zeros(nmodes,2*run);
exampletitle('2 Modes Prony Series')

%%Equilibrate
for n = 2:1:100
    vhalf = v(1,:) + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + sum(S(1:nmodes,:)));
    r(1,:) = r(1,:) + dt*vhalf;
    S(1:nmodes,:) = theta.*S(1:nmodes,:) - (1 - theta).*gamma.*repmat(vhalf,[nmodes,1]) + alpha.*sqrt(2*T*gamma).*W(nmodes*(n-2)+1:nmodes*(n-1),:);
    Favg = mean(sum(S(1:nmodes,:))-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(1,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + sum(S(1:2,:)) - Favg);
end
W = randn(nmodes*N,2*run); % noise
%%Sample
tic
OnTheFly_CoarseGrainACF(v(1,:),dt,20,0,'P',64,'M',16);
for n = 2:1:N
    vhalf = v(n-1,:) + dt/(2)*(-repmat(k,1,run).*(r(n-1,:)-repmat(r_eq,1,run)) + sum(S(nmodes*(n-2)+1:nmodes*(n-1),:)));
    r(n,:) = r(n-1,:) + dt*vhalf;
    S(nmodes*(n-1)+1:nmodes*n,:) = theta.*S(nmodes*(n-2)+1:nmodes*(n-1),:) - (1 - theta).*gamma.*vhalf + alpha.*sqrt(2*T*gamma).*W(nmodes*(n-2)+1:nmodes*(n-1),:);
    Favg = mean(sum(S(nmodes*(n-1)+1:nmodes*n,:))-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(n,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run)) + sum(S(nmodes*(n-1)+1:nmodes*n,:)) - Favg);
    OnTheFly_CoarseGrainACF(v(n,:),dt,20,1);
end
toc
[vacfCG,t_vacfCG,vacfCG_err]=OnTheFly_CoarseGrainACF(v(1,:),dt,20,2,'TauMax',10);
figure(2)
errorbar(t_vacfCG,vacfCG,vacfCG_err/sqrt(run),'-o') %Errors are reported as standard errors
hold on
drawnow
disp(' ')
%% Power Law, 4 Modes

% Pre-calculation coefficients
nmodes=4;
dt = 1e-2; %Time step
tau_rel = logspace(log10(tau_min),log10(tau_max),nmodes);
X = exp(-t./tau_rel);
gamma = (X\y).*tau_rel';
theta = exp(-dt./tau_rel)';
alpha = (1-theta)./sqrt(dt)';

% Inizialization

r = zeros(N,2*run);
v = zeros(N,2*run);
S = zeros(nmodes*N,2*run);
W = randn(nmodes*N,2*run); % noise


r(1,:) = repmat(r0,1,run);
v(1,:) = zeros(1,2*run);
S(1:nmodes,:) = zeros(nmodes,2*run);
exampletitle('4 Modes Prony Series')

%%Equilibrate
for n = 2:1:100
    vhalf = v(1,:) + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + sum(S(1:nmodes,:)));
    r(1,:) = r(1,:) + dt*vhalf;
    S(1:nmodes,:) = theta.*S(1:nmodes,:) - (1 - theta).*gamma.*repmat(vhalf,[nmodes,1]) + alpha.*sqrt(2*T*gamma).*W(nmodes*(n-2)+1:nmodes*(n-1),:);
    Favg = mean(sum(S(1:nmodes,:))-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(1,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + sum(S(1:2,:)) - Favg);
end
W = randn(nmodes*N,2*run); % noise
%%Sample
tic
OnTheFly_CoarseGrainACF(v(1,:),dt,20,0,'P',64,'M',16);
for n = 2:1:N
    vhalf = v(n-1,:) + dt/(2)*(-repmat(k,1,run).*(r(n-1,:)-repmat(r_eq,1,run)) + sum(S(nmodes*(n-2)+1:nmodes*(n-1),:)));
    r(n,:) = r(n-1,:) + dt*vhalf;
    S(nmodes*(n-1)+1:nmodes*n,:) = theta.*S(nmodes*(n-2)+1:nmodes*(n-1),:) - (1 - theta).*gamma.*vhalf + alpha.*sqrt(2*T*gamma).*W(nmodes*(n-2)+1:nmodes*(n-1),:);
    Favg = mean(sum(S(nmodes*(n-1)+1:nmodes*n,:))-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(n,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run)) + sum(S(nmodes*(n-1)+1:nmodes*n,:)) - Favg);
    OnTheFly_CoarseGrainACF(v(n,:),dt,20,1);
end
toc
[vacfCG,t_vacfCG,vacfCG_err]=OnTheFly_CoarseGrainACF(v(1,:),dt,20,2,'TauMax',10);

errorbar(t_vacfCG,vacfCG,vacfCG_err/sqrt(run),'-o') %Errors are reported as standard errors
drawnow
disp(' ')
%% Power Law, 6 Modes

% Pre-calculation coefficients
nmodes=6;
dt = 1e-2; %Time step
tau_rel = logspace(log10(tau_min),log10(tau_max),nmodes);
X = exp(-t./tau_rel);
gamma = (X\y).*tau_rel';
theta = exp(-dt./tau_rel)';
alpha = (1-theta)./sqrt(dt)';

% Inizialization

r = zeros(N,2*run);
v = zeros(N,2*run);
S = zeros(nmodes*N,2*run);
W = randn(nmodes*N,2*run); % noise


r(1,:) = repmat(r0,1,run);
v(1,:) = zeros(1,2*run);
S(1:nmodes,:) = zeros(nmodes,2*run);
exampletitle('6 Modes Prony Series')

%%Equilibrate
for n = 2:1:100
    vhalf = v(1,:) + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + sum(S(1:nmodes,:)));
    r(1,:) = r(1,:) + dt*vhalf;
    S(1:nmodes,:) = theta.*S(1:nmodes,:) - (1 - theta).*gamma.*repmat(vhalf,[nmodes,1]) + alpha.*sqrt(2*T*gamma).*W(nmodes*(n-2)+1:nmodes*(n-1),:);
    Favg = mean(sum(S(1:nmodes,:))-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(1,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + sum(S(1:2,:)) - Favg);
end
W = randn(nmodes*N,2*run); % noise
%%Sample
tic
OnTheFly_CoarseGrainACF(v(1,:),dt,20,0,'P',64,'M',16);
for n = 2:1:N
    vhalf = v(n-1,:) + dt/(2)*(-repmat(k,1,run).*(r(n-1,:)-repmat(r_eq,1,run)) + sum(S(nmodes*(n-2)+1:nmodes*(n-1),:)));
    r(n,:) = r(n-1,:) + dt*vhalf;
    S(nmodes*(n-1)+1:nmodes*n,:) = theta.*S(nmodes*(n-2)+1:nmodes*(n-1),:) - (1 - theta).*gamma.*vhalf + alpha.*sqrt(2*T*gamma).*W(nmodes*(n-2)+1:nmodes*(n-1),:);
    Favg = mean(sum(S(nmodes*(n-1)+1:nmodes*n,:))-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(n,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run)) + sum(S(nmodes*(n-1)+1:nmodes*n,:)) - Favg);
    OnTheFly_CoarseGrainACF(v(n,:),dt,20,1);
end
toc
[vacfCG,t_vacfCG,vacfCG_err]=OnTheFly_CoarseGrainACF(v(1,:),dt,20,2,'TauMax',10);

errorbar(t_vacfCG,vacfCG,vacfCG_err/sqrt(run),'-o') %Errors are reported as standard errors
drawnow
disp(' ')
%% Power Law, 8 Modes

% Pre-calculation coefficients
nmodes=8;
dt = 1e-2; %Time step
tau_rel = logspace(log10(tau_min),log10(tau_max),nmodes);
X = exp(-t./tau_rel);
gamma = (X\y).*tau_rel';
theta = exp(-dt./tau_rel)';
alpha = (1-theta)./sqrt(dt)';

% Inizialization

r = zeros(N,2*run);
v = zeros(N,2*run);
S = zeros(nmodes*N,2*run);
W = randn(nmodes*N,2*run); % noise


r(1,:) = repmat(r0,1,run);
v(1,:) = zeros(1,2*run);
S(1:nmodes,:) = zeros(nmodes,2*run);
exampletitle('8 Modes Prony Series')

%%Equilibrate
for n = 2:1:100
    vhalf = v(1,:) + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + sum(S(1:nmodes,:)));
    r(1,:) = r(1,:) + dt*vhalf;
    S(1:nmodes,:) = theta.*S(1:nmodes,:) - (1 - theta).*gamma.*repmat(vhalf,[nmodes,1]) + alpha.*sqrt(2*T*gamma).*W(nmodes*(n-2)+1:nmodes*(n-1),:);
    Favg = mean(sum(S(1:nmodes,:))-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(1,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + sum(S(1:2,:)) - Favg);
end
W = randn(nmodes*N,2*run); % noise
%%Sample
tic
OnTheFly_CoarseGrainACF(v(1,:),dt,20,0,'P',64,'M',16);
for n = 2:1:N
    vhalf = v(n-1,:) + dt/(2)*(-repmat(k,1,run).*(r(n-1,:)-repmat(r_eq,1,run)) + sum(S(nmodes*(n-2)+1:nmodes*(n-1),:)));
    r(n,:) = r(n-1,:) + dt*vhalf;
    S(nmodes*(n-1)+1:nmodes*n,:) = theta.*S(nmodes*(n-2)+1:nmodes*(n-1),:) - (1 - theta).*gamma.*vhalf + alpha.*sqrt(2*T*gamma).*W(nmodes*(n-2)+1:nmodes*(n-1),:);
    Favg = mean(sum(S(nmodes*(n-1)+1:nmodes*n,:))-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(n,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run)) + sum(S(nmodes*(n-1)+1:nmodes*n,:)) - Favg);
    OnTheFly_CoarseGrainACF(v(n,:),dt,20,1);
end
toc
[vacfCG,t_vacfCG,vacfCG_err]=OnTheFly_CoarseGrainACF(v(1,:),dt,20,2,'TauMax',10);

errorbar(t_vacfCG,vacfCG,vacfCG_err/sqrt(run),'-o') %Errors are reported as standard errors
plot(t_vacfCG,GLE_VACF(t_vacfCG,omega,lambda,gamma_l))
hold off
legend('2 Modes','4 Modes','6 Modes','8 Modes','Analytical')
title('GLE with Power law memory kernel')
xlabel('Time')
ylabel('VACF')