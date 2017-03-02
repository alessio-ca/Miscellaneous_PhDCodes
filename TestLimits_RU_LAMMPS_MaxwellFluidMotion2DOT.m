%The second order integrator of the Generalized Langevin Equation (GLE) is
%based on Baczweski and Bond (J. Chem. Phys. 129, 044107 (2013))


close all; clc;

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
k=[0 0]; %Trap elastic constant
r_eq = [0 0]; %Trap Equilibrium position
N = 1e+5; %Number of steps
run = 100; %Number of trajectories
r0 = [0 0]; %Initial position

%% Underdamped regime

% Pre-calculation coefficients
dt = 1e-2; %Time step
tau_rel = 1;
gamma = 1;
theta = exp(-dt/tau_rel);
alpha = (1-theta)/sqrt(dt);

% Inizialization

r = zeros(N,2*run);
v = zeros(N,2*run);
S = zeros(N,2*run);
W = randn(N,2*run); % noise

r(1,:) = repmat(r0,1,run);
v(1,:) = zeros(1,2*run);
S(1,:) = zeros(1,2*run);
Omega = sqrt(gamma/tau_rel - 1/(4*tau_rel.^2));
exampletitle('Underdamped regime')
disp(['Omega = ', num2str(Omega)])
disp(['Tau_p = ', num2str(1/gamma)])
disp(['Tau_M = ', num2str(1/(tau_rel.^(-1) + gamma))])

%%Equilibrate
for n = 2:1:100
    vhalf = v(1,:) + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + S(1,:));
    r(1,:) = r(1,:) + dt*vhalf;
    S(1,:) = theta*S(1,:) - (1 - theta)*gamma*vhalf + alpha*sqrt(2*T*gamma)*W(n-1,:);
    Favg = mean(S(1,:)-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(1,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + S(1,:) - Favg);
end
W = randn(N,2*run); % noise
%%Sample
tic
OnTheFly_CoarseGrainACF(v(1,:),dt,20,0);
for n = 2:1:N
    vhalf = v(n-1,:) + dt/(2)*(-repmat(k,1,run).*(r(n-1,:)-repmat(r_eq,1,run)) + S(n-1,:));
    r(n,:) = r(n-1,:) + dt*vhalf;
    S(n,:) = theta*S(n-1,:) - (1 - theta)*gamma*vhalf + alpha*sqrt(2*T*gamma)*W(n-1,:);
    Favg = mean(S(n,:)-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(n,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run)) + S(n,:) - Favg);
    OnTheFly_CoarseGrainACF(v(n,:),dt,20,1);
end
toc
[vacfCG,t_vacfCG,vacfCG_err]=OnTheFly_CoarseGrainACF(v(1,:),dt,20,2,'TauMax',10);
[vacf,t_vacf,vacf_err]=acf_routine(dt,v,'TauMax',10);
figure(1)
subplot(2,2,1)
errorbar(t_vacfCG,vacfCG,vacfCG_err/sqrt(run),'o')
hold on
plot(t_vacf,exp(-t_vacf/(2*tau_rel)).*(cos(Omega*t_vacf) + 1/(2*tau_rel*Omega)*sin(Omega*t_vacf)))
hold off
legend('Simulation','Analytical')
title('Underdamped regime')
xlabel('Time')
ylabel('VACF')
drawnow
disp(' ')
%% Critically damped regime

% Pre-calculation coefficients
dt = 1e-2; %Time step
tau_rel = .5;
gamma = .5;
theta = exp(-dt/tau_rel);
alpha = (1-theta)/sqrt(dt);

% Inizialization

r = zeros(N,2*run);
v = zeros(N,2*run);
S = zeros(N,2*run);
W = randn(N,2*run); % noise

r(1,:) = repmat(r0,1,run);
v(1,:) = zeros(1,2*run);
S(1,:) = zeros(1,2*run);
Omega = sqrt(gamma/tau_rel - 1/(4*tau_rel.^2));
exampletitle('Critically damped regime')
disp(['Omega = ', num2str(Omega)])
disp(['Tau_p = ', num2str(1/gamma)])
disp(['Tau_M = ', num2str(1/(tau_rel.^(-1) + gamma))])

%%Equilibrate
for n = 2:1:100
    vhalf = v(1,:) + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + S(1,:));
    r(1,:) = r(1,:) + dt*vhalf;
    S(1,:) = theta*S(1,:) - (1 - theta)*gamma*vhalf + alpha*sqrt(2*T*gamma)*W(n-1,:);
    Favg = mean(S(1,:)-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(1,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + S(1,:) - Favg);
end
W = randn(N,2*run); % noise
%%Sample
tic
OnTheFly_CoarseGrainACF(v(1,:),dt,20,0);
for n = 2:1:N
    vhalf = v(n-1,:) + dt/(2)*(-repmat(k,1,run).*(r(n-1,:)-repmat(r_eq,1,run)) + S(n-1,:));
    r(n,:) = r(n-1,:) + dt*vhalf;
    S(n,:) = theta*S(n-1,:) - (1 - theta)*gamma*vhalf + alpha*sqrt(2*T*gamma)*W(n-1,:);
    Favg = mean(S(n,:)-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(n,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run)) + S(n,:) - Favg);
    OnTheFly_CoarseGrainACF(v(n,:),dt,20,1);
end
toc
[vacfCG,t_vacfCG,vacfCG_err]=OnTheFly_CoarseGrainACF(v(1,:),dt,20,2,'TauMax',10);
figure(1)
subplot(2,2,2)
errorbar(t_vacfCG,vacfCG,vacfCG_err/sqrt(run),'o')
hold on
plot(t_vacf,exp(-t_vacf/(2*tau_rel)).*(1+t_vacf/(2*tau_rel)))
hold off
legend('Simulation','Analytical')
title('Critically damped regime')
xlabel('Time')
ylabel('VACF')
drawnow
disp(' ')
%% Overdamped regime

% Pre-calculation coefficients
dt = 1e-2; %Time step
tau_rel = .25;
gamma = .25;
theta = exp(-dt/tau_rel);
alpha = (1-theta)/sqrt(dt);

% Inizialization

r = zeros(N,2*run);
v = zeros(N,2*run);
S = zeros(N,2*run);
W = randn(N,2*run); % noise

r(1,:) = repmat(r0,1,run);
v(1,:) = zeros(1,2*run);
S(1,:) = zeros(1,2*run);
Omega = sqrt(gamma/tau_rel - 1/(4*tau_rel.^2));
exampletitle('Overdamped regime')
disp(['Omega = ', num2str(Omega)])
disp(['Tau_p = ', num2str(1/gamma)])
disp(['Tau_M = ', num2str(1/(tau_rel.^(-1) + gamma))])

%%Equilibrate
for n = 2:1:100
    vhalf = v(1,:) + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + S(1,:));
    r(1,:) = r(1,:) + dt*vhalf;
    S(1,:) = theta*S(1,:) - (1 - theta)*gamma*vhalf + alpha*sqrt(2*T*gamma)*W(n-1,:);
    Favg = mean(S(1,:)-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(1,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + S(1,:) - Favg);
end
W = randn(N,2*run); % noise
%%Sample
tic
OnTheFly_CoarseGrainACF(v(1,:),dt,20,0);
for n = 2:1:N
    vhalf = v(n-1,:) + dt/(2)*(-repmat(k,1,run).*(r(n-1,:)-repmat(r_eq,1,run)) + S(n-1,:));
    r(n,:) = r(n-1,:) + dt*vhalf;
    S(n,:) = theta*S(n-1,:) - (1 - theta)*gamma*vhalf + alpha*sqrt(2*T*gamma)*W(n-1,:);
    Favg = mean(S(n,:)-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(n,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run)) + S(n,:) - Favg);
    OnTheFly_CoarseGrainACF(v(n,:),dt,20,1);
end
toc
[vacfCG,t_vacfCG,vacfCG_err]=OnTheFly_CoarseGrainACF(v(1,:),dt,20,2,'TauMax',10);
figure(1)
subplot(2,2,3)
errorbar(t_vacfCG,vacfCG,vacfCG_err/sqrt(run),'o')
hold on
plot(t_vacf,exp(-t_vacf/(2*tau_rel)).*(cos(Omega*t_vacf) + 1/(2*tau_rel*Omega)*sin(Omega*t_vacf)))
hold off
legend('Simulation','Analytical')
title('Overdamped regime')
xlabel('Time')
ylabel('VACF')
drawnow
disp(' ')
%% Langevin regime

% Pre-calculation coefficients
dt = 1e-2; %Time step
tau_rel = 1e-8;
gamma = 1;
theta = exp(-dt/tau_rel);
alpha = (1-theta)/sqrt(dt);

% Inizialization

r = zeros(N,2*run);
v = zeros(N,2*run);
S = zeros(N,2*run);
W = randn(N,2*run); % noise

r(1,:) = repmat(r0,1,run);
v(1,:) = zeros(1,2*run);
S(1,:) = zeros(1,2*run);
t = (0:1:N-1)'*dt;
Omega = sqrt(gamma/tau_rel - 1/(4*tau_rel.^2));
exampletitle('Langevin regime')
disp(['Omega = ', num2str(Omega)])
disp(['Tau_p = ', num2str(1/gamma)])
disp(['Tau_M = ', num2str(1/(tau_rel.^(-1) + gamma))])

%%Equilibrate
for n = 2:1:100
    vhalf = v(1,:) + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + S(1,:));
    r(1,:) = r(1,:) + dt*vhalf;
    S(1,:) = theta*S(1,:) - (1 - theta)*gamma*vhalf + alpha*sqrt(2*T*gamma)*W(n-1,:);
    Favg = mean(S(1,:)-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(1,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + S(1,:) - Favg);
end
W = randn(N,2*run); % noise
%%Sample
tic
OnTheFly_CoarseGrainACF(v(1,:),dt,20,0);
for n = 2:1:N
    vhalf = v(n-1,:) + dt/(2)*(-repmat(k,1,run).*(r(n-1,:)-repmat(r_eq,1,run)) + S(n-1,:));
    r(n,:) = r(n-1,:) + dt*vhalf;
    S(n,:) = theta*S(n-1,:) - (1 - theta)*gamma*vhalf + alpha*sqrt(2*T*gamma)*W(n-1,:);
    Favg = mean(S(n,:)-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run))); %Center of mass acceleration, to be substracted
    v(n,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run)) + S(n,:) - Favg);
    OnTheFly_CoarseGrainACF(v(n,:),dt,20,1);
    if mod(n,1000000)==0
        disp(n)
    end
end
toc
[vacfCG,t_vacfCG,vacfCG_err]=OnTheFly_CoarseGrainACF(v(1,:),dt,20,2,'TauMax',10);
figure(1)
subplot(2,2,4)
errorbar(t_vacfCG,vacfCG,vacfCG_err/sqrt(run),'o')
hold on
plot(t_vacf,exp(-gamma*t_vacf))
hold off
legend('Simulation','Analytical')
title('Langevin regime')
xlabel('Time')
ylabel('VACF')
drawnow
