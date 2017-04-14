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
k=[0 0 0]; %Trap elastic constant
r_eq = [0 0 0]; %Trap Equilibrium position
N = 5e+4; %Number of steps
run = 1000; %Number of trajectories
r0 = [0 0 0]; %Initial position

%%
% Maxwell's model coefficients
tau_rel = 1;
gamma = 1;
Omega = sqrt(gamma/tau_rel - 1/(4*tau_rel.^2));
disp(['Omega = ', num2str(Omega)])
disp(['Tau_p = ', num2str(1/gamma)])
disp(['Tau_M = ', num2str(1/(tau_rel.^(-1) + gamma))])

% Inizialization
dt = 1e-2; %Time step
theta = exp(-dt/tau_rel);
alpha = (1-theta)/sqrt(dt);
r = zeros(N,3*run);
v = zeros(N,3*run);
S = zeros(N,3*run);
W = randn(N,3*run); % noise
r(1,:) = repmat(r0,1,run);
v(1,:) = zeros(1,3*run);
S(1,:) = zeros(1,3*run);


%Equilibrate
for n = 2:1:1000
    vhalf = v(1,:) + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + S(1,:));
    r(1,:) = r(1,:) + dt*vhalf;
    S(1,:) = theta*S(1,:) - (1 - theta)*gamma*vhalf + alpha*sqrt(2*T*gamma)*W(n-1,:);
    Favg = [mean(S(1,1:3:end)-repmat(k(1),1,run).*(r(1,1:3:end)-repmat(r_eq(1),1,run))),...
            mean(S(1,2:3:end)-repmat(k(2),1,run).*(r(1,2:3:end)-repmat(r_eq(2),1,run))),...
            mean(S(1,3:3:end)-repmat(k(3),1,run).*(r(1,3:3:end)-repmat(r_eq(3),1,run)))]; %Center of mass acceleration, to be substracted
    v(1,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(1,:)-repmat(r_eq,1,run)) + S(1,:) - repmat(Favg,1,run));
end
W = randn(N,3*run); % noise
%%Sample
OnTheFly_CoarseGrainACF(v(1,:),dt,20,0);
OnTheFly_CoarseGrainACF(v(1,:),dt,20,1);
for n = 2:1:N
    vhalf = v(n-1,:) + dt/(2)*(-repmat(k,1,run).*(r(n-1,:)-repmat(r_eq,1,run)) + S(n-1,:));
    r(n,:) = r(n-1,:) + dt*vhalf;
    S(n,:) = theta*S(n-1,:) - (1 - theta)*gamma*vhalf + alpha*sqrt(2*T*gamma)*W(n-1,:);
    Favg = [mean(S(n,1:3:end)-repmat(k(1),1,run).*(r(n,1:3:end)-repmat(r_eq(1),1,run))),...
            mean(S(n,2:3:end)-repmat(k(2),1,run).*(r(n,2:3:end)-repmat(r_eq(2),1,run))),...
            mean(S(n,3:3:end)-repmat(k(3),1,run).*(r(n,3:3:end)-repmat(r_eq(3),1,run)))]; %Center of mass acceleration, to be substracted
    v(n,:) = vhalf + dt/(2)*(-repmat(k,1,run).*(r(n,:)-repmat(r_eq,1,run)) + S(n,:) - repmat(Favg,1,run));
    OnTheFly_CoarseGrainACF(v(n,:),dt,20,1);
end
toc
[vacfCG,t_vacfCG,vacfCG_err]=OnTheFly_CoarseGrainACF(v(1,:),dt,20,2,'TauMax',10);
[vacf,t_vacf,vacf_err]=acf_routine(dt,v,'TauMax',10);

figure(1)
errorbar(t_vacfCG,vacfCG,vacfCG_err/sqrt(run),'o')
hold on
plot(t_vacf,exp(-t_vacf/(2*tau_rel)).*(cos(Omega*t_vacf) + 1/(2*tau_rel*Omega)*sin(Omega*t_vacf)))
plot(t_vacf,exp(-gamma*t_vacf))
hold off
xlabel('Time')
ylabel('VACF')
