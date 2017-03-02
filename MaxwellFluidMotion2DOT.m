%The second order integrator of the Generalized Langevin Equation (GLE) is
%based on Gordon et al. (Mol. Phys., 106:11 (2008))

tic

close all; clc;

%Parameter declaration

dt = 1e-10; %Time step
R = 1e-6; %Radius of bead
rho = 2650; %Density of the bead
eta = 0.001; %Viscosity (fictituous)
tau = 1e-6; %Relaxation time
T = 300; %Temperature
k = [1e-6 1e-6]; %Trap elastic constant
k=[0 0];
r_eq = [0 0]; %Trap Equilibrium position
N = 1e+5; %Number of steps

r0 = [0 0]; %Initial position
v0 = [0 0]; %Initial velocity

% Pre-calculation coefficients
m = (4./3.)*pi*R.^3*rho; 
kB=1.38e-23; %Boltzmann constant
kBT=kB*T;
D = kBT/(6*pi*eta*R);
gamma = kBT/D;
K0 = gamma/tau;

%Useful functions
f1 = @(x) 1 - exp(2*x);
f2 = @(x) cosh(x) - 1;
f3 = @(x) sinh(x) - x;
f4 = @(x) exp(2*x) - 1 - 2*x - 4*(exp(x) - 1 - x);
f5 = @(x) cosh(x) - 1 - 0.5*x.^2;
f6 = @(x) 4*x*(exp(x) - 1 - x - 0.5*x.^2) - (exp(2*x) - 1 - (2*x) - 0.5*(2*x).^2 - (1./6.)*(2*x).^3);


%Covariance matrices
WVZ = (K0*kBT/m)*[f1(-dt/(2*tau))                   -2*f2(dt/(2*tau))   2*exp(-dt/(2*tau))*f3(dt/(2*tau));...
    -2*f2(dt/(2*tau))                 f4(dt/(2*tau))      -2*f5(dt/(2*tau))                ;...
    2*exp(-dt/(2*tau))*f3(dt/(2*tau)) -2*f5(dt/(2*tau))   f6(-dt/(2*tau))                     ];
WYX = [WVZ(1,1)                      -exp(-dt/(2*tau))*WVZ(1,2)  exp(dt/(2*tau))*WVZ(1,3);...
    -exp(-dt/(2*tau))*WVZ(1,2)   -(K0*kBT/m)*f4(-dt/(2*tau)) -WVZ(2,3)               ;...
    exp(dt/(2*tau))*WVZ(1,3)     -WVZ(2,3)                    -(K0*kBT/m)*f6(dt/(2*tau))];

[L1,p1] = chol(WVZ,'lower');
[L2,p2] = chol(WYX,'lower');
if any([p1 p2] > [0 0])
    disp('Warning: low argument in asymptotic expansions. Using Taylor expansion for the f* functions...')
    f1 = @(x) -2*x;
    f2 = @(x) 0.5*x.^2;
    f3 = @(x) (1./6.)*x.^3;
    f4 = @(x) (2./3.)*x.^3;
    f5 = @(x) (1./24.)*x.^4;
    f6 = @(x) (-1./10.)*x.^5;
    %Covariance matrices
    WVZ = (K0*kBT/m)*[f1(-dt/(2*tau))                   -2*f2(dt/(2*tau))   2*exp(-dt/(2*tau))*f3(dt/(2*tau));...
        -2*f2(dt/(2*tau))                 f4(dt/(2*tau))      -2*f5(dt/(2*tau))                ;...
        2*exp(-dt/(2*tau))*f3(dt/(2*tau)) -2*f5(dt/(2*tau))   f6(-dt/(2*tau))                     ];
    WYX = [WVZ(1,1)                      -exp(-dt/(2*tau))*WVZ(1,2)  exp(dt/(2*tau))*WVZ(1,3);...
        -exp(-dt/(2*tau))*WVZ(1,2)   -(K0*kBT/m)*f4(-dt/(2*tau)) -WVZ(2,3)               ;...
        exp(dt/(2*tau))*WVZ(1,3)     -WVZ(2,3)                    -(K0*kBT/m)*f6(dt/(2*tau))];
    L1 = chol(WVZ,'lower');
    L2 = chol(WYX,'lower');
end
   
%%
% Inizialization
r = zeros(N,2);
v = zeros(N,2);
W = zeros(2*N,2);
Y = zeros(N,2);
X = zeros(N,2);
V = zeros(N,2);
Z = zeros(N,2);


r(1,:) = r0;
D0 = -k.*(r(1,:) - r_eq);
F0 = -gamma*v0;
S0 = sqrt(K0*kBT/m)*randn(1,2);
v(1,:) = v0 - 0.5*dt*(D0 + F0 + S0);
A=(L1*randn(2,3)')';
W(1,:)=A(:,1);
V(1,:)=A(:,2);
Z(1,:)=A(:,3);
t = (0:1:N-1)'*dt;
%%

% Main loop
% It simulates a Brownian Particle in a solvent according to the
% non-overdamped Generalized Langevin Equation.
for n = 2:1:N
    A=(L2*randn(2,3)')';
    W(2*n-2,:)=A(:,1);
    Y(n-1,:)=A(:,2);
    X(n-1,:)=A(:,3);
    A=(L1*randn(2,3)')';
    W(2*n-1,:)=A(:,1);
    V(n,:)=A(:,2);
    Z(n,:)=A(:,3);
    v(n,:) = v(n-1,:) + dt.*(D0 + F0 + S0) + tau*(V(n-1,:)+Y(n-1,:));
    S1 = S0*exp(-dt./(2*tau)) +  W(2*n-2,:);
    r(n,:) = r(n-1,:) + dt*v(n,:) - dt.^3/(24.*tau)*S1 + tau.^2*(X(n-1,:) + Z(n,:));
    F0 = exp(-dt/tau)*F0 - K0*tau*(1-exp(-dt/tau))*v(n,:);
    S0 = S1*exp(-dt/(2*tau)) + W(2*n-1,:);
    D0 = -k.*(r(n,2) - r_eq);
    if mod(n,1000000)==0
        disp(n)
    end
end

[vacf,t_vacf]=acf_routine(dt,v,'TauMax',dt*length(t)/4);
semilogx(t_vacf,vacf)
hold on
semilogx(t_vacf,exp(-gamma*t_vacf/m))
hold off

