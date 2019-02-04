%% CONTIN_CRIE Routine Tests

clear all
close all

T = 300;
eta = 1e-3;
n = 1.33;
R = 115e-9;
beta = 1;
Ng = 161;
cut = 0;
Nnoise = 6;
tail = [0,0];

%% Constants
kB = 1.38*10^(-23);
lambda = 633*10^(-9);
theta = 173*2*pi/360;
q = 4*pi*n*sin(theta/2)/lambda;

D = kB*T/(6*pi*eta*R);
Gamma = q^2 * D;

%% Generate Data
% t = logspace(-6,0)';
% y = 0.5*exp(-0.5*Gamma*t) + 0.5*exp(-2*Gamma*t);
% for i = 1:length(y)
%     corr = min(1/y(i),10)/1000;
%     y(i) = y(i) + corr*(2 * rand(1,1) - 1);
% end

%% Load Data
load('F:\CONTIN_CRIE_TestFolder\tBiMod.mat')
load('F:\CONTIN_CRIE_TestFolder\yBiMod.mat')
t = tBiMod;
y = yBiMod;


clear DLS_Input
%Input parameters for CRIE run
DLS_Input.Iquad = 2;
DLS_Input.Igrid = 2;
DLS_Input.Kernel = 1;
DLS_Input.Nnq = 0;
DLS_Input.Anq = 0;
DLS_Input.Neq = 0;
DLS_Input.Aeq = 0;
DLS_Input.Nneg = 1;
DLS_Input.Ny0 = 1;
DLS_Input.iwt = 4;        % No weights
DLS_Input.wt = y.^2;
DLS_Input.Nalpha = 40;
DLS_Input.alpha_lims = [1e-6,1e6];

%CRIE run
[s,~,yfit,~,info] = CRIE(t,y,Ng,0,DLS_Input);
%%
index = 27;
loglog(info.sol,info.Lx,'k--o')
plotDet = gca;
hold on
loglog([plotDet.XLim(1),info.sol(index)],[info.Lx(index),info.Lx(index)],':r',...
    [info.sol(index),info.sol(index)],[plotDet.YLim(1),info.Lx(index)],':r')
hold off
xlabel('residual norm || A x - b ||')
ylabel('solution (semi)norm || L x ||');
%%
Dh = q^2 * (kB*T)./(3 * pi * eta * info.g);
semilogx(1e9*Dh(50:end),s(50:end)./(-trapz(Dh*1e9,s)),'--*')
%% Cumulant Analysis
A = [ones(25,1),t(1:25),t(1:25).^2];
Coeff = A\log(y(1:25));
semilogy(t(1:25),y(1:25),'*')
hold on
semilogy(t(1:25),exp(A(1:25,:)*(Coeff)))
hold off
ZAvg = q^2 * (kB*T)./(3 * pi * eta * (-Coeff(2)))
PDI = 2*Coeff(3)/(Coeff(2)^2)
            

