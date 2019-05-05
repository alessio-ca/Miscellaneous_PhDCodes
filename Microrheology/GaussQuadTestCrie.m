clear variables;close all
%% Constants
T = 298;
eta = 0.0008872; %Fluid viscosity (in Pa.s)
kB = 1.38*10^(-23);
R = 0.115e-6; %Bead radius
rho_p = 1.04; %Bead density (in g/cm^3)
rho_f = 1.00; %Fluid intensity (in g/cm^3)
m_p = 1e3*rho_p*(4/3)*pi*R^3; %Bead mass (in Kg)
m_f = 1e3*rho_f*(4/3)*pi*R^3; %Displaced fluid mass (in Kg)
D0 = (kB*T)/(6*pi*R*eta);
J0 = (pi*R)/(kB*T);
t = logspace(-6,1,300)';
eta_Max = 500*eta;
G_Max = 1e2;
m_star = m_p + 0.5*m_f;
tau_p = m_star/(6*pi*eta*R);
tau_Max = eta_Max/G_Max;
Omega=sqrt(6*pi*G_Max*R/m_star - 0.25*(1/tau_Max - 1/tau_p)^2);
Damp_const = 0.5*(1/tau_Max + 1/tau_p); 
[MSD,~,dMSD] = Gen_Maxwell_fun(t,Omega,Damp_const,tau_Max,tau_p);
MSD = (3*kB*T/m_star)*MSD;

%% Fit to MSD (for smoothness)

% Parameters specification:
clear ILT_Input
Ng = 161;
Nl = 1;
ILT_Input.Iquad = 2;      % Simpson's rule
ILT_Input.Igrid = 2;      % Log grid
ILT_Input.Kernel = 2;     % ILT Kernel
ILT_Input.Nnq = 0;
ILT_Input.Anq = 0;
ILT_Input.Neq = 0;
ILT_Input.Aeq = 0;
ILT_Input.Nneg = 1;       % Non-negativity constraint
ILT_Input.Ny0 = 0;        % y(0)=1 condition
ILT_Input.iwt = 1;        % Unweighted
%ILT_Input.iwt = 4;       % User weights
%ILT_Input.wt = dataTemp.^2;
ILT_Input.alpha_lims = [0.0001,0.01];
ILT_Input.Nbg = 2;

[s,g,yfit,lambda,info] = CRIE_MSD(t,MSD,Ng,Nl,ILT_Input);
%%
% Parameters specification:
clear ILT_Input
Ng = 161;
Nl = 2;
ILT_Input.Iquad = 3;      % Gauss-Laguerre
ILT_Input.Igrid = 3;      % Gauss-Laguerre grid
ILT_Input.Kernel = 2;     % ILT Kernel
ILT_Input.Nnq = 0;
ILT_Input.Anq = 0;
ILT_Input.Neq = 0;
ILT_Input.Aeq = 0;
ILT_Input.Nneg = 1;       % Non-negativity constraint
ILT_Input.Ny0 = 0;        % y(0)=1 condition
ILT_Input.iwt = 1;        % Unweighted
%ILT_Input.iwt = 4;       % User weights
%ILT_Input.wt = dataTemp.^2;
ILT_Input.alpha_lims = [0.0001,0.01];
ILT_Input.Nbg = 3;

[s1,g1,yfit1,lambda1,info1] = CRIE_MSD(t,MSD,Ng,Nl,ILT_Input);