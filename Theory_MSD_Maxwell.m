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



%% Test 1: H2O
MSD = 6 * D0 * t; % 3D diffusivity
dMSD = 6 * D0 * ones(length(t),1);
figure(1)
yyaxis left
loglog(t,MSD)
yyaxis right
loglog(t,dMSD)

dMSD_Fit = dMSD;
t_Fit = t;
t_Fit(dMSD_Fit<1e-16)=[];
dMSD_Fit(dMSD_Fit<1e-16)=[];

% Parameters specification:
Ng = 161;
Nl = 1;
ILT_Input.Iquad = 2;      % Simpson's rule
ILT_Input.Igrid = 2;      % Log grid
ILT_Input.Kernel = 1;     % ILT Kernel
ILT_Input.Nnq = 0;        
ILT_Input.Anq = 0;
ILT_Input.Neq = 0;
ILT_Input.Aeq = 0;
ILT_Input.Nneg = 1;       % Non-negativity constraint
ILT_Input.Ny0 = 0;        % y(0)=1 condition
ILT_Input.iwt = 1;        % Unweighted analysis
ILT_Input.alpha_lims = [0.01,100];
ILT_Input.Nbg = 1;

[s,g,yfit,lambda,info] = CRIE(t_Fit,dMSD_Fit,Ng,Nl,ILT_Input);

bg_term = s(end);
fit_term = info.A(:,1:end-1)*s(1:end-1);

if all(abs(fit_term)./bg_term<1e-3) 
    s(1:end-1,:) = 0;
end
s = abs(s);
Om = logspace(log10(1/t(end-3)),log10(1/t(1)),100)';
[gM,~] = meshgrid(info.g,Om);
[cM,OmM] = meshgrid(info.c,Om);
A = cM.*(log(10)*10.^(gM)).*1./(1i*OmM.*(10.^gM + 1i*OmM));
y = A*s(1:end-1) - s(end)./Om.^2;
G = 3*kB*T./(3*pi*R*1i*Om.*y);
%%
[OmMas,GMas]=MSDtoG_Mason(t,MSD,'R',115e-9,'CG',1.01,'T',T);
%%
[OmEv,GEv]=MSDtoG_Evans_oversampling(t,MSD,1e6,'R',115e-9,'CG',1.01,'T',T,'Beta',1,'Jfactor',J0);
%%
Colours = get(groot,'DefaultAxesColorOrder');
h1=loglog(Om,real(G),'o','MarkerEdgeColor',Colours(1,:),'MarkerFaceColor',Colours(1,:));
hold on
h2=loglog(Om,imag(G),'o','MarkerEdgeColor',Colours(1,:));
h3=loglog(OmMas,real(GMas),'s','MarkerEdgeColor',Colours(2,:),'MarkerFaceColor',Colours(2,:));
h4=loglog(OmMas,imag(GMas),'s','MarkerEdgeColor',Colours(2,:));
h5=loglog(OmEv,real(GEv),'h','MarkerEdgeColor',Colours(3,:),'MarkerFaceColor',Colours(3,:));
h6=loglog(OmEv,imag(GEv),'h','MarkerEdgeColor',Colours(3,:));
h7=loglog(Om,zeros(length(Om),1),'k-');
h8=loglog(Om,eta*Om,'k--');
hold off
xlabel('\omega [rad/s]')
ylabel('G'',G'''' [Pa]')
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(get(get(h5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h7,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('G''''  - KC','G''''  - Mason','G''  - Evans','G''''  - Evans','G''''  - Theory','Location','northwest')

%% Test 2: Viscoelastic Material
eta_Max = 500*eta;
G_Max = 1e2;
m_star = m_p + 0.5*m_f;
tau_p = m_star/(6*pi*eta*R);
tau_Max = eta_Max/G_Max;
Omega=sqrt(6*pi*G_Max*R/m_star - 0.25*(1/tau_Max - 1/tau_p)^2);
Damp_const = 0.5*(1/tau_Max + 1/tau_p); 
[MSD,~,dMSD] = Gen_Maxwell_fun(t,Omega,Damp_const,tau_Max,tau_p);
MSD = (3*kB*T/m_star)*MSD;
dMSD = (3*kB*T/m_star)*dMSD;
figure(1)
yyaxis left
loglog(t,MSD)
yyaxis right
loglog(t,dMSD)

dMSD_Fit = dMSD;
t_Fit = t;
t_Fit(dMSD_Fit<1e-16)=[];
dMSD_Fit(dMSD_Fit<1e-16)=[];

% Parameters specification:
Ng = 161;
Nl = 1;
ILT_Input.Iquad = 2;      % Simpson's rule
ILT_Input.Igrid = 2;      % Log grid
ILT_Input.Kernel = 1;     % ILT Kernel
ILT_Input.Nnq = 0;        
ILT_Input.Anq = 0;
ILT_Input.Neq = 0;
ILT_Input.Aeq = 0;
ILT_Input.Nneg = 1;       % Non-negativity constraint
ILT_Input.Ny0 = 0;        % y(0)=1 condition
ILT_Input.iwt = 1;        % Unweighted analysis
ILT_Input.alpha_lims = [0.01,100];
ILT_Input.Nbg = 1;

[s,g,yfit,lambda,info] = CRIE(t_Fit,dMSD_Fit,Ng,Nl,ILT_Input);

bg_term = s(end);
fit_term = info.A(:,1:end-1)*s(1:end-1);

if all(abs(fit_term)./bg_term<1e-3) 
    s(1:end-1,:) = 0;
end
s = abs(s);
Om = logspace(log10(1/t(end-3)),log10(1/t(1)),100)';
[gM,~] = meshgrid(info.g,Om);
[cM,OmM] = meshgrid(info.c,Om);
A = cM.*(log(10)*10.^(gM)).*1./(1i*OmM.*(10.^gM + 1i*OmM));
y = A*s(1:end-1) - s(end)./Om.^2;
G = 3*kB*T./(3*pi*R*1i*Om.*y);
%%
[OmMas,GMas]=MSDtoG_Mason(t,MSD,'R',115e-9,'CG',1.01,'T',T);
%%
[OmEv,GEv]=MSDtoG_Evans_oversampling(t,MSD,1e6,'R',115e-9,'CG',1.01,'T',T,'Beta',1,'Jfactor',J0);
%%
GJ = @(x,etaI,GI,tauM) -1i.*x.*(etaI + (GI*tauM)./(1 - 1i*x*tauM));
Colours = get(groot,'DefaultAxesColorOrder');
h1=loglog(Om,real(G),'o','MarkerEdgeColor',Colours(1,:),'MarkerFaceColor',Colours(1,:));
hold on
h2=loglog(Om,imag(G),'o','MarkerEdgeColor',Colours(1,:));
h3=loglog(OmMas,real(GMas),'s','MarkerEdgeColor',Colours(2,:),'MarkerFaceColor',Colours(2,:));
h4=loglog(OmMas,imag(GMas),'s','MarkerEdgeColor',Colours(2,:));
h5=loglog(OmEv,real(GEv),'h','MarkerEdgeColor',Colours(3,:),'MarkerFaceColor',Colours(3,:));
h6=loglog(OmEv,imag(GEv),'h','MarkerEdgeColor',Colours(3,:));
h7=loglog(Om,real(GJ(Om,eta,G_Max,tau_Max)),'k-');
h8=loglog(Om,-imag(GJ(Om,eta,G_Max,tau_Max)),'k--');

OmGuideLin1 = logspace(-1,1);
OmGuideLin2 = logspace(0,2);

loglog(OmGuideLin1,10000*eta*OmGuideLin1,'k--');
loglog(OmGuideLin2,0.1*eta*OmGuideLin2.^2,'k--');


hold off
xlabel('\omega [rad/s]')
ylabel('G'',G'''' [Pa]')
%set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(get(get(h5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(get(get(h7,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('G''  - KC','G''''  - KC','G'' - Mason','G''''  - Mason','G''  - Evans','G''''  - Evans','G''  - Theory','G''''  - Theory','Location','northwest')
%% Test 3: H2O (with noise)
load('C:\Users\ac2014\Documents\MATLAB\LDdiff_H2O_10runs.mat')

lambda = 633*10^(-9);
theta = 173*2*pi/360;
n = 1.33;
q = 4*pi*n*sin(theta/2)/lambda;
g1 = exp(-q^2*msd/2);
tau = tau(1:floor(length(tau)/10));
g1 = g1(1:length(tau));
tau(1)=[];
msd(1)=[];
g1(1)=[];
LogInt = unique(round(logspace(0,log10(length(g1)),300)))';
msd=msd(LogInt);
tau=tau(LogInt);
g1=g1(LogInt);

% Parameters specification:
Ng = 161;
Nl = 0;
ILT_Input.Iquad = 2;      % Simpson's rule
ILT_Input.Igrid = 2;      % Log grid
ILT_Input.Kernel = 1;     % ILT Kernel
ILT_Input.Nnq = 0;        
ILT_Input.Anq = 0;
ILT_Input.Neq = 0;
ILT_Input.Aeq = 0;
ILT_Input.Nneg = 1;       % Non-negativity constraint
ILT_Input.Ny0 = 1;        % y(0)=1 condition
ILT_Input.iwt = 1;        % Unweighted analysis
ILT_Input.alpha_lims = [0.01,100];
ILT_Input.Nbg = 1;

[s,g,yfit,lambda,info] = CRIE(tau,g1,Ng,Nl,ILT_Input);
%%
s_Ref = s;  
s_Ref(s_Ref < 1e-6*max(s))=0; %Eliminate numerical noise in s (by clipping to 0)
s_Ref_Ex = s_Ref;
s_Ref_Ex(s_Ref_Ex < max(s_Ref))=0;
dMSD_KC = (6/q^2) * info.A * (s_Ref.*g) ./ (info.A * (s_Ref));
dMSD_KC_Ex = (6/q^2) * info.A * (s_Ref_Ex.*g) ./ (info.A * (s_Ref_Ex));
dMSD_KC = rmmissing(dMSD_KC); %Remove NaNs (due to dividing per zero)
dMSD_KC(dMSD_KC==0)=[]; %Remove zeros (unphysical)
dMSD_KC_Ex = rmmissing(dMSD_KC_Ex); %Remove NaNs (due to dividing per zero)
dMSD_KC_Ex(dMSD_KC_Ex==0)=[]; %Remove zeros (unphysical)
%%
tau_KC = tau(1:length(dMSD_KC));
tau_KC_Ex = tau(1:length(dMSD_KC_Ex));

figure(1)
semilogx(tau,yfit,tau,info.A * abs(s_Ref),'*',tau,info.A * abs(s_Ref_Ex),'o')
figure(2)
yyaxis left
loglog(tau,3*msd,'-',tau,6*D0*tau,'k--')
yyaxis right
loglog(tau_KC,dMSD_KC,tau_KC_Ex,dMSD_KC_Ex,tau(1:end-1),3*diff(msd)./diff(tau))
%%
% Parameters specification:
Ng = 161;
Nl = 2;
ILT_Input.Iquad = 2;      % Simpson's rule
ILT_Input.Igrid = 2;      % Log grid
ILT_Input.Kernel = 1;     % ILT Kernel
ILT_Input.Nnq = 0;        
ILT_Input.Anq = 0;
ILT_Input.Neq = 0;
ILT_Input.Aeq = 0;
ILT_Input.Nneg = 1;       % Non-negativity constraint
ILT_Input.Ny0 = 0;        % y(0)=1 condition
ILT_Input.iwt = 1;        % Unweighted analysis
ILT_Input.alpha_lims = [0.01,100];
ILT_Input.Nbg = 1;

[s,g,yfit,lambda,info] = CRIE(tau_KC,dMSD_KC,Ng,Nl,ILT_Input);
[s_Ex,g_Ex,yfit_Ex,lambda_Ex,info_Ex] = CRIE(tau_KC_Ex,dMSD_KC_Ex,Ng,Nl,ILT_Input);
%%
Om = logspace(log10(1/t(end-3)),log10(1/t(1)),100)';

bg_term = s(end);
fit_term = info.A(:,1:end-Nl)*s(1:end-Nl);

if all(abs(fit_term)./bg_term<1e-3) 
    s(1:end-Nl,:) = 0;
end
s = abs(s);
[gM,~] = meshgrid(info.g,Om);
[cM,OmM] = meshgrid(info.c,Om);

A = cM.*(log(10)*10.^(gM)).*1./(1i*OmM.*(10.^gM + 1i*OmM));
y = A*s(1:end-Nl) - s(end)./Om.^2;
y_Sum = A*s(1:end-Nl);
y_Const =  - s(end)./Om.^2;
G = 3*kB*T./(3*pi*R*1i*Om.*y);
G_Sum = 3*kB*T./(3*pi*R*1i*Om.*y_Sum);
G_Const = 3*kB*T./(3*pi*R*1i*Om.*y_Const);



bg_term_Ex = s_Ex(end);
fit_term_Ex = info_Ex.A(:,1:end-Nl)*s_Ex(1:end-Nl);

if all(abs(fit_term_Ex)./bg_term_Ex<1e-3)
    s_Ex(1:end-Nl,:) = 0;
end
s_Ex = abs(s_Ex);
[gM_Ex,~] = meshgrid(info_Ex.g,Om);
[cM_Ex,OmM_Ex] = meshgrid(info_Ex.c,Om);

A_Ex = cM_Ex.*(log(10)*10.^(gM_Ex)).*1./(1i*OmM_Ex.*(10.^gM_Ex + 1i*OmM_Ex));
y_Ex = A_Ex*s_Ex(1:end-1) - s_Ex(end)./Om.^2;
G_Ex = 3*kB*T./(3*pi*R*1i*Om.*y_Ex);
%%
[OmMas,GMas]=MSDtoG_Mason(tau,3*msd,'R',115e-9,'CG',1.01,'T',T);
%%
[OmEv,GEv]=MSDtoG_Evans_oversampling(tau,3*msd,1e6,'R',115e-9,'CG',1.01,'T',T,'Beta',1,'Jfactor',J0);
%%
Colours = get(groot,'DefaultAxesColorOrder');
h1=loglog(Om,real(G),'o','MarkerEdgeColor',Colours(1,:),'MarkerFaceColor',Colours(1,:));
hold on
h2=loglog(Om,imag(G),'o','MarkerEdgeColor',Colours(1,:));
h3=loglog(OmMas,real(GMas),'s','MarkerEdgeColor',Colours(2,:),'MarkerFaceColor',Colours(2,:));
h4=loglog(OmMas,imag(GMas),'s','MarkerEdgeColor',Colours(2,:));
h5=loglog(OmEv,real(GEv),'h','MarkerEdgeColor',Colours(3,:),'MarkerFaceColor',Colours(3,:));
h6=loglog(OmEv,imag(GEv),'h','MarkerEdgeColor',Colours(3,:));
h9=loglog(Om,real(G_Sum),'o','MarkerEdgeColor',Colours(4,:),'MarkerFaceColor',Colours(4,:));
h10=loglog(Om,imag(G_Sum),'o','MarkerEdgeColor',Colours(4,:));
h11=loglog(Om,real(G_Const),'o','MarkerEdgeColor',Colours(5,:),'MarkerFaceColor',Colours(5,:));
h12=loglog(Om,imag(G_Const),'o','MarkerEdgeColor',Colours(5,:));
h7=loglog(Om,zeros(length(Om),1),'k-');
h8=loglog(Om,eta*Om,'k--');
hold off
xlabel('\omega [rad/s]')
ylabel('G'',G'''' [Pa]')
%set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(get(get(h5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(get(get(h7,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('G''  - KC','G''''  - KC','G'' - Mason','G''''  - Mason','G''  - Evans','G''''  - Evans','G''  - Theory','G''''  - Theory','Location','northwest')


