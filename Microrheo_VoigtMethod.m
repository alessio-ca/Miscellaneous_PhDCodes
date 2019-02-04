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
%% Test 1: H2O (with noise, Voigt-fluid fit method)
load('C:\Users\ac2014\Documents\MATLAB\LDdiff_H2O_10runs_err.mat')
Nerr = (length(msd):-1:1)';
tau = tau(1:floor(length(tau)/10));
tau(1)=[];
msd(1)=[];
msderr(1)=[];
Nerr(1)=[];
LogInt = unique(round(logspace(0,log10(length(tau)),300)))';
msd=msd(LogInt);
tau=tau(LogInt);
msderr=msderr(LogInt);
Nerr=Nerr(LogInt);
msdstd = msd./sqrt(Nerr);
%%
% Parameters specification:
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
ILT_Input.iwt = 4;        % User weights
ILT_Input.wt = 1./msderr.^2;
ILT_Input.alpha_lims = [0.01,100];
ILT_Input.Nbg = 2;

[s,g,yfit,lambda,info] = CRIE(tau,msd,Ng,Nl,ILT_Input);

clear ILT_Input
% Parameters specification:
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
ILT_Input.alpha_lims = [0.01,100];
ILT_Input.Nbg = 2;

[s_ne,g_ne,yfit_ne,lambda_ne,info_ne] = CRIE(tau,msd,Ng,Nl,ILT_Input);
%%
close all
subplot(1,2,1)
plot(tau,yfit,tau,yfit_ne,tau,msd,tau,2*D0*tau,'k--')
subplot(1,2,2)
loglog(tau,abs(yfit./(2*D0*tau) - 1),tau,abs(yfit_ne./(2*D0*tau) - 1),tau,abs(msd./(2*D0*tau) - 1))
%%
Om = logspace(log10(1/tau(end-3)),log10(1/tau(1)),100)';

s_Ref = abs(s);
s_Ref_ne = abs(s_ne);
s_Ref(s_Ref < 1e-6*max(s))=0; %Eliminate numerical noise in s (by clipping to 0)
s_Ref_ne(s_Ref_ne < 1e-6*max(s_ne))=0; %Eliminate numerical noise in s (by clipping to 0)

[gM,~] = meshgrid(info.g,Om);
[cM,OmM] = meshgrid(info.c,Om);
[gM_ne,~] = meshgrid(info_ne.g,Om);
[cM_ne,OmM_ne] = meshgrid(info_ne.c,Om);

A_Om = cM.*(log(10)*10.^(gM)).*1./(1 + 1i*OmM./(10.^gM));
y_Om = A_Om*s_Ref(1:end-Nl) - 1i*s_Ref(end)./Om;
A_Om_ne = cM_ne.*(log(10)*10.^(gM_ne)).*1./(1 + 1i*OmM./(10.^gM_ne));
y_Om_ne = A_Om_ne*s_Ref_ne(1:end-Nl) - 1i*s_Ref_ne(end)./Om;

G = kB*T./(3*pi*R*y_Om);
G_ne = kB*T./(3*pi*R*y_Om_ne);
%%
loglog(Om,imag(G),Om,imag(G_ne),Om,8.9e-4*Om,'k--')
Gfit = Om;
eta = Gfit\imag(G);
eta_ne = Gfit\imag(G_ne);

%%
[OmMas,GMas]=MSDtoG_Mason(tau,3*msd,'R',115e-9,'CG',1.01,'T',T,'cutoff',0);
%%
[OmEv,GEv]=MSDtoG_Evans_oversampling(tau,msd,1e6,'R',115e-9,'CG',1.01,'T',T,'Beta',1,'Jfactor',3*J0);
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
%set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(get(get(h5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%set(get(get(h7,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('G''  - KC','G''''  - KC','G'' - Mason','G''''  - Mason','G''  - Evans','G''''  - Evans','G''  - Theory','G''''  - Theory','Location','northwest')



