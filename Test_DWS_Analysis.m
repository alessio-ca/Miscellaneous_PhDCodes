%% Constants
T = 273+48;
eta = 0.0005679; %Fluid viscosity (in Pa.s)
n = 1.33;
kB = 1.38*10^(-23);
lambda = 685*10^(-9); %Laser wavelength
L = 2e-3; %Cuvette thickness
l_star = 401.87e-6; %Mean free path
k0 = 2*pi*n/lambda;
R = 0.115e-6; %Bead radius
rho_p = 1.04; %Bead density (in g/cm^3)
rho_f = 1.00; %Fluid intensity (in g/cm^3)
m_p = 1e3*rho_p*(4/3)*pi*R^3; %Bead mass (in Kg)
m_f = 1e3*rho_f*(4/3)*pi*R^3; %Displaced fluid mass (in Kg)
D = kB*T/(6*pi*eta*R);

%% TEST ON WATER: DWS ROUTINE
g1_an = @(x) (L/l_star + 4/3)/(5/3) .* (sinh(x) + (2/3).*x.*cosh(x)) ./ ...
    ((1 + (4/9).*x.^2).*sinh((L/l_star).*x) + (4/3).*x.*cosh((L/l_star).*x));
tDWS = logspace(-6,1,300)';
g2DWS = g1_an(sqrt(k0^2 * 6 * D * tDWS)).^2;
[tau_CON,MSD_CON]=DWS_Analysis(tDWS,g2DWS,'T',T,'R',R,'fit','CON');
[tau_RAT,MSD_RAT]=DWS_Analysis(tDWS,g2DWS,'T',T,'R',R,'fit','RAT');
[tau_SPL,MSD_SPL]=DWS_Analysis(tDWS,g2DWS,'T',T,'R',R,'fit','SPL');

close all
loglog(tau_CON,MSD_CON*1e+12,'o',tau_RAT,MSD_RAT*1e+12,'x',tau_SPL,MSD_SPL*1e+12,'*',tau_CON,6*D*tau_CON*1e+12,'k--')
legend('CON','RAT','SPL','Pure diffusive','Location','northwest')
xlabel('Time (s)')
ylabel('MSD (\mum^2)')

%% TEST ON SINGLE RELAXATION MAXWELL WITH SOLVENT DISSIPATION: DWS ROUTINE
g1_an = @(x) (L/l_star + 4/3)/(5/3) .* (sinh(x) + (2/3).*x.*cosh(x)) ./ ...
    ((1 + (4/9).*x.^2).*sinh((L/l_star).*x) + (4/3).*x.*cosh((L/l_star).*x));
tDWS = logspace(-8,1,300)';
eta_Max = 10*eta;
G_Max = 2e2;
m_star = m_p + 0.5*m_f;
tau_p = m_star/(6*pi*eta*R);
tau_f = 1e3*rho_f*R^2/eta;
tau_Max = eta_Max/G_Max;
Omega=sqrt(6*pi*G_Max*R/m_star - 0.25*(1/tau_Max - 1/tau_p)^2);
Damp_const = 0.5*(1/tau_Max + 1/tau_p); 
MSDDWS = (3*kB*T/m_star)*Gen_Maxwell_fun(tDWS,Omega,Damp_const,tau_Max,tau_p);
g1DWS = g1_an(sqrt(k0^2 * MSDDWS));
g2DWS = g1DWS.^2;
[tau_CON,MSDMAX_CON,ACFMAX_CON]=DWS_Analysis(tDWS,g2DWS,'T',T,'R',R,'fit','CON');
[tau_RAT,MSDMAX_RAT,ACFMAX_RAT]=DWS_Analysis(tDWS,g2DWS,'T',T,'R',R,'fit','RAT');
[tau_SPL,MSDMAX_SPL,ACFMAX_SPL]=DWS_Analysis(tDWS,g2DWS,'T',T,'R',R,'fit','SPL');

close all
subplot(1,2,1)
semilogx(tDWS,g1DWS,tau_CON,ACFMAX_CON,'o',tau_RAT,ACFMAX_RAT,'x',tau_SPL,ACFMAX_SPL,'*')
legend('Theory','CON','RAT','SPL')
xlabel('Time (s)')
ylabel('g_1(t)')

subplot(1,2,2)
loglog(tDWS,MSDDWS*1e+12,tau_CON,MSDMAX_CON*1e+12,'o',tau_RAT,MSDMAX_RAT*1e+12,'x',tau_SPL,MSDMAX_SPL*1e+12,'*',tDWS,6*D*tDWS*1e+12,'k--')
legend('Theory','CON','RAT','SPL','Pure diffusive','Location','northwest')
xlabel('Time (s)')
ylabel('MSD (\mum^2)')

%% TEST ON EXPERIMENTAL DATA: DWS ROUTINE
% Import data
close all
filename_ICF = 'F:\AC_2017_10_16\cooling_48C_ICF.dat';
filename_MR = 'F:\AC_2017_10_16\cooling_48C_MR.dat';
ICF = DWSread(filename_ICF,'raw');
MR = DWSread(filename_MR,'MR');

%%
% Preliminary plot (individuate the tail, if present)
semilogx(ICF(:,1),ICF(:,2))
xlabel('Time (s)')
ylabel('g_2-1 (t)')

%% 
% Analysis section (here the tail is taken from 0.1 to 1)
tail=[1e-1,1];
[tau_CON,MSD_CON,ACF_CON]=DWS_Analysis(ICF(:,1),ICF(:,2),'T',T,'tail',tail,'fit','CON');
[tau_RAT,MSD_RAT,ACF_RAT]=DWS_Analysis(ICF(:,1),ICF(:,2),'T',T,'tail',tail,'fit','RAT');
[tau_SPL,MSD_SPL,ACF_SPL]=DWS_Analysis(ICF(:,1),ICF(:,2),'T',T,'tail',tail,'fit','SPL');

%%
% Intermediate plot: g1 and MSD
close all
subplot(1,2,1)
semilogx(ICF(:,1),sqrt(ICF(:,2)),tau_CON,ACF_CON,'o',tau_RAT,ACF_RAT,'x',tau_SPL,ACF_SPL,'*')
legend('Raw','CON','RAT','SPL')
xlabel('Time (s)')
ylabel('g_1(t)')

subplot(1,2,2)
loglog(MR(:,1),MR(:,2),tau_CON,MSD_CON*1e+12,'o',tau_RAT,MSD_RAT*1e+12,'x',tau_SPL,MSD_SPL*1e+12,'*',MR(:,1),6*D*MR(:,1)*1e+12,'k--')
legend('Raw','CON','RAT','SPL','Pure diffusive','Location','northwest')
xlabel('Time (s)')
ylabel('MSD (\mum^2)')

%%
%For the spline, discard the first two order of magnitude (due to noise)
MSD_SPL_Acc=MSD_SPL(tau_SPL>100*tau_SPL(1));
tau_SPL_Acc=tau_SPL(tau_SPL>100*tau_SPL(1));

%%
% First MSD-to-G conversion: identify artifacts etc.
Jfactor = 1 / (kB*T/(pi*R));
[omega_in,G_in,~,~,DeltaA]=MSDtoG_Evans_oversampling(tau_CON,MSD_CON,1/tau_CON(1),'CG',1.001,'Beta',1,'Nf0',10,'Nfinf',4,'Jfactor',Jfactor);

%%
close all
subplot(1,2,1)
loglog(MR(:,5),MR(:,6),'b',MR(:,5),MR(:,7),'r',omega_in,real(G_in),'*-',omega_in,imag(G_in),'x-')
legend('G'',Raw','G'''',Raw','G'', Prelim','G'''', Prelim','location','northwest')
xlabel('\omega (rad/s)') 
ylabel('Complex Modulus (Pa)')
subplot(1,2,2)
loglog(DeltaA)
xlabel('Pts')
ylabel('\DeltaA')

%%
%From Delta_A plot, we see that residuals (Delta_A) go to zero after ~10^4
%points but display an artificial plateau at Delta_A~5e-5. A filter Eps is
%therefore needed. Plus, the last residual display a sudden jump. "Black
%magic" will be required.
Eps = 5e-5;
[omega_fin,G_fin]=MSDtoG_Evans_oversampling(tau_CON,MSD_CON,1/tau_CON(1),'CG',1.001,'Beta',10,'Nf0',10,'Nfinf',4,'Jfactor',Jfactor,'Eps',Eps,'BM',3);

%%
close all
loglog(MR(:,5),MR(:,6),'b',MR(:,5),MR(:,7),'r',omega_in,real(G_in),'k-',omega_in,imag(G_in),'k--',omega_fin,real(G_fin),'*-',omega_fin,imag(G_fin),'x-')
legend('G'', LS','G'''', LS','G'', Prelim','G'''', Prelim','G'', Fin','G'''', Fin','location','northwest')
xlabel('\omega (rad/s)')
ylabel('Complex Modulus (Pa)')

%%
%For comparison, the G* is calculated for all the series
[omega_RAT,G_RAT]=MSDtoG_Evans_oversampling(tau_RAT,MSD_RAT,1/tau_RAT(1),'CG',1.001,'Beta',10,'Nf0',10,'Nfinf',4,'Jfactor',Jfactor,'BM',1);
[omega_SPL,G_SPL]=MSDtoG_Evans_oversampling(tau_SPL_Acc,MSD_SPL_Acc,1/tau_SPL_Acc(1),'CG',1.001,'Beta',10000,'Nf0',10,'Nfinf',4,'Jfactor',Jfactor,'Eps',2e-6,'BM',1);
%%
loglog(MR(:,5),MR(:,6),'b',MR(:,5),MR(:,7),'r')
hold on
loglog(omega_fin,real(G_fin),'o-',omega_fin,imag(G_fin),'x-')
loglog(omega_RAT,real(G_RAT),'o-',omega_RAT,imag(G_RAT),'x-')
loglog(omega_SPL,real(G_SPL),'o-',omega_SPL,imag(G_SPL),'x-')
hold off
legend('G'', LS','G'''', LS','G'', CON','G'''', CON','G'', RAT','G'''', RAT','G'', SPL','G'''', SPL','location','northwest')
xlabel('\omega (rad/s)')
ylabel('Complex Modulus (Pa)')

%Conclusions: the CON and SPL pretty much agree, with the CON covering a
%larger bandwidth. The RAT for this case is totally off, so not trustworthy
%(we could have reached this conclusion by judging the ACF and MSD plots
%earlier). The SPL is (often) correct, but in this case the CON is 
%equivalent and spans a larger frequenct range. The CON result is therefore 
%chosen as the FINAL one. The G' result is accurate up to w=10^6 rad/s 
%(before the little dip). The G''result is accurate up to the whole range.
