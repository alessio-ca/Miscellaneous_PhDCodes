clear all
close all
clc

Coeff_A = xlsread('L:\OT_ParticleWall\FitCoeff.xlsx','B3:G6');
Coeff_B = xlsread('L:\OT_ParticleWall\FitCoeff.xlsx','B10:G13');

A = 2.414*1e-5;
B = 247.8;
C = 140;

kB = 1.38*1e-23;
R = 0.75*1e-6;
eta = @(x) A*10.^(B./(x-C));

T = 298;
GammaTh = (6*pi*eta(T)*R);
Dth = kB*T/GammaTh;

Far_A_x = Coeff_A(1,[1,3,5]);
Far_A_y = Coeff_A(1,[2,4,6]);
Far_B_x = Coeff_B(1,[1,3,5]);
Far_B_y = Coeff_B(1,[2,4,6]);
Med_A_x = Coeff_A(2,[1,3,5]);
Med_A_y = Coeff_A(2,[2,4,6]);
Med_B_x = Coeff_B(2,[1,3,5]);
Med_B_y = Coeff_B(2,[2,4,6]);
MCl_A_x = Coeff_A(3,[1,3,5]);
MCl_A_y = Coeff_A(3,[2,4,6]);
MCl_B_x = Coeff_B(3,[1,3,5]);
MCl_B_y = Coeff_B(3,[2,4,6]);
Clo_A_x = Coeff_A(4,[1,3,5]);
Clo_A_y = Coeff_A(4,[2,4,6]);
Clo_B_x = Coeff_B(4,[1,3,5]);
Clo_B_y = Coeff_B(4,[2,4,6]);

%% Temperature estimation
D_x = 10^(Far_A_x(1))/(2*Far_B_x(1));
D_y = 10^(Far_A_y(1))/(2*Far_B_y(1));

Tspace = linspace(273,323);
plot(Tspace,Tspace./eta(Tspace))
hold on
plot(Tspace,(6*pi*R*D_y/kB)*ones(size(Tspace)))
hold off
xlabel('T [K]')
ylabel('T/\eta [K/(Pa*s)]')
%%
T = 298.5;
GammaTh = (6*pi*eta(T)*R);
Dth = kB*T/GammaTh;
%% Trap constants VS T
Kx_Far = GammaTh./Far_B_x;
Ky_Far = GammaTh./Far_B_y;
Kx_Far_A = 2*kB*T./10.^Far_A_x;
Ky_Far_A = 2*kB*T./10.^Far_A_y;

G = [0.0025,0.005,0.01];
plot(G,Kx_Far,'--*',G,Ky_Far,'--*')
xlabel('Gain')
ylabel('k_{trap} [N/m]')
%% Temperature correction for same height series: Far
Tspace = linspace(273,323);
TEtaX = ones(size(Tspace))'*(3*pi*R*10.^Far_A_x./(kB*Far_B_x));
TEtaY = ones(size(Tspace))'*(3*pi*R*10.^Far_A_y./(kB*Far_B_y));

plot(Tspace,Tspace./eta(Tspace))
hold on
plot(Tspace,TEtaX(:,1),'b-',Tspace,TEtaX(:,2),'b--',Tspace,TEtaX(:,3),'b-.')
plot(Tspace,TEtaY(:,1),'r-',Tspace,TEtaY(:,2),'r--',Tspace,TEtaY(:,3),'r-.')
hold off
xlabel('T [K]')
ylabel('T/\eta [K/(Pa*s)]')
%%
TCorr = [T, 298.75 , 297.75];
GammaCorr = [GammaTh,(6*pi*eta(TCorr(2))*R),(6*pi*eta(TCorr(3))*R)];

Kx_A_C = 2*kB*TCorr./10.^Far_A_x;
Ky_A_C = 2*kB*TCorr./10.^Far_A_y;
Kx_C = GammaCorr./Far_B_x;
Ky_C = GammaCorr./Far_B_y;

disp(['Ky ratios (Corrected): ',num2str(Ky_C./Ky_C(1))])
disp(['Ky_A ratios (Corrected): ',num2str(Ky_A_C./Ky_A_C(1))])

%% Gamma values (assuming T constant)
T=298;
Kx_Far_A = 2*kB*T./10.^Far_A_x;
Ky_Far_A = 2*kB*T./10.^Far_A_y;
GammaFar_x = Kx_Far_A.*Far_B_x; 
GammaFar_y = Ky_Far_A.*Far_B_y; 


Kx_Med_A = 2*kB*T./10.^Med_A_x;
Ky_Med_A = 2*kB*T./10.^Med_A_y;
GammaMed_x = Kx_Med_A.*Med_B_x; 
GammaMed_y = Ky_Med_A.*Med_B_y; 

Kx_MCl_A = 2*kB*T./10.^MCl_A_x;
Ky_MCl_A = 2*kB*T./10.^MCl_A_y;
GammaMCl_x = Kx_MCl_A.*MCl_B_x; 
GammaMCl_y = Ky_MCl_A.*MCl_B_y; 

Kx_Clo_A = 2*kB*T./10.^Clo_A_x;
Ky_Clo_A = 2*kB*T./10.^Clo_A_y;
GammaClo_x = Kx_Clo_A.*Clo_B_x; 
GammaClo_y = Ky_Clo_A.*Clo_B_y;

%% Faxen
distList = [15,10,3,1.5] + R*1e6;

GammaWeak_x = [GammaFar_x(1),GammaMed_x(1),GammaMCl_x(1),GammaClo_x(1)];
GammaWeak_y = [GammaFar_y(1),GammaMed_y(1),GammaMCl_y(1),GammaClo_y(1)];
GammaMod_x = [GammaFar_x(2),GammaMed_x(2),GammaMCl_x(2),GammaClo_x(2)];
GammaMod_y = [GammaFar_y(2),GammaMed_y(2),GammaMCl_y(2),GammaClo_y(2)];
GammaStrong_x = [GammaFar_x(3),GammaMed_x(3),GammaMCl_x(3),GammaClo_x(3)];
GammaStrong_y = [GammaFar_y(3),GammaMed_y(3),GammaMCl_y(3),GammaClo_y(3)];

FaxList = logspace(0,1.5);
Fax = @(x) 1./(1-(9/16)*(R*1e6./x) +(1/8)*(R*1e6./x).^3);

loglog(distList/(R*1e6),GammaWeak_x./GammaTh,'bo--','MarkerSize',8)
hold on
loglog(distList/(R*1e6),GammaWeak_y./GammaTh,'bo--','MarkerSize',8,'MarkerFaceColor','b')
loglog(distList/(R*1e6),GammaMod_x./GammaTh,'ro--','MarkerSize',8)
loglog(distList/(R*1e6),GammaMod_y./GammaTh,'ro--','MarkerSize',8,'MarkerFaceColor','r')
loglog(distList/(R*1e6),GammaStrong_x./GammaTh,'o--','Color',[0,0.5,0],'MarkerSize',8)
loglog(distList/(R*1e6),GammaStrong_y./GammaTh,'o--','Color',[0,0.5,0],'MarkerSize',8,'MarkerFaceColor',[0,0.5,0])
loglog(FaxList/(R*1e6),Fax(FaxList),'k-')
hold off
xlabel('S/R')
xlim([1 30])
ylim([0.9 1.5])
ylabel('\gamma^{exp}/\gamma^{Stokes}')
legend('Weak, X','Weak, Y','Medium, X','Medium, Y','Strong, X','Strong, Y')
set(gca,'FontSize',14)















