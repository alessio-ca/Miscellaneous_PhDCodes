% Microrheology with OT using QPD.
%
% Calculates the complex modulus G from QPD-acquired position data of a
% single bead.

% Clear everything
clear all;
close all;
clc;
%%
%Load data & substract mean
%x = dataTot(:,1:5) - mean(dataTot(:,1:5));
%y = dataTot(:,6:10) - mean(dataTot(:,6:10));
x = reff(:,1:10) - mean(reff(:,1:10));
y = reff(:,11:20) - mean(reff(:,11:20));

%Coefficients
dt=1e-4;
fps=1/dt;

%Calculates MSD & PACF & plot for quality assessment
[acf_x,~]=acf_routine(dt,x);
[acf_y,tau]=acf_routine(dt,y);
subplot(1,2,1)
semilogx(tau,acf_x,'bo')
hold on
loglog(tau,acf_y,'ro')
hold off
legend('PACF,x','PACF,y')


[msd_x,~]=msd_routine(dt,x);
[msd_y,tau]=msd_routine(dt,y);
subplot(1,2,2)
loglog(tau,msd_x,'bo')
hold on
loglog(tau,msd_y,'ro')
hold off
legend('MSD,x','MSD,y')
%%
close all
%Calculates Gx: select proper parameters!
TauMax_x = 1;
CG_x = 1.45;
Beta_x = 500;
[omega,Gx]=PACFtoG_Evans_oversampling(tau,acf_x,fps,'TauMax',TauMax_x,'CG',CG_x,'Beta',Beta_x);
%%
%Calculates Gy: select proper parameters!

TauMax_y = 1;
CG_y = 1.45;
Beta_y = 500;
[omega,Gy]=PACFtoG_Evans_oversampling(tau,acf_y,fps,'TauMax',TauMax_y,'CG',CG_y,'Beta',Beta_y);

%% Plot everything
loglog(omega,real(Gx),'bo')
hold on
loglog(omega,imag(Gx),'bs')
loglog(omega,real(Gy),'ro')
loglog(omega,imag(Gy),'rs')
hold off

legend('G'',x','G'''',x','G'',y','G'''',y')

%% Optional: export data (change names!)
Tau_syn_eff_halved = tau;
Acf_X_syn_eff_halved = acf_x;
Acf_Y_syn_eff_halved = acf_y;
Msd_X_syn_eff_halved = msd_x;
Msd_Y_syn_eff_halved = msd_y;
Omega_syn_eff_halved = omega;
G_X_syn_eff_halved = Gx;
G_Y_syn_eff_halved = Gy;

save('Microrheo_10kHz_syn_eff_halved.mat',...
    'Tau_syn_eff_halved', 'Omega_syn_eff_halved', ...
    'Acf_X_syn_eff_halved','Acf_Y_syn_eff_halved', ...
    'Msd_X_syn_eff_halved','Msd_Y_syn_eff_halved',...
    'G_X_syn_eff_halved','G_Y_syn_eff_halved')


