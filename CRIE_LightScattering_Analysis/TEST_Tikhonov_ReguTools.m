clear all; close all

%Noisy data
[A,b_bar,x] = shaw(32);
randn('seed',41997);
e = 1e-3*randn(size(b_bar)); 
b = b_bar + e;

%CSVD decomposition
[U,s,V] = csvd(A);

%%
%L-curve (Tikhonov regularization)
[lambda_Ti,rho_Ti,eta_Ti,reg_param_Ti] = l_curve(U,s,b);   
axis([1e-3,1,1,1e3])
[lambda_Tik,rho_Tik,eta_Tik] = l_corner(rho_Ti,eta_Ti,reg_param_Ti,U,s,b,'tikh');  

%%
%L-curve (TSVD)
[k_T,rho_T,eta_T,reg_param_T] = l_curve(U,s,b,'tsvd'); 
axis([1e-3,1,1,1e3])
[k_TSVD,rho_TSVD,eta_TSVD] = l_corner(rho_T,eta_T,reg_param_T,U,s,b,'tsvd'); 

%%
%L-curve (general corner)
[k_c,rho_c,eta_c] = l_corner(rho_Ti,eta_Ti);
[k_l,rho_l,eta_l] = l_corner(rho_T,eta_T);

%%
%L-curve: CRIE_Corner_New
[k_cr,rho_rc,eta_lc] =  CRIE_corner_New(rho_T,eta_T);

%%
%L-curve (general corner - hyperbole shape, many points)
rho_Parab = rho_Ti;
eta_Parab = exp(1./log(rho_Parab/(4.5*10^(-3))));
[k_p,rho_p,eta_p] = l_corner(rho_Parab,eta_Parab);
figure(1)
plot_lc(rho_Parab,eta_Parab,'o')
ax = axis;
HoldState = ishold; hold on;
loglog([min(rho_Parab)/100,rho_p],[eta_p,eta_p],':r',...
    [rho_p,rho_p],[min(eta_Parab)/100,eta_p],':r')
title(['L-curve corner at ',num2str(k_p)]);
axis(ax)
if (~HoldState), hold off; end

%%
%L-curve (general corner - hyperbole shape, scattered points)
rho_Parab_s = rho_Parab(1:4:end);
eta_Parab_s = exp(1./log(rho_Parab_s/(4.5*10^(-3))));
[k_p_s,rho_p_s,eta_p_s] = l_corner(rho_Parab_s,eta_Parab_s);
figure(2)
plot_lc(rho_Parab_s,eta_Parab_s,'o')
ax = axis;
HoldState = ishold; hold on;
loglog([min(rho_Parab_s)/100,rho_p_s],[eta_p_s,eta_p_s],':r',...
    [rho_p_s,rho_p_s],[min(eta_Parab_s)/100,eta_p_s],':r')
title(['L-curve corner at ',num2str(k_p_s)]);
axis(ax)
if (~HoldState), hold off; end

%%
%L-curve (CRIE_corner_New - hyperbole shape, scattered points)
[k_p_s,rho_p_s,eta_p_s] = CRIE_corner_New(rho_Parab_s,eta_Parab_s,(1:length(rho_Parab_s)),10^3);
figure(2)
plot_lc(rho_Parab_s,eta_Parab_s,'o')
ax = axis;
HoldState = ishold; hold on;
loglog([min(rho_Parab_s)/100,rho_p_s],[eta_p_s,eta_p_s],':r',...
    [rho_p_s,rho_p_s],[min(eta_Parab_s)/100,eta_p_s],':r')
title(['L-curve corner at ',num2str(k_p_s)]);
axis(ax)
if (~HoldState), hold off; end

%%
%L-curve (CRIE_corner_New - DLS data)
load('F:\CONTIN_CRIE_TestFolder\L-curve_Fix\DLSdata_Lcurve.mat')
rho_DLS = flipud(DLS_Lcurve(:,1));
eta_DLS = flipud(DLS_Lcurve(:,2));
reg_param_DLS = flipud(DLS_Lcurve(:,3));
[k_D,rho_D,eta_D] = CRIE_corner_New(rho_DLS,eta_DLS,reg_param_DLS);
plot_lc(rho_DLS,eta_DLS,'o')
ax = axis;
HoldState = ishold; hold on;
loglog([min(rho_DLS)/100,rho_D],[eta_D,eta_D],':r',...
    [rho_D,rho_D],[min(eta_DLS)/100,eta_D],':r')
title(['L-curve corner at ',num2str(k_D)]);
axis(ax)
if (~HoldState), hold off; end







