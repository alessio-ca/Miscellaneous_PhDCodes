clear variables;close all
%% Data loading
uiopen('F:\SAXS_Fitting\Grenoble_SAXS_data analysis_2.xlsx',1)
%% Assignment
q = data(:,1);
ff = data(:,2);
loglog(q,ff)
%% Parameter declaration
ql_u = 0.9e-1;
ql_d = 0.45e-1;
Rmin = 31;
Rmax = 131;
P = @(q,R) (4*pi*(R.^3).*(sin(q.*R) - q.*R.*cos(q.*R))./(q.*R).^3).^2;  %Single particle form factor
%% Black magic
ff_an = ff(q<ql_u);
q_an = q(q<ql_u);
ff_an = ff_an(q_an>ql_d);
q_an = q_an(q_an>ql_d);
[RM,qM] = meshgrid((Rmin:Rmax)',q_an);
Kernel = P(qM,RM);

clear ILT_Input
Ng = 100;
Nl = 0;
ILT_Input.Iquad = 2;      % Simpson's rule
ILT_Input.Igrid = 1;      % Lin grid
ILT_Input.Kernel = Kernel;     % ILT Kernel
ILT_Input.g_lims = [Rmin,Rmax];
ILT_Input.Nnq = 0;        
ILT_Input.Anq = 0;
ILT_Input.Neq = 0;
ILT_Input.Aeq = 0;
ILT_Input.Nneg = 1;       % Non-negativity constraint
ILT_Input.Ny0 = 0;        % y(0)=1 condition
ILT_Input.iwt = 1;        % no weights
ILT_Input.alpha_lims = [1e2,1e13];
ILT_Input.Nalpha = 30;
ILT_Input.Nbg = 0;

[s,g,yfit,lambda,info] = CRIE(q_an,ff_an,Ng,Nl,ILT_Input);
%%
clear ILT_Input
Ng = 100;
Nl = 0;
ILT_Input.Iquad = 2;      % Simpson's rule
ILT_Input.Igrid = 1;      % Lin grid
ILT_Input.Kernel = Kernel;     % ILT Kernel
ILT_Input.g_lims = [Rmin,Rmax];
ILT_Input.Nnq = 0;        
ILT_Input.Anq = 0;
ILT_Input.Neq = 0;
ILT_Input.Aeq = 0;
ILT_Input.Nneg = 1;       % Non-negativity constraint
ILT_Input.Ny0 = 0;        % y(0)=1 condition
ILT_Input.iwt = 4;        % User weights
ILT_Input.wt = 1./(ff_an*0.01).^2;
ILT_Input.alpha_lims = [1e2,1e14];
ILT_Input.Nalpha = 30;
ILT_Input.Nbg = 0;

[s_w,g_w,yfit_w,lambda_w,info_w] = CRIE(q_an,ff_an,Ng,Nl,ILT_Input);
%%
close all
figure(1)
loglog(q_an,yfit,'b',q_an,yfit_w,'r',q_an,ff_an,'k*')
figure(2)
plot(g,s,'-b*',g_w,s_w,'-r*')
%%
ff_th_93 = 1e-10*P(q,93);
ff_th_67 = 1e-10*P(q,67);
loglog(q,ff,q,ff_th_93,q,ff_th_67,q,0.442*ff_th_93 + 0.195*ff_th_67)
%%
g_67 = g(25:50);
N_67 = s(25:50);
g_93 = g(50:75);
N_93 = s(50:75);

f_67 = fit(g_67,N_67,'gauss1')
f_93 = fit(g_93,N_93,'gauss1')
f = fit(g,s,'gauss2')
plot(f_67,g_67,N_67)
hold on
plot(f_93,g_93,N_93)
plot(f,g,s)
hold off



