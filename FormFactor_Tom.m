clear variables;close all
%% Data loading
uiopen('F:\SAXS_Fitting\FF_data_for_Tom.xlsx',1)
%% Assignment
q = data(:,1);
ff = data(:,2);
loglog(q,ff)
%% Parameter declaration
ql_d = 0.3e-1;
ql_u = 0.8e-1;
Rmin = 65;
Rmax = 100;
P = @(q,R) q.^4.*(4*pi*(R.^3).*(sin(q.*R) - q.*R.*cos(q.*R))./(q.*R).^3).^2;  %Single particle form factor
%% Black magic
ff_an = ff(q<ql_u);
q_an = q(q<ql_u);
ff_an = ff_an(q_an>ql_d);
q_an = q_an(q_an>ql_d);
[RM,qM] = meshgrid(linspace(Rmin,Rmax,101)',q_an);
Kernel = P(qM,RM);

clear ILT_Input
Ng = 101;
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
ILT_Input.alpha_lims = [1e4,1e7];
ILT_Input.Nalpha = 30;
ILT_Input.Nbg = 0;

[s,g,yfit,lambda,info] = CRIE(q_an,ff_an.*q_an.^4,Ng,Nl,ILT_Input);

%%
c = [4*ones(1,(Ng-3)/2);2*ones(1,(Ng-3)/2)]; %Simpson rule for coeff
c = c(:);
h = (info.g(end) - info.g(1))/(Ng-1);
c = (h/3)*[1;c;4;2];
t = logspace(-2.5,-0.9,200)';
[cM,tM] = meshgrid(c,t);
%%
[RM,qM] = meshgrid(linspace(Rmin,Rmax,101)',t);
Kernel = P(qM,RM);
A = cM.*Kernel;
%%
loglog(q,ff,t,(A*s)./t.^4)
%%
qTh = t;
ffTh = (A*s)./t.^4;
filename = 'formfactorTom.xlsx';
xlswrite(filename,{'Wavenumber','Form Factor (0.1 %)','Wavenumber Th','Form Factor Th'},'A1:D1')
xlswrite(filename,q,1,'A2')
xlswrite(filename,ff,1,'B2')
xlswrite(filename,qTh,1,'C2')
xlswrite(filename,ffTh,1,'D2')
%%




