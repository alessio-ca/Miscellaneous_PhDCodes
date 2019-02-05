t = t';
data = data';

%% Constants
T = 298 + 5;
n = 1.33;
R = 115e-9;
Beta = 1;
cut = 0.1;
tail = [6e5,3e6];

kB = 1.38*10^(-23);
lambda = 633*10^(-9);
theta = 173*2*pi/360;
q = 4*pi*n*sin(theta/2)/lambda;
%% g1(tau) calc. 

%Delete first 6 points (usually noise) & format data
Time_ToFit = t(6:end);
ACF_ToFit = data(6:end);
if all(tail)
    tailFit = ACF_ToFit(Time_ToFit > tail(1) & Time_ToFit < tail(2));
    tailFit = mean(tailFit);
    ACF_ToFit = (ACF_ToFit - tailFit)/(1+tailFit);
end
CritPts = find(ACF_ToFit< cut*ACF_ToFit(1));
tTemp = Time_ToFit(1:CritPts(1));

ACF_ToFit_Sqrt = sqrt(ACF_ToFit(1:CritPts(1)));
%g1(tau) calculation with  optional head treatment (least-Squares)
switch Beta
    case 1
        dataTempLog = log(ACF_ToFit_Sqrt(1:10));
        A = [ones(length(tTemp(1:10)),1),tTemp(1:10)];
        Coeff = A\dataTempLog;
        beta = exp(Coeff(1));
        
    case 0
        beta = 1;
        
    otherwise
        error('Error: Unrecognized value for beta.');
        
end
dataTemp = ACF_ToFit_Sqrt/beta;
MSDTemp = (6/q^2)*(-log(dataTemp));

%Plotting section
f=figure(1);
semilogx(t,data(:,1),'o')
hold on
semilogx(tTemp,ACF_ToFit_Sqrt,'d')
semilogx(tTemp,dataTemp,'s')
semilogx(t,cut*ones(length(t),1),'--k');
semilogx(linspace(0,tTemp(10))',exp([ones(100,1),linspace(0,tTemp(10))']*Coeff),'-.k')
semilogx(t,ones(length(t),1),'-k');
legend('g2 (Raw)','g1 (Raw)','g1(Normalized)','Location','northeast')
xlabel('Time (us)')
ylabel('ACF')
hold off
drawnow
%% Fit to g1(tau) (for smoothness)
splineg1_Over = csape(tTemp,dataTemp,[2 2]);
splineMSD_Over = csape(tTemp,MSDTemp,[2 2]);
t_Over = logspace(log10(tTemp(1)),log10(tTemp(end)),1000)';
g1_Over = fnval(splineg1_Over,t_Over);
MSD_Over = fnval(splineMSD_Over,t_Over);
%%
% Parameters specification:
clear ILT_Input
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
ILT_Input.iwt = 4;        % User weights
ILT_Input.wt = dataTemp.^2;
ILT_Input.alpha_lims = [1e-10,1e-8];
ILT_Input.Nbg = 0;

[sg,gg,yfitg,lambdag,infog] = CRIE(tTemp,dataTemp,Ng,Nl,ILT_Input);
spg1 = csape(tTemp,yfitg,[2 2]);
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
ILT_Input.iwt = 4;        % User weights
ILT_Input.wt = dataTemp.^2;
ILT_Input.alpha_lims = [0.0001,0.01];
ILT_Input.Nbg = 2;

[s,g,yfit,lambda,info] = CRIE_MSD(tTemp,MSDTemp,Ng,Nl,ILT_Input);
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
ILT_Input.iwt = 4;        % User weights
ILT_Input.wt = g1_Over.^2;
ILT_Input.alpha_lims = [0.0001,0.01];
ILT_Input.Nbg = 2;

[sO,gO,yfitO,lambdaO,infoO] = CRIE_MSD(t_Over,MSD_Over,Ng,Nl,ILT_Input);
%%
MSDg1 = (6/q^2)*(-log(yfitg));
[nM,xM] = meshgrid((1:1)',yfitg);
MSDg1Log =(6/q^2)*(-sum(((-1).^(nM-1)).*((xM-1).^nM)./nM,2));
loglog(tTemp,MSDg1,tTemp,MSDg1Log)
%%
MSDmsd = yfit;
MSDmsdO = yfitO;
MSDg1=MSDg1/3; % Correction factor
MSDmsd = MSDmsd/3;
MSDmsdO = MSDmsdO/3;

splineMSD_Down = csape(t_Over,yfitO,[2 2]);

subplot(1,2,1)
loglog(tTemp,MSDTemp/3)
hold on
loglog(tTemp,MSDg1)
loglog(tTemp,MSDmsd)
loglog(t_Over,MSDmsdO)
hold off
subplot(1,2,2)
loglog(tTemp,(MSDTemp/3 - MSDg1).^2)
hold on
loglog(tTemp,(MSDTemp/3 - MSDmsd).^2)
loglog(tTemp,(MSDTemp/3 - fnval(splineMSD_Down,tTemp)/3).^2)
hold off
%%
s_Ref = abs(s);
s_Ref(s_Ref < 1e-6*max(s))=0; %Eliminate numerical noise in s (by clipping to 0)
s_Ref(end) = s(end);

Om = 1./tTemp;

[gM,OmM] = meshgrid(info.g,Om);

A_Om = 1./(1 + 1i*OmM./gM)*info.Coeff;
A_Om_Modes = (1./(1 + 1i*OmM./gM))*diag(info.Coeff);
y_Om = A_Om - 1i*s_Ref(end)./Om;
y_Om_Modes = A_Om_Modes - repmat(1i*s_Ref(end)./(Ng*Om),1,Ng);
Om = 1e6 * Om;

G_UR = kB*T./(3*pi*R*y_Om);
G_UR_Modes = kB*T./(3*pi*R*y_Om_Modes);
%%
subplot(1,2,1)
loglog(Om,real(y_Om_Modes),'k')
hold on
loglog(Om,real(y_Om),'r')
hold off

subplot(1,2,2)
loglog(Om,-imag(y_Om_Modes),'k')
hold on
loglog(Om,-imag(y_Om),'r')
hold off
%%
[OmMasg,GMasg]=MSDtoG_Mason(tTemp*1e-6,MSDg1,'R',115e-9,'CG',1.01,'T',T,'cutoff',0.03);
[OmMasgO,GMasgO]=MSDtoG_Mason(t_Over*1e-6,(2/q^2)*(-log(fnval(spg1,t_Over))),'R',115e-9,'CG',1.01,'T',T,'cutoff',0.03,'Gpoints',5000);
%%
loglog(OmMasg,real(GMasg),OmMasg,imag(GMasg))
%loglog(OMRF,real(GMRF),OMRF,imag(GMRF))
hold on
loglog(OmMasgO,real(GMasgO),OmMasgO,imag(GMasgO))
%loglog(OBRE,GBRE(:,1),OBRE,GBRE(:,2))
loglog(Om,9*real(G_UR),Om,9*imag(G_UR))
hold off
%% KV model
G_KV = 1;
tau_KV = 1e-4;
eta0_KV = 1e-2;
t_KV = logspace(-6,0);
CC_KV = (1/G_KV) * (1 - exp(-t_KV/tau_KV)) + t_KV/eta0_KV;
G_KV = 1./(1./(eta0_KV*1i./t_KV) + (1/G_KV) * (1./(1 + 1i*tau_KV./t_KV)));
subplot(1,2,1)
loglog(t_KV,CC_KV)
subplot(1,2,2)
loglog(1./t_KV,real(G_KV),1./t_KV,imag(G_KV))


