t = t';
data = data';

%% Constants
kB = 1.38*10^(-23);
lambda = 633*10^(-9);
theta = 173*2*pi/360;
q = 4*pi*n*sin(theta/2)/lambda;

T = 298 + 5;
n = 1.33;
R = 115e-9;
Beta = 1;
cut = 0.1;
tail = [6e5,3e6];
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
%%
MSDg1 = (6/q^2)*(-log(yfitg));
MSDmsd = yfit;
MSDg1=MSDg1/3; % Correction factor
MSDmsd = MSDmsd/3;
tTemp = tTemp; % Correction factor
subplot(1,2,1)
loglog(tTemp,MSDTemp/3)
hold on
loglog(tTemp,MSDg1)
loglog(tTemp,MSDmsd)
hold off
subplot(1,2,2)
loglog(tTemp,(MSDTemp/3 - MSDg1).^2)
hold on
loglog(tTemp,(MSDTemp/3 - MSDmsd).^2)
hold off
%%

s_Ref = abs(s);
s_Ref(s_Ref < 1e-6*max(s))=0; %Eliminate numerical noise in s (by clipping to 0)
s_Ref(end) = s(end);

Om = 1./tTemp;

[gM,OmM] = meshgrid(info.g,Om);

A_Om = 1./(1 + 1i*OmM./(10.^gM))*info.Coeff;
y_Om = A_Om - 1i*s_Ref(end)./Om;
Om = 1e6 * Om;

G_UR = kB*T./(3*pi*R*y_Om);
[OmMasg,GMasg]=MSDtoG_Mason(tTemp*1e-6,MSDg1,'R',115e-9,'CG',1.01,'T',T,'cutoff',0.03);
%[OmMas,GMas]=MSDtoG_Mason(tTemp*1e-6,yfit,'R',115e-9,'CG',1.01,'T',T,'cutoff',0.03);
%%
loglog(OmMasg,real(GMasg),OmMasg,imag(GMasg))
%loglog(OMRF,real(GMRF),OMRF,imag(GMRF))
hold on
%loglog(OBRE,GBRE(:,1),OBRE,GBRE(:,2))
loglog(Om,9*real(G_UR),Om,9*imag(G_UR))
hold off


