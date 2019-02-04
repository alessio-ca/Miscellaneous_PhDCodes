t = t';
data = data';

%% Constants
n = 1.33;
R = 115e-9;
Beta = 1;
tail = [0,0];

kB = 1.38*10^(-23);
lambda = 633*10^(-9);
theta = 173*2*pi/360;
q = 4*pi*n*sin(theta/2)/lambda;
%% g1(tau) calc. 
for i=1:size(data,2)
    
    %Delete first 6 points (usually noise) & format data
    Time_ToFit = t(6:end,i);
    ACF_ToFit = data(6:end,i);
    if all(tail)
        tailFit = ACF_ToFit(Time_ToFit > tail(1) & Time_ToFit < tail(2));
        tailFit = mean(tailFit);
        ACF_ToFit = (ACF_ToFit - tailFit)/(1+tailFit);
    end
    cut = 0.1*max(ACF_ToFit);
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
    save(['ACF_',num2str(i),'.mat'],'tTemp','dataTemp')
    save(['MSD_',num2str(i),'.mat'],'tTemp','MSDTemp')
end

%%
import_folder = uigetdir;
%%
initialColorOrder = get(gca,'ColorOrder');
newDefaultColors = flip(cool(11));
set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');
hold on
for i = 1:size(data,2)
    load([import_folder,'/','ACF_',num2str(i)]);
    plot(tTemp,dataTemp);
end
hold off
h=gca;
h.XScale='log';
xlabel('Time [s]')
ylabel('g_1 (\tau)')
%%
initialColorOrder = get(gca,'ColorOrder');
newDefaultColors = flip(cool(9));
set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');
hold on
for i = 1:size(data,2)
    load([import_folder,'/','MSD_',num2str(i)]);
    plot(tTemp,MSDTemp);
end
hold off
h=gca;
h.XScale='log';
h.YScale='log';
xlabel('Time [s]')
ylabel('MSD')

