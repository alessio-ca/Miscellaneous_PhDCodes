%% DLS_g1ExtractionScript

%Extract correlation curves from DLS data file
%Batch version
%Assumes that field names are:
%   - t
%   - data
%   - SNR

%Assumes decreasing temperatures with Trend = -1
%Assumes increasing temperatures with Trend = +1


clear variables

NM = 1; %Number of measurements for each temperature
NT = 8; %Number of temperatures
Tin = 60;
Tstep = 5;

Trend = -1;

[FileName,PathName] = uigetfile;
load([PathName,FileName]);
if size(data,1) ~= NM*NT
    error('Wrong datafile length!');
end
%%
for i=1:NT
    for j=1:NM
        %[tau,MSD,ACF]=DLS_Analysis(t((i-1)*NM + j,:)',data((i-1)*NM +j,:)','T',273 + Tin - (i-1)*Tstep,'Nrun',6);
        [tau,MSD,ACF,tau_Long,ACF_Long]=DLS_CRIEAnalysis_Simple(t((i-1)*NM + j,:)',data((i-1)*NM +j,:)','T',273 + Tin + Trend*(i-1)*Tstep,'cut',0.05);
        save([PathName,'T_',num2str(Tin + Trend*(i-1)*Tstep),'_',num2str(j)],'tau','ACF','tau_Long','ACF_Long','MSD');
    end
end
