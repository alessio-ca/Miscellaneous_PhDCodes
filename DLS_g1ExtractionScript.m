%% DLS_g1ExtractionScript

%Extract correlation curves from DLS data file
%Batch version
%Assumes that field names are:
%   - t
%   - data
%   - SNR

%Assumes decreasing temperatures

clear variables

NM = 5; %Number of measurements for each temperature
NT = 4; %Number of temperatures
Tin = 25;
Tstep = 5;

[FileName,PathName] = uigetfile;
load([PathName,FileName]);
if size(data,1) ~= NM*NT
    error('Wrong datafile length!');
end
%%
for i=1:NT
    for j=1:NM
        [tau,MSD,ACF]=DLS_Analysis(t((i-1)*NM + j,:)',data((i-1)*NM +j,:)','T',273 + Tin - (i-1)*Tstep,'Nrun',6);
        save([PathName,'2_T_',num2str(Tin - (i-1)*Tstep),'_',num2str(j)],'tau','ACF','MSD');
    end
end
