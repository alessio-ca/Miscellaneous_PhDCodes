%% DLS_g1ExtractionScript_PreRun

%PreRun to fine-tune parameters before extraction of correlation curves 
%from DLS data file

%Assumes that field names are:
%   - t
%   - data
%   - SNR

%Assumes decreasing temperatures with Trend = -1
%Assumes increasing temperatures with Trend = +1


clear variables

NM = 3; %Number of measurements for each temperature
NT = 8; %Number of temperatures
Tin = 60;
Tstep = 5;

Trend = -1;

[FileName,PathName] = uigetfile;
load([PathName,FileName]);
if size(data,1) ~= NM*NT
    error('Wrong datafile length!');
end

i = 1;
j = 1;
[tau,MSD,ACF]=DLS_CRIEAnalysis_Simple(t((i-1)*NM + j,:)',data((i-1)*NM +j,:)','T',273 + Tin + Trend*(i-1)*Tstep,'cut',0.1);
