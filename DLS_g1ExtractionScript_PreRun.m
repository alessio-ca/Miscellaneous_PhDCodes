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

NM = 2; %Number of measurements for each temperature
NT = 7; %Number of temperatures
Tin = 29;
Tstep = 3;

Trend = -1;

[FileName,PathName] = uigetfile;
load([PathName,FileName]);
if size(data,1) ~= NM*NT
    error('Wrong datafile length!');
end

i = 1;
j = 1;
[tau,MSD,ACF]=DLS_CRIEAnalysis(t((i-1)*NM + j,:)',data((i-1)*NM +j,:)','T',273 + Tin + Trend*(i-1)*Tstep,'batch',0,'cut',0.1);
