%% DLS_g1Average

%Script to average the ACFs from the data from DLS_Analysis if multiple
%files are present.

%Input data are variables with names
% - tau (output from DLS_Analysis)
% - ACF (output from DLS_Analysis)
% - SNR (raw from acquisition data)

%Assumes decreasing temperatures with Trend = -1
%Assumes increasing temperatures with Trend = +1


clear variables
close all

NM = 3; %Number of measurements for each temperature
NT = 2; %Number of temperatures
Tin = 20;
Tstep = 2;

Trend = +1;

[FileName,PathName] = uigetfile;
load([PathName,FileName]);

if size(data,1) ~= NM*NT
    error('Wrong datafile length!');
end
C = cell(NT,NM);
%% Matrix building
for i=1:NT
    for j=1:NM
        load([PathName,'T_',num2str(Tin - (i-1)*Tstep),'_',num2str(j)]);
        C{i,j} = [tau,ACF];
    end
end
%% Extraction of min and max times
minT = cellfun(@(x) min(x(:,1)),C);
minT = max(minT,[],2);
maxT = cellfun(@(x) max(x(:,1)),C);
maxT = min(maxT,[],2);

tFit = struct;
for i = 1:NT
    for j=1:NM
        tFit(i,j).time = logspace(log10(minT(i)),log10(maxT(i)),100)';
    end
end
%% Data spline interpolations
data_spline=cellfun(@(x) csape(x(:,1)',x(:,2)'),C);
ACF_vec=arrayfun(@(x,y) fnval(x.time,y),tFit,data_spline,'UniformOutput',false);
ACF_vec = permute(reshape([ACF_vec{:}], [], size(ACF_vec, 1), size(ACF_vec, 2)), [1 2 3]);
SNR_mat = reshape(SNR,[NM,NT])';
SNR_mat = permute(repmat(SNR_mat,1,1,100),[3,1,2]);
ACF_avg = sum(ACF_vec.*SNR_mat,3)./sum(SNR_mat,3);
%% Fitting the weighted average ACF
for i = 1:NT
    [tau,MSD,ACF]=DLS_Analysis(tFit(i,1).time,ACF_avg(:,i).^2,'T',273 + Tin - (i-1)*Tstep,'Nrun',6,'beta',0);
    save([PathName,'Avg_T_',num2str(Tin - (i-1)*Tstep)],'tau','ACF','MSD');
end

