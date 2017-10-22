% Data preparation for microrheology with OM. Imports data from multiple OM
% files (csv files from AnalyseTrack).

% Clear everything
clear all; 
close all;
clc;

%Input: enter directory name
folder_name=uigetdir('/Users/AleCac/Desktop/');
listings=dir([folder_name,'/','*.csv']);
numFiles = length(listings);

%Input: enter um-to-pix conversion factor
um_to_pix = 0.1;

%Add imports to a cell array
dataTot = cell(1,numFiles);

parfor fileNum = 1:numFiles
    fileName = listings(fileNum).name;
    dataTot{fileNum} = OM_import_Analyse([folder_name,'/',fileName],1);
end

disp('Data imported. Post-processing...')
%%
%Check if there are files with non-optimal data (e.g. data with jumps, data
%not starting from t=1).
for k=1:length(dataTot)
    %Starting point
    first_t = dataTot{k}(1,1);
    if first_t ~= 1
        disp(['Warning: File ', num2str(k),' starts from ', num2str(first_t)]);
    end
    
    %Jumps
    jumps = diff(dataTot{k}(:,1));
    if any(jumps - 1)
        disp(['Warning: File ', num2str(k),' presents jumps.']);
    end
end
%%
%Get minimum array size (if there are different lengths)
%Resize the array
dataSize = min(cellfun('length',dataTot));
for k=1:length(dataTot)
    dataTot{k} = dataTot{k}(1:dataSize,:);
end
%%
%Convert to array & shuffle to account for indexing
dataTot = cell2mat(dataTot);
dataTot = um_to_pix*[dataTot(:,2:3:end),dataTot(:,3:3:end)];
%%
%Eliminate drift by fitting to a 0th-order polynomial
t = (1:size(dataTot,1))';
A = [ones(length(t),1)];%,t];%t.^2,t.^3];%,t,t.^2,t.^3];
bc = A\(dataTot-mean(dataTot));
dataTot = (dataTot-mean(dataTot)) - A*bc;

disp('Data processed. Saving...')
%%
%Save array in mat form in parent directory
[path,name]=fileparts(folder_name);
save([path,'/',name,'.mat'],'dataTot','-v6')


disp('Done!');


    


