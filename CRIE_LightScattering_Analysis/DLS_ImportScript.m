clear variables

file = '\\np-nobelium\OE_Personal\ac2014\Documents\Malvern Instruments\Zetasizer\Export Data\exported.txt';
%Assumes a one-row header and five-column specifics
%Assumes a structure:
%   - [time,data,SNR]

exported=dlmread(file,'\t',1,5); 
t = exported(:,1:floor(length(exported)/2));
data = exported(:,floor(length(exported)/2) + 1:end-1);
SNR = exported(:,end);

export_folder = uigetdir;
export_file = fullfile(export_folder, 'DLS_data.mat');
save(export_file,'t','data','SNR')