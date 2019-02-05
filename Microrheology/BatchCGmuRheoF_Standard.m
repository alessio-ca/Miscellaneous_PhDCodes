%This routine execute the standard coarse-graining microrheology



close all;
%Input: enter directory for data analysis (remember to put \ at the end!)
dirIn='/Users/AleCac/Documents/Python_stuff/';
%Input: enter name(s) of file to be analyzed. Wildcard (*) accepted.
dirSel = 'track_particle_output60x_1.20_Si02_1000FPS_1.5um_Tweezer_G_0.05.13Aug2016_15.27.31*';
%Input: enter name of directory where out files will be saved
Out = 'CG_Microrheo_Analysis';


dirOut=[dirIn,Out];
if ~exist(dirOut, 'dir')
    % Folder does not exist so create it.
    mkdir(dirOut);
end
dirList = dir([dirIn,dirSel]);
if isempty(dirList)
    % No files to be analyzed
    error('There are no input files, fool!');
end
for i = 1:size(dirList,1)
    [~,Outfilename,~] = fileparts(dirList(i).name);
    try
        %Inputs are (N,F,M,scale,Dim,ShutT,srate,Temp,R,SwitchFifo)
        %N - number of points per array
        %F - factor by which arrays are coarse-grained
        %M - number of arrays
        %scale - factor to convert displacement series to physical units
        %Dim - number of columns of data after index and time
        %   Data should be formatted
        %ShutT  - Shutter exposure time for camera, if visualising
        %       - can be set to 0
        %srate - Sampling rate
        %Temp - Temperature
        %R - radius of probe
        %SwitchFifo - Set to 1 if reading from a named pipe, 0 otherwise
        %ParticleN - index of particle to be studied.
        %filename - File to be analyzed
        %header - Set to 1 if your file has a one-line header (more is NOT allowed)
        %out - Output filename
        %Counter - optional for batch analysis
        CGMuRheoF_CSV(100,4,5,100e-9,2,0,1000,298,1.5e-6 / 2,0,1,strcat(dirIn,dirList(i).name),0,strcat(dirOut,Outfilename),i)
    catch ME
        if (strcmp(ME.identifier,'MATLAB:badsubscript'))
            sprintf('Track %d Done',i)
        else
            rethrow(ME)
        end
    end
end

disp('Moduli calculated');
beep;