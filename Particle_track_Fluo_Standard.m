%This routine executes a single or multi-particle tracking for bright/dark field
%images in which the particle has a defined blob appearance (for instance,
%FLUORESCENCE MICROSCOPY pictures or BRIGHT FIELD MICROSCOPY pictures on
%small particles). It applies a standard centroid location algorithm with
%sub-pixel resolution (see Track_Fin.m documentation for more details).

%INPUT: requires knowledge of the approximate DIAMETER of the blobs.
%REQUIRES (STRONGLY SUGGESTED) pretracking routines SignalToNoise.m and
%ImageFiltering.m to get a good value of the blob's diameter in order to
%maximize the signal-to-noise ratio.


close all
%Input: enter mother directory name (it must end with a "\" sign)
folder1='F:\2016-10-07 PS_stuck\';
%Input: enter child directory name
folder2='PS_PFSon_MicroscopeON_1000fps_LED0_LampON.07Oct2016_18.16.29';
%Input: enter filename stamp
stamp='PS_PFSon_MicroscopeON_1000fps_LED0_LampON.07Oct2016_18.16.29_';
%Input: enter the diameter of the blobs (should be the same as chosen
%during the pretracking routines
featdiam=5; 

featsize = floor(featdiam/2);
if mod(featsize,2) == 0
    featsize = featsize + 1;
end

listlen = length(dir([folder1,folder2,'\','*.tiff']));
%Input: by default, the program process all the pictures. Modify the 2nd
%and 3rd field if necessary
poslist=Ftrack_Fin([folder1,folder2,'\',stamp],1,listlen-1,featsize); 

%Input: by default, the program will create an output file in the mother
%directory. Modify if necessary
fullFileName=fullfile(folder1,['track_particle_output',folder2,'.csv']);
if exist(fullFileName,'file')==2
    %Input to specify how to deal with pre-existing file
    loop='A';
    while (loop ~= 'Y') && (loop ~= 'N')
        loop=input('Output track file already exists. Do you want to overwrite it? [Y/N] ','s');
        if loop=='Y'
            delete(fullFileName);
        elseif loop=='N'
            while exist(fullFileName,'file')==2
                fullFileName = [fullFileName(1:end-4),'bis.csv'];
            end
        else
            disp('Incorrect input from user.')
            loop = 'A';
        end
    end
end
    
csvwrite(fullFileName,poslist)

% param =struct('mem',0,'good',2,'dim',2,'quiet',0);
% pos=poslist(:,[2:3,1]);
% linked =track(pos,10,param);

% traj_length=tabulate(linked(:,4))
% msd=MSD(linked);
% loglog(msd(:,1),msd(:,2:3))
% fullFileName=fullfile(folder,'msd.csv');
% csvwrite(fullFileName,msd)